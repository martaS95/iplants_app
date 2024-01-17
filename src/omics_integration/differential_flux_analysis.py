import os
import scipy as sp
from scipy.stats import sem, hypergeom
import statsmodels.api
from cobra.sampling import ACHRSampler
from cobra.io import read_sbml_model
from utils.config import PROJECT_PATH
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class DFA:
    """
    Class to perform differential flux analysis between metabolic models based on the paper from
    ...
    """
    def __init__(self, modelid, datasetid, specific_models: dict, models_objective: dict, pathways_map: str = None):
        """
        Parameters
        ----------
        modelid: str
            model used to reconstruct the specific models
        datasetid: str
            datasets used to reconstruct the specific models
        specific_models: dict
            name of the model for each tissue
        models_objective: dict
            objective reactions to maximize in simulations
        pathways_map: str, Optional
            path of the csv file containing the pathway of each reaction
        """

        self.models_folder = os.path.join(PROJECT_PATH, 'reconstruction_results', modelid, 'results_troppo', datasetid,
                                          'reconstructed_models')

        self.results_folder = os.path.join(PROJECT_PATH, 'reconstruction_results', modelid, 'results_troppo', datasetid,
                                           'dfa')

        self.specific_models = specific_models
        self.objectives = models_objective
        self.sampled_fluxes = None
        self.pathways = pathways_map
        self.results = None

    def sampling(self, thinning: int = 100, n_jobs: int = 4):
        """
        Performs flux sampling
        Parameters
        -------
        thinning: int (default = 100)
            the thinning factor of the generated sampling chain. A thinning of
            10 means samples are returned every 10 steps
        n_jobs: int (default = 4)
            number of threads to use for flux sampling
        Returns
        -------
        self.sampled_fluxes: dict
            sampled fluxes created (pandas dataframe) as values
        """

        sampling_dic = {}

        for tissue in self.specific_models:
            modelname = self.specific_models[tissue]
            print(modelname, 'starting flux sampling')
            try:
                df_sampling = pd.read_csv(os.path.join(self.results_folder, '%s_sampling.csv' % modelname),
                                          index_col=0)
            except FileNotFoundError:
                model_path = os.path.join(self.models_folder, '%s.xml' % modelname)
                model_obj = read_sbml_model(model_path)
                model_obj.objective = self.objectives[tissue]

                df_sampling = ACHRSampler(model_obj, thinning=thinning, n_jobs=n_jobs)

            sampling_dic[modelname] = df_sampling

        self.sampled_fluxes = sampling_dic
        return self.sampled_fluxes

    def kstest(self):
        """
        Calculate the K-S test to detect significantly altered reactions fluxes.
        Results are saved in a csv file.

        Returns
        -------
        list: The results of the K-S test for each reaction.

        """

        modelnames = '_'. join(list(self.sampled_fluxes.keys()))
        sampled_fluxes1, sampled_fluxes2 = list(self.sampled_fluxes.values())
        rxns1 = set(sampled_fluxes1.columns)
        rxns2 = set(sampled_fluxes2.columns)

        rxns1_unique = rxns1.difference(rxns2)
        rxns2_unique = rxns2.difference(rxns1)

        for reac in rxns1_unique:
            sampled_fluxes2[reac] = [0.0000000] * sampled_fluxes2.shape[0]

        for reac in rxns2_unique:
            sampled_fluxes1[reac] = [0.0000000] * sampled_fluxes2.shape[0]

        rxns_common = set(sampled_fluxes1.columns)

        pvals = []
        rxnid = []
        fc = []

        for rxn in rxns_common:
            data1 = sampled_fluxes1[rxn].round(decimals=4)
            data2 = sampled_fluxes2[rxn].round(decimals=4)

            data1 = data1.sample(n=1000)
            data2 = data2.sample(n=1000)

            if (data1.std() != 0 and data1.mean() != 0) or (data2.std() != 0 and data2.mean() != 0):
                kstat, pval = sp.stats.ks_2samp(data1, data2)

                foldc = (data1.mean() - data2.mean()) / abs(data1.mean() + data2.mean())

                pvals.append(pval)
                rxnid.append(rxn)
                fc.append(foldc)

        data_mwu = pd.DataFrame({'Reaction': rxnid, 'Pvalue': pvals})
        data_mwu = data_mwu.set_index('Reaction')

        reject, padj, _, _ = statsmodels.stats.multitest.multipletests(data_mwu['Pvalue'], alpha=0.05, method='fdr_bh',
                                                                       is_sorted=False, returnsorted=False)

        data_mwu['Padj'] = padj
        data_mwu['Reject'] = reject
        data_mwu['FC'] = fc

        data_sigFC = data_mwu.loc[(abs(data_mwu['FC']) > 0.82) & (data_mwu['Padj'] < 0.05), :]

        rxns1 = set(sampled_fluxes1.columns)
        rxns2 = set(sampled_fluxes2.columns)

        rxn_in1 = rxns1.difference(rxns2)
        rxn_in2 = rxns2.difference(rxns1)

        act = []
        rep = []

        for rx in rxn_in1:  # Activated reactions
            sig = bootstrapCI(sampled_fluxes1[rx])

            if sig == 1:
                act.append(rx)

        for rx in rxn_in2:  # Activated reactions
            sig = bootstrapCI(sampled_fluxes2[rx])

            if sig == 1:
                rep.append(rx)

        df_abs = pd.DataFrame({'Reaction': act + rep, 'Padj': np.zeros(len(act + rep))})
        df_abs = df_abs.set_index('Reaction')
        data_return = data_sigFC + df_abs

        file = os.path.join(self.results_folder, '%s_DFA_reaction_result_all.csv' % modelnames)
        data_sigFC.to_csv(file)

        self.results = data_sigFC.index.to_list()

        return self.results

    def pathway_enrichment(self):
        """
        Maps significantly altered reactions to pathways using the pathways self.pathways.
        Results are saved in csv and jpg files.
        """
        if not self.pathways:
            print('No pathway information is available, not possible to perform pathway enrichment analysis')
            return None

        subs = pd.read_csv(self.pathways, dtype=str)

        dataset = pd.DataFrame()

        for path in subs.columns:
            reaction_set = subs[path]

            rxn = reaction_set.reset_index(drop=True)

            df_temp = pd.DataFrame({path: rxn})

            dataset = pd.concat([dataset, df_temp], axis=1)

        listrxnSize = []
        setSize = []

        d = [g for g in self.results]

        for col in dataset.columns:
            df = pd.DataFrame({'Reaction': dataset[col]})
            df.dropna()

            out = []

            for reac in df['Reaction']:
                if reac in self.results:
                    out.append(reac)
                    if reac in d:
                        d.remove(reac)

            listrxnSize.append(len(out))
            setSize.append(len(dataset[col].dropna()))

        hyperdata = pd.DataFrame({'Pathways': dataset.columns, 'ListReactions': listrxnSize, 'SetSize': setSize})

        hits = hyperdata['ListReactions']
        pool = hyperdata['SetSize']

        allrxns = hyperdata['SetSize'].sum()
        targetrxns = hyperdata['ListReactions'].sum()

        pvalList = []

        for h, p in zip(hits, pool):
            rv = hypergeom(allrxns - p, p, targetrxns)

            pval = rv.pmf(h)

            pvalList.append(pval)

        hyperdata['P-value'] = pvalList

        reject, padj, _, _ = statsmodels.stats.multitest.multipletests(hyperdata['P-value'], alpha=0.05,
                                                                       method='fdr_bh',
                                                                       is_sorted=False, returnsorted=False)

        hyperdata['P-value_adj'] = padj
        hyperdata['Reject'] = reject

        hyperdata_sig = hyperdata[(hyperdata['Reject']) & (hyperdata['ListReactions'] != 0)]

        hyperdata_sorted = hyperdata_sig.sort_values(by='P-value_adj', ascending=False)
        hyperdata_sorted = hyperdata_sorted[hyperdata_sorted['Pathways'] != 'transporters']
        hyperdata_sorted = hyperdata_sorted[hyperdata_sorted['Pathways'] != 'drains']

        plt.figure(figsize=(12, 10))

        sc = plt.scatter(hyperdata_sorted['P-value_adj'], np.arange(0, len(hyperdata_sorted['Pathways'])),
                         s=hyperdata_sorted['ListReactions'], color=(0.9, 0.3, 0.1, 0.9))

        plt.xlabel('Adjusted p-value')

        plt.yticks(np.arange(0, len(hyperdata_sorted['Pathways'])), labels=hyperdata_sorted['Pathways'])

        try:
            handles, labels = sc.legend_elements(prop="sizes", alpha=0.8)
            plt.legend(handles, labels, bbox_to_anchor=(1.6, 1.02), loc='upper right', title="Reactions")

            plt.tight_layout()

            names = '_'.join(list(self.specific_models.keys()))

            plt.savefig(os.path.join(self.results_folder, '%s_DFA_pathway_result_all.png' % names), dpi=600)

            hyperdata_sorted.to_csv(os.path.join(self.results_folder, '%s_DFA_pathway_result_all.csv' % names))
        except ValueError:
            return None

    def run_complete_dfa(self, thinning: int = 100, n_jobs: int = 4):
        """
        Run the three steps of dfa: sampling, ks test and pathway enrichment
        Parameters
        -------
        thinning: int (default = 100)
            the thinning factor of the generated sampling chain. A thinning of
            10 means samples are returned every 10 steps
        n_jobs: int (default = 4)
            number of threads to use for flux sampling
        """
        print('Sampling is starting...')
        self.sampling(thinning=thinning, n_jobs=n_jobs)
        print('Sampling is complete!')

        print('Starting KS test')
        self.kstest()
        print('KS is complete!')

        print('Starting Pathway enrichment analysis test')
        self.pathway_enrichment()
        print('Pathway enrichment analysis is complete!')

        print('Differential flux analysis is complete!')


def bootstrapCI(rxn: pd.Series):
    """
    Calculate the confidence interval of a reaction

    Parameters
    ----------
    rxn: pd.Series
        values for the reaction

    Returns
    -------
    int: 1 if the reaction is significantly different from the mean, 0 otherwise.

    """
    bsci = []

    for i in range(1000):
        bt_samp = rxn.sample(1000, replace=True)
        bsci.append(bt_samp.mean())

    ci_low = np.percentile(bsci, 2.5)
    ci_high = np.percentile(bsci, 97.5)

    if ci_low > 0 or ci_high < 0:
        return 1
    else:
        return 0


def get_sampling_fluxes():
    folder = os.path.join(PROJECT_PATH, 'reconstruction_results', 'vvinif2023', 'results_troppo', 'ALL_BERRY',
                          'dfa')

    files = [f for f in os.listdir(folder) if 'sampling.csv' in f]

    all_fluxes = []

    for f in files:
        df = pd.read_csv(os.path.join(folder, f), index_col=0)
        subdf = df
        subdf.index = [f + '_' + str(x) for x in df.index]
        all_fluxes.append(subdf)

    all_df = pd.concat(all_fluxes)
    all_df = all_df.fillna(0)
    all_df.to_csv(os.path.join(folder, 'all_sampling.csv'))


def complete_reactions_dfa(base_model_file, dfa_reacs_dir):

    base_model = read_sbml_model(base_model_file)

    files = [f for f in os.listdir(dfa_reacs_dir) if f.endswith('reaction_result_all.csv')]
    for f in files:
        df = pd.read_csv(os.path.join(dfa_reacs_dir, f))
        reactions = df['Reaction'].to_list()
        paths_by_reactions = []
        for reac in reactions:
            reac_obj = base_model.reactions.get_by_id(reac)
            groups = base_model.get_associated_groups(reac_obj)
            groups_list = []
            for g in groups:
                groups_list.append((g.id, g.name))
            paths_by_reactions.append(groups_list)
        df['pathways'] = paths_by_reactions
        df.to_csv(os.path.join(dfa_reacs_dir, f.replace('.csv', '_paths.csv')), index=False)


if __name__ == '__main__':
    model_id = 'vvinif2023'
    dataset_id = 'ALL_BERRY'

    basefolder = os.path.join(PROJECT_PATH, 'reconstruction_results', model_id)
    base_file = os.path.join(basefolder, 'vvinif2023_FINAL.xml')

    dfa_dir = os.path.join(PROJECT_PATH, 'reconstruction_results', model_id, 'results_troppo', dataset_id,
                           'dfa')

    paths = os.path.join(PROJECT_PATH, 'reconstruction_results', model_id, 'pathways.csv')

    models_folder = os.path.join(PROJECT_PATH, 'reconstruction_results', model_id, 'results_troppo', dataset_id,
                                 'reconstructed_models')

    model_files = os.listdir(models_folder)
    model_files = [f[:-4] for f in model_files]

    models_dic = {x: x for x in model_files}

    obj = {x: 'e-Biomass_vvinif2023__cyto' for x in model_files}

    for model in model_files:

        obj_cls = DFA(modelid=model_id, datasetid=dataset_id, specific_models=models_dic, models_objective=obj,
                      pathways_map=paths)
        obj_cls.sampling()

    # complete_reactions_dfa(base_file, dfa_dir)

    # get_dfa_reactions()
    # get_sampling_fluxes()
