import os
import pandas as pd
import numpy as np
from utils.config import OMICS_DATADIR
from cobra.io import read_sbml_model


class OmicsDataset:
    """
    Class to define and preprocess the omics_integration dataset to use for integration with models
    Some methods of this class are specific for the dataset under analysis (GSE36128) - mapping probes to gene ids and
    different annotation versions
    Others are specific for the same dataset but extract from GREAT database
    """

    def __init__(self, dataset_id: str, samples: list = None, genes: list = None):
        """
        Initializes an Omics Dataset instance
        Parameters
        ----------
        dataset_id: str
            id of the dataset to use. Folder and files will be named with this id
        samples: list, optional (default = None)
            list of samples in the dataset
        genes: list, optional (default = None)
            list of genes in the dataset
        """

        self.dataset_id = dataset_id
        self.samples = samples
        self.genes = genes
        self.original_data = None
        self._data = None

        if not os.path.exists(os.path.join(OMICS_DATADIR, dataset_id, 'UPDATE')):
            os.makedirs(os.path.join(OMICS_DATADIR, dataset_id))

        self.data_file = os.path.join(OMICS_DATADIR, dataset_id, 'UPDATE', dataset_id + '.csv')

        # if the file exists, loads it, otherwise preprocessing methods should be run
        # to get the final omics_integration dataset

        if os.path.isfile(self.data_file):
            self.data = pd.read_csv(self.data_file, index_col=0)

            if not self.samples:
                self.samples = self.data.columns
            if not self.genes:
                self.genes = self.data.index

    @property
    def data(self) -> pd.DataFrame:
        return self._data

    @data.setter
    def data(self, value):
        self._data = value

    def load_data_files(self) -> pd.DataFrame:
        """
        Loads data from tsv files for each sample with normalized data from GEO and writes a CSV file
        files with the folling structure:
        'ID_REF'    'VALUE'
        ...         ...
        :return:
        self.data: pd.DataFrame
            expression data for all samples
        """
        folder_data = os.path.join(OMICS_DATADIR, self.dataset_id, 'expr_files')

        files = os.listdir(folder_data)
        dic_info = {}
        rows = []
        for f in files:
            path = os.path.join(folder_data, f)
            table = pd.read_csv(path, sep='\t', index_col=0)
            if not rows:
                rows = list(table.index)
            dic_info[f[:-4]] = table['VALUE']

        self.original_data = pd.DataFrame(dic_info, index=rows)

        original_file = os.path.join(OMICS_DATADIR, self.dataset_id, self.dataset_id + '_original.csv')
        self.original_data.to_csv(original_file)

        # remove genes with nan values

        self.data = self.original_data.dropna()
        self.data.to_csv(self.data_file)

        return self.original_data

    def replace_probe_ids_by_genes(self, mapping_file: str) -> pd.DataFrame:
        """
        Reads a pair file with microarray data that has the genes corresponding to each probe and
        gets the mean expression value for each gene (there are four probes for each gene)
        Parameters
        ----------
        mapping_file: str
            file with the probes and the corresponding gene

        Returns
        -------
        self.data: pd.Dataframe
            mean expression value for each gene (index is genes_id and not probes ids)
        """
        full_path = os.path.join(OMICS_DATADIR, self.dataset_id, mapping_file)
        df = pd.read_csv(full_path, header=1, sep='\t')
        df = df.loc[:, ["SEQ_ID", "PROBE_ID"]]
        df = df[df['SEQ_ID'] != 'RANDOM_GC48_DUAL']

        gene_probe_map = {}

        for ind, row in df.iterrows():
            gene = row['SEQ_ID']
            probe = row['PROBE_ID']

            gene_probe_map[probe] = gene

            if gene not in gene_probe_map:
                gene_probe_map[gene] = [probe]
            else:
                gene_probe_map[gene].append(probe)

        new_df = self.data.copy()
        new_df = new_df.filter(items=[x for x in new_df.index if not x.startswith('RANDOM')], axis=0)
        new_df['genes_id'] = new_df.index.map(gene_probe_map)

        new_df = new_df.groupby(['genes_id']).mean()

        new_df.to_csv(self.data_file)

        self.data = new_df
        self.genes = list(self.data.index)
        self.samples = list(self.data.columns)

        return self.data

    def replace_genes_by_v1ids(self, mapping_file: str) -> pd.DataFrame:
        """
        Get the consensus id for each gene using the ids in the raw data.
        Parameters
        ----------
        mapping_file: str
            name of the file with the mapping of gene ids (VitisNet_manual_annotation.tsv)

        Returns
        -------
        self.data: pd.Dataframe
            dataset with v1 ids
        """
        old_names = self.data.index
        df_map = pd.read_excel(os.path.join(OMICS_DATADIR, self.dataset_id, mapping_file), index_col=0)
        df_map = df_map.loc[:, ['dataset annotation', 'reference annotation']]
        df_map = df_map.drop_duplicates()

        map_dic = {}

        for name in old_names:
            match = df_map['reference annotation'][df_map['dataset annotation'] == name].item()
            map_dic[name] = match

        new_index = []

        for ind in self.data.index:
            new_index.append(map_dic[ind])

        self.data.index = new_index
        self.genes = list(self.data.index)

        self.data.to_csv(self.data_file)

        return self.data

    def replace_v1_v4(self, mapping_file: str) -> pd.DataFrame:
        """
        Get the 'v4 annotation' for each gene using the consensus gene annotation. Some data will be lost.
        Parameters
        ----------
        mapping_file: str
            xlsx file mapping the different annotation versions

        Returns
        -------
        self.data: pd.Dataframe
            dataset with v4 ids as index
        """
        df_map = pd.read_csv(os.path.join(OMICS_DATADIR, self.dataset_id, mapping_file),
                             index_col=0)
        df_map = df_map['v1']
        df_map = df_map.dropna()

        dic = {}

        for ind in df_map.index:
            if '|' in df_map[ind]:
                v1_ids = df_map[ind].split(' | ')
                for ide in v1_ids:
                    dic[ide] = ind
            else:
                dic[df_map[ind]] = ind

        inds = []
        rows = []
        for old_ind, row in self.data.iterrows():
            if '<br>' not in old_ind:
                try:
                    new_ind = dic[old_ind]
                    inds.append(new_ind)
                    rows.append(row)
                except KeyError:
                    pass
            else:
                old_inds = str(old_ind).split('<br>')
                for i in old_inds:
                    try:
                        new_ind = dic[i]
                        inds.append(new_ind)
                        rows.append(row)
                    except KeyError:
                        pass

        new_df = pd.DataFrame(rows, columns=self.data.columns)
        new_df.index = inds

        new_df.to_csv(self.data_file)
        self.data = new_df

        return self.data

    def replicate_same_ids(self, proteins_file: str = None) -> pd.DataFrame:
        """
        Get the v4 annotation for all the proteins for each gene
        Parameters
        ----------
        proteins_file: str, optional
            faa file with the protein sequences in the genome

        Returns
        -------
        self.data: pd.Dataframe
            dataset with v4 ids as index
        """

        model = read_sbml_model('C:/Users/BiSBII/Documents/plantdb/reconstruction_results/vvinif2023/'
                                'vvinif2023_FINAL.xml')

        dic_genes = {}

        for gene in model.genes:
            new_gene_id, prot_id = gene.id.split('_')
            if new_gene_id not in dic_genes:
                dic_genes[new_gene_id] = [prot_id]
            else:
                dic_genes[new_gene_id].append(prot_id)

        new_rows = []
        new_index = []

        for ind, row in self.data.iterrows():
            replicates = dic_genes[ind]

            for rep in replicates:
                new_ind = str(ind) + '_' + rep
                new_index.append(new_ind)
                new_rows.append(row)

        new_df = pd.DataFrame(new_rows, columns=self.data.columns)
        new_df.index = new_index

        new_df.to_csv(self.data_file)

        self.data = new_df

        return self.data

    def replace_v1ids_locus(self, mapping_file: str) -> pd.DataFrame:
        """
        Get the 'locus tag' for each gene using the consensus gene annotation. Some data will be lost.
        Parameters
        ----------
        mapping_file: str
            xlsx file mapping the different annotation versions

        Returns
        -------
        self.data: pd.Dataframe
            dataset with locustag ids as index
        """
        df_map = pd.read_excel(os.path.join(OMICS_DATADIR, self.dataset_id, mapping_file),
                               sheet_name='in_refseq', index_col=0)
        df_map = df_map.loc[:, ['reference annotation', 'locus tag']]

        new_rows = []
        new_index = []

        for ind, row in self.data.iterrows():
            if '<br>' not in ind:
                try:
                    locus = df_map['locus tag'][df_map['reference annotation'] == ind].to_list()
                    for x in locus:
                        new_row = row.tolist()
                        new_rows.append(new_row)
                        new_index.append(x)
                except ValueError:
                    pass
            else:
                inds = str(ind).split('<br>')
                for i in inds:
                    try:
                        locus = df_map['locus tag'][df_map['reference annotation'] == i].to_list()
                        for x in locus:
                            new_row = row.tolist()
                            new_rows.append(new_row)
                            new_index.append(x)
                    except ValueError:
                        pass

        new_df = pd.DataFrame(new_rows, columns=self.data.columns)
        new_df['new_locus'] = new_index

        new_df = new_df.groupby(['new_locus']).mean()

        new_df.to_csv(self.data_file)

        self.data = new_df
        self.genes = list(self.data.index)
        return self.data

    def replace_locusid_proteins(self, mapping_file: str) -> pd.DataFrame:
        """
        Get the 'protein id' for each gene from the locus tag (protein id is the annotation in the model)
        Parameters
        ----------
        mapping_file: str
            csv file mapping genes and proteins for vitis vinifera genome (mapping_locus_proteinid.csv)

        Returns
        -------
        self.data: pd.Dataframe
            omics_integration dataset with protein ids as index
        """
        new_rows = []
        new_index = []

        df_map = pd.read_csv(os.path.join(OMICS_DATADIR, self.dataset_id, mapping_file), index_col=0)

        for ind, row in self.data.iterrows():
            try:
                proteins = df_map['protein'][df_map['gene'] == ind].to_list()
                for p in proteins:
                    new_row = row.tolist()
                    new_rows.append(new_row)
                    new_index.append(p)
            except ValueError:
                pass

        new_df = pd.DataFrame(new_rows, index=new_index, columns=self.data.columns)
        new_df.to_csv(self.data_file)

        self.data = new_df
        self.genes = list(self.data.index)
        return self.data

    def remove_replicates(self) -> pd.DataFrame:
        """
        Removes the replicates for each sample by calculating the mean.
        Returns
        -------
        self.data: pd.Dataframe
           mean values of expression data for each sample without replicates
        """
        self.samples = ['leaf', 'stem', 'berry_green', 'berry_mature']

        self.data = pd.DataFrame(np.log2(self.data.values + 1.1), index=self.data.index, columns=self.data.columns)

        metadata = pd.read_excel(os.path.join(OMICS_DATADIR, self.dataset_id, self.dataset_id + '_metadata.xlsx'))
        tissues = set(metadata['tissue'])
        samples = []
        col_names = []
        for tissue in tissues:
            sample_ids = metadata['ENA'][metadata['tissue'] == tissue].to_list()
            descp = metadata['description'][metadata['tissue'] == tissue].to_list()
            year = metadata['year'][metadata['tissue'] == tissue].to_list()
            i = 0
            while i < len(sample_ids):
                chunk = sample_ids[i:i + 3]
                samples.append(chunk)
                new_name = tissue + '_' + descp[i][:-5] + '_' + str(year[i])
                col_names.append(new_name)
                i += 3

        mean_rep_cols = []
        for cols in samples:
            mean_col = self.data[cols].mean(axis=1)
            mean_rep_cols.append(mean_col)

        df_mean_noreps = pd.concat(mean_rep_cols, axis=1)
        df_mean_noreps.columns = col_names
        df_mean_noreps.index = self.data.index

        df_mean_noreps.to_csv(os.path.join(OMICS_DATADIR, self.dataset_id, 'UPDATE',
                                           self.dataset_id + 'GREAT_LOGTPM100_noreps.csv'))

        join_samples = []
        for sample in self.samples:
            cols_sample = [x for x in df_mean_noreps.columns if x.startswith(sample)]
            mean_col = df_mean_noreps[cols_sample].mean(axis=1)
            join_samples.append(mean_col)

        df_mean_final = pd.concat(join_samples, axis=1)
        df_mean_final.columns = self.samples
        df_mean_final.index = self.data.index

        df_mean_final.to_csv(self.data_file)

        self.data = df_mean_final

        return self.data

    def join_berry_tissues(self) -> pd.DataFrame:
        """
        Joins the expressions for all berry tissues (flesh, skin, pericarp) into one sample by calculating the mean.
        Returns
        -------
        self.data: pd.Dataframe
           mean values of expression data for berry samples
        """
        other_tissues = [x for x in self.samples if 'berry' not in x]

        df_new = self.data.loc[:, other_tissues]

        berry_samples = ['_PFS', '_R']

        for sample in berry_samples:
            cols_sample = [x for x in self.data.columns if x.startswith('berry') and x.endswith(sample)]
            mean_col = self.data[cols_sample].mean(axis=1)
            df_new['berry' + sample] = mean_col

        df_new.to_csv(self.data_file)

        self.data = df_new
        self.samples = list(self.data.columns)
        return self.data

    def calculate_log2(self):
        """
        Converts the expression values to log2 values
        Returns
        -------
        self.data: pd.Dataframe
            log2 of expression data for all samples
        """
        log2 = np.log2(self.data)

        log2.to_csv(self.data_file)
        self.data = log2

        return self.data

    def select_samples(self):
        all_data = pd.read_csv(os.path.join(OMICS_DATADIR, self.dataset_id, self.dataset_id + '_vespucci.csv'),
                               index_col=0)
        samples = pd.read_excel(os.path.join(OMICS_DATADIR, self.dataset_id, self.dataset_id + '_metadata.xlsx'))
        sample_names = [x + '.ch1' for x in samples['sample'].to_list()]

        df_filtered = all_data[sample_names]
        df_filtered.columns = [c.replace('.ch1', '') for c in df_filtered.columns]
        df_filtered.to_csv(self.data_file)


def join_feature_counts():
    datafolder = os.path.join(OMICS_DATADIR, 'RNASeq', 'UPDATE', 'RAW_counts')
    project1 = os.path.join(datafolder, 'FeatureCounts_PRJNA386889')
    project2 = os.path.join(datafolder, 'FeatureCounts_PRJNA383160')

    metafile = os.path.join(OMICS_DATADIR, 'RNASeq', 'UPDATE', 'RNASeq_metadata.xlsx')
    metadata = pd.read_excel(metafile, index_col=0)
    allsamples = metadata.index

    tpm_file = pd.read_csv(os.path.join(OMICS_DATADIR, 'RNASeq', 'UPDATE', 'GREAT_TPM_MODEL_GENES.csv'),
                           index_col=0)
    genes = tpm_file.index

    files1 = os.listdir(project1)
    files2 = os.listdir(project2)

    all_counts = pd.DataFrame()

    for f in files1:
        sample = f[15:25]
        if sample in allsamples:
            counts = pd.read_csv(os.path.join(project1, f), sep="\t", skiprows=1, index_col=0)
            all_counts[sample] = counts.iloc[:, 5]

    for f in files2:
        sample = f[15:25]
        if sample in allsamples:
            counts = pd.read_csv(os.path.join(project2, f), sep="\t", skiprows=1, index_col=0)
            all_counts[sample] = counts.iloc[:, 5]

    all_counts.index = counts.index

    to_keep = set(all_counts.index).intersection(set(genes))

    all_counts = all_counts.loc[to_keep, :]
    all_counts.to_csv(os.path.join(datafolder, 'raw_counts_model_genes.csv'))


def replicate_same_ids() -> pd.DataFrame:
    """
    Get the v4 annotation for all the proteins for each gene

    Returns
    -------
    self.data: pd.Dataframe
        dataset with v4 ids as index
    """

    model = read_sbml_model('C:/Users/BiSBII/Documents/plantdb/reconstruction_results/vvinif2023/'
                            'vvinif2023_FINAL.xml')

    data = pd.read_csv('C:/Users/BiSBII/Documents/MM_ML/data/GREAT_LOG_TPM_GSE98923_MODEL_GENES.csv',
                       index_col=0)

    dic_genes = {}

    for gene in model.genes:
        new_gene_id, prot_id = gene.id.split('_')
        if new_gene_id not in dic_genes:
            dic_genes[new_gene_id] = [prot_id]
        else:
            dic_genes[new_gene_id].append(prot_id)

    new_rows = []
    new_index = []

    for ind, row in data.iterrows():
        replicates = dic_genes[ind]

        for rep in replicates:
            new_ind = str(ind) + '_' + rep
            new_index.append(new_ind)
            new_rows.append(row)

    new_df = pd.DataFrame(new_rows, columns=data.columns)
    new_df.index = new_index

    new_df.to_csv('C:/Users/BiSBII/Documents/MM_ML/data/GREAT_LOG_TPM_GSE98923_MODEL_PROTEINS.csv')

    data = new_df

    return data


def remove_replicates() -> pd.DataFrame:
    """
    Removes the replicates for each sample by calculating the mean.
    Returns
    -------
    self.data: pd.Dataframe
       mean values of expression data for each sample without replicates
    """
    data = pd.read_csv('C:/Users/BiSBII/Documents/MM_ML/data/GREAT_LOG_TPM_GSE98923_MODEL_GENES.csv',
                       index_col=0)

    metadata = pd.read_excel('C:/Users/BiSBII/Documents/MM_ML/data/GSE98923_metadata.xlsx')

    mean_col = data[['SRR5560551', 'SRR5560552']].mean(axis=1)
    data['SRR5560551'] = mean_col

    metadata = metadata[metadata['ENA'] != 'SRR5560552']

    tissues = set(metadata['tissue'])
    samples = []
    col_names = []
    for tissue in tissues:
        sample_ids = metadata['ENA'][metadata['tissue'] == tissue].to_list()
        descp = metadata['description'][metadata['tissue'] == tissue].to_list()
        cultivar = metadata['cultivar_abv'][metadata['tissue'] == tissue].to_list()
        year = metadata['year'][metadata['tissue'] == tissue].to_list()
        i = 0

        while i < len(sample_ids):
            chunk = sample_ids[i:i + 3]
            samples.append(chunk)
            new_name = cultivar[i] + '_' + descp[i][:-5].replace(' ', '') + '_' + str(year[i])
            col_names.append(new_name)
            i += 3

    mean_rep_cols = []
    for cols in samples:
        for col in cols:
            if col not in data.columns:
                cols.remove(col)

        mean_col = data[cols].mean(axis=1)
        mean_rep_cols.append(mean_col)
    print(samples)
    df_mean_noreps = pd.concat(mean_rep_cols, axis=1)
    df_mean_noreps.columns = col_names
    df_mean_noreps.index = data.index

    df_mean_noreps.to_csv('C:/Users/BiSBII/Documents/MM_ML/data/GREAT_LOG_TPM_GSE98923_MODEL_PROTEINS_NOREPS2.csv')
    return df_mean_noreps
