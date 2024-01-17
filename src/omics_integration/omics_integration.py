import os.path
from typing import List, Optional, Literal
import pandas as pd
import cobra.io

from utils.config import PROJECT_PATH
from omics_dataset import OmicsDataset
from omics_processing import thresholding_filter
from troppo_integration import troppo_omics_integration


class OmicsIntegration:
    """
    Class to perform omics integration in metabolic models

    This class was changed to reconstruct one model at a time, as each model had a different biomass equation and
    medium conditions.

    All tissue samples are included in the omics dataset to calculate local filters, but then each sample is used
    to create each model individually, defining the objective and medium of each tissue

    If all tissues have the same objective and medium, reconstructing all models at once is way faster.
    """

    def __init__(self, dataset_id: str, model_id: str, model_file: str, obj_reaction: dict, medium: dict,
                 integration_thresholds: List[int], algorithm: Literal['fastcore', 'tinit'] = 'fastcore',
                 integration_strategy: Literal['default', 'Global', 'Local1', 'Local2'] = 'default',
                 protected_reacs: dict = None, global_threshold_lower: int = 0, samples: list = None,
                 global_threshold_upper: int = 3, local_threshold: int = 4, n_threads=4, macros: dict = None):
        """
        Initializes an Omics integration instance
        Parameters
        ----------
        dataset_id: str
            id of the dataset to use

        model_id: str
            id of the model to use

        obj_reaction: dict
            name of the reaction to use as objective for each sample

        medium: dict
            dictionary defining the bounds for the medium reactions (drains)

        integration_thresholds: list
            thresholds to be used by the integration algorithms

        algorithm: str (default: fastcore)
            algorithm for omics_integration integration. There are two possible algorithms: fastcore and tinit.

        integration_strategy: str (default: default)
            Thresholding strategy for omics_integration integration.
             - default: Does not use any local thresholding strategy to filter the omics_integration data;
             - Global: Uses only a global threshold;
             - Local1: Uses a global and a local threshold;
             - Local2: Uses two global thresholds, one lower and one upper, and a local threshold.

        protected_reacs: dict, optional
            list of reactions to keep in the model for each sample

        local_threshold: int (default = 4)
            index for the local threshold used by local strategies.
            The local threshold will correspond to a percentil value of the dataset.
            Currently, the index can get an int value between 0 and 4, corresponding to
            the percentil options [0.1, 0.25, 0.5, 0.75, 0.9].

        global_threshold_lower: int (default = 0)
            index for the global threshold used by global and local strategies.
            Same percentil options as for local_threshold

        global_threshold_upper: int (default = 3)
            index for the global threshold used by the local2 strategy.
            Same percentil options as for local_threshold

        n_threads: int (default = 4)
            number of threads to run the algorithm

        samples: list, optional
            samples to build a model for. If None, all samples in the dataset are considered

        macros: list, optional
            macro reactions that are different for each sample
        """

        self.dataset_id = dataset_id
        self.model_id = model_id
        self.algorithm = algorithm
        self.objective = obj_reaction
        self.medium = medium
        self.int_thresholds = integration_thresholds
        self.int_strategy = integration_strategy
        self.protected_reacs = protected_reacs
        self.global_threshold_lower = global_threshold_lower
        self.global_threshold_upper = global_threshold_upper
        self.local_threshold = local_threshold
        self.n_threads = n_threads
        self.macros = macros

        dataset_object = OmicsDataset(dataset_id=dataset_id)

        if samples:
            self.samples = samples
        else:
            self.samples = dataset_object.samples

        self.omics_data = dataset_object.data.transpose()

        self.model_file = os.path.join(PROJECT_PATH, 'reconstruction_results', self.model_id, model_file + '.xml')

        self.results_path = os.path.join(PROJECT_PATH, 'reconstruction_results', self.model_id, 'results_troppo',
                                         self.dataset_id)

        if not os.path.exists(self.results_path):
            os.mkdir(self.results_path)

    def load_model(self, obj_function, model_medium, sample) -> Optional[cobra.Model]:
        """
        This method is used to load the model and update the medium.
        Parameters
        ---------
        obj_function: str
            objective function to the model
        model_medium: dict
            model_medium exchange reactions for this sample model
        sample: str
            name of the sample
        Returns
        -------
        model : cobra.Model
            The loaded model.
        """
        if os.path.isfile(self.model_file):
            cobra_model = cobra.io.read_sbml_model(self.model_file)
        else:
            print('The model file does not exists')
            return None

        sample_model = cobra_model.copy()

        sample_model.objective = obj_function

        for reaction_id, bounds in model_medium.items():
            try:
                r_obj = sample_model.reactions.get_by_id(reaction_id)
                r_obj.bounds = bounds
            except KeyError:
                pass

        to_remove = [v for k, v in self.objective.items() if k != sample]

        if sample not in self.macros:
            for reac_list in self.macros.values():
                to_remove.extend(reac_list)
        else:
            keep = self.macros[sample]
            for reac_list in self.macros.values():
                for reac in reac_list:
                    if reac not in keep:
                        to_remove.append(reac)

            other = [r.replace(sample, '') for r in keep]
            to_remove.extend(other)

        sample_model.remove_reactions(to_remove)

        return sample_model

    def local_filtering(self) -> pd.DataFrame:
        """
        This method performs local threshold filtering of the omics_integration dataset
        Returns
        -------
        self.omics_data: pd.Dataframe
            omics_integration data filtered
        """
        self.omics_data = thresholding_filter(omics_dataframe=self.omics_data,
                                              thresholding_strategy=self.int_strategy,
                                              global_threshold_lower=self.global_threshold_lower,
                                              global_threshold_upper=self.global_threshold_upper,
                                              local_threshold=self.local_threshold)
        return self.omics_data

    def model_reconstruction(self, model_temp: cobra.Model, sample: str, integration_result_dict: dict):
        """
        This function is used to reconstruct the model based on the troppo integration results.

        Parameters
        ----------
        model_temp: cobra.Model
            The COBRA model template.
        sample: str
            The sample name.
        integration_result_dict: dict
            The integration results from troppo.
        """
        temp_model = model_temp.copy()

        reactions_to_deactivate = [reaction for reaction, value in
                                   integration_result_dict[sample].items() if value is False]

        for reaction in reactions_to_deactivate:
            temp_model.remove_reactions([reaction])

        model_name = (f'{sample}_{self.int_strategy}_{self.algorithm}_{self.global_threshold_lower}_'
                      f'{self.global_threshold_upper}_{self.local_threshold}.xml')

        new_folder = os.path.join(self.results_path, 'reconstructed_models')

        if not os.path.exists(new_folder):
            os.mkdir(new_folder)

        output_model = os.path.join(self.results_path, 'reconstructed_models', model_name)
        cobra.io.write_sbml_model(temp_model, output_model)

        print(f'Model reconstruction for {sample} finished.')
        print_model_details(temp_model)

    def reconstruction_pipeline(self):
        """
        This function is used to run the reconstruction pipeline.

        In the end, this function generates a SBML file of the tissue-specific models that resulted from the
        omics_integration data integration with Troppo.
        """

        print('-------------------------------------------------------------------------------------------------------')
        print('-------------------------------------- Processing Omics Dataset. --------------------------------------')
        print('-------------------------------------------------------------------------------------------------------')

        if self.int_strategy != 'default':
            self.omics_data = self.local_filtering()

            print(f'{self.int_strategy} threshold filter applied.')

        print('-------------------------------------------------------------------------------------------------------')
        print('----------------------- Starting Omics Integration with Troppo for each sample. -----------------------')
        print('-------------------------------------------------------------------------------------------------------')

        troppo_result_dict = {}

        sample_names = []
        for sample in self.samples:
            objective_function = self.objective[sample]
            sample_medium = self.medium[sample]
            template_model = self.load_model(obj_function=objective_function, model_medium=sample_medium, sample=sample)

            drains = [d.id for d in template_model.reactions if d.id.startswith('EX_')]

            transporters = [t.id for t in template_model.reactions if t.id.startswith('T_')]

            sample_protected = self.protected_reacs[sample] + drains + transporters

            sample_omics = self.omics_data.loc[sample, :]
            sample_omics_df = pd.DataFrame(sample_omics, columns=[sample])
            sample_omics_df = sample_omics_df.transpose()

            if not template_model:
                return None

            print(f'Omics integration with {self.algorithm} started for sample {sample}.')
            print('---------------------------------------------------------------------------------------------------')

            for threshold in self.int_thresholds:
                print(f'Threshold = {threshold}')
                print('-----------------------------------------------------------------------------------------------')
                troppo_result = troppo_omics_integration(model=template_model, algorithm=self.algorithm,
                                                         threshold=threshold, thread_number=self.n_threads,
                                                         omics_dataset=sample_omics_df, dataset=self.dataset_id,
                                                         protected_reacs=sample_protected)

                th = str(round(threshold, 2)).replace('.', '_')
                troppo_result_dict[f'{sample}_t{th}'] = troppo_result[sample]
                if sample not in sample_names:
                    sample_names.append(sample)

            print('---------------------------------------------------------------------------------------------------')

            for key in troppo_result_dict:
                if sample in key:

                    print(f'Context-specific model reconstruction for {key} started.')
                    print('-------------------------------------------------------------------------------------------')

                    self.model_reconstruction(model_temp=template_model, sample=key,
                                              integration_result_dict=troppo_result_dict)

        sample_names = '_'.join(sample_names)

        if self.int_strategy == 'default':
            result_file = os.path.join(self.results_path,
                                       f'{sample_names}_{self.algorithm}_{self.int_strategy}.csv')
        else:
            file_path = f'{sample_names}_{self.algorithm}_{self.int_strategy}_' \
                        f'{self.global_threshold_lower}_{self.global_threshold_upper}_{self.local_threshold}.csv'
            result_file = os.path.join(self.results_path, file_path)

        troppo_result_dataframe = pd.DataFrame.from_dict(troppo_result_dict, orient='index')
        troppo_result_dataframe.to_csv(result_file)

        print('------------------------------------------ Pipeline Finished ------------------------------------------')
        print('-------------------------------------------------------------------------------------------------------')


def print_model_details(cobra_model: cobra.Model):
    """
    Function to print the details of the currently loaded COBRA model.

    Parameters
    ----------
    cobra_model : cobra.Model

    """
    transporters = []

    for reac in cobra_model.reactions:
        if len(reac.compartments) == 2:
            transporters.append(reac.id)

    print('Total Reactions:', len(cobra_model.reactions))
    print('Reactions:', (len(cobra_model.reactions)) - len(transporters) - len(cobra_model.exchanges))
    print('Transporters:', len(transporters))
    print('Exchanges:', len(cobra_model.exchanges))


if __name__ == '__main__':
    dataset = 'RNAseq'
    modelid = 'vvinif2023'
    modelfile = 'vvinif2023_FINAL'

    obj = {'leaf': 'e-Biomass_vvinif2023_leaf__cyto', 'stem': 'e-Biomass_vvinif2023_stem__cyto',
           'berry_green': 'e-Biomass_vvinif2023__cyto', 'berry_mature': 'e-Biomass_vvinif2023_berry_mature__cyto'}

    protec = {'leaf': ['e-Biomass_vvinif2023_leaf__cyto'], 'stem': ['e-Biomass_vvinif2023_stem__cyto'],
              'berry_green': ['e-Biomass_vvinif2023__cyto'],
              'berry_mature': ['e-Biomass_vvinif2023_berry_mature__cyto']}

    macro = {'leaf': ['e-Cofactor_vvinif2023_leaf__cyto'],
             'berry_mature': ['e-Carbohydrates_vvinif2023_berry_mature__cyto']}

    all_samples = ['berry_green', 'stem', 'berry_mature', 'leaf']

    med_photo = {'EX_Light_drain': (-100, 10000), 'EX_CARBON-DIOXIDE_drain': (-10000, 10000),
                 'EX_SUCROSE_drain': (0, 10000), 'EX_FE+2_drain': (-10000, 10000),
                 'EX_MG+2_drain': (-10000, 10000), 'EX_WATER_drain': (-10000, 10000),
                 'EX_OXYGEN-MOLECULE_drain': (-10000, 10000), 'EX_PROTON_drain': (-10000, 10000),
                 'EX_Pi_drain': (-10000, 10000), 'EX_NITRATE_drain': (-10000, 10000),
                 'EX_SULFATE_drain': (-10000, 10000), 'EX_ARG_drain': (0, 10000),
                 'EX_L-ASPARTATE_drain': (0, 10000), 'EX_TYR_drain': (0, 10000),
                 'EX_PHE_drain': (0, 10000), 'EX_TRP_drain': (0, 10000), 'EX_ASN_drain': (0, 10000),
                 'EX_SER_drain': (0, 10000), 'EX_GLN_drain': (0, 10000), 'EX_GLY_drain': (0, 10000),
                 'EX_GLT_drain': (0, 10000), 'EX_THR_drain': (0, 10000), 'EX_LEU_drain': (0, 10000),
                 'EX_VAL_drain': (0, 10000), 'EX_ILE_drain': (0, 10000), 'EX_HIS_drain': (0, 10000),
                 'EX_MET_drain': (0, 10000), 'EX_LYS_drain': (0, 10000),
                 'EX_L-ALPHA-ALANINE_drain': (0, 10000)}

    med_respiration = {'EX_Light_drain': (0, 10000), 'EX_CARBON-DIOXIDE_drain': (0, 10000),
                       'EX_SUCROSE_drain': (-1, 10000), 'EX_FE+2_drain': (-10000, 10000),
                       'EX_MG+2_drain': (-10000, 10000), 'EX_WATER_drain': (-10000, 10000),
                       'EX_OXYGEN-MOLECULE_drain': (-10000, 10000), 'EX_PROTON_drain': (-10000, 10000),
                       'EX_Pi_drain': (-10000, 10000), 'EX_NITRATE_drain': (-10000, 10000),
                       'EX_SULFATE_drain': (-10000, 10000), 'EX_ARG_drain': (0, 10000),
                       'EX_L-ASPARTATE_drain': (0, 10000), 'EX_TYR_drain': (0, 10000),
                       'EX_PHE_drain': (0, 10000), 'EX_TRP_drain': (0, 10000), 'EX_ASN_drain': (0, 10000),
                       'EX_SER_drain': (0, 10000), 'EX_GLN_drain': (0, 10000), 'EX_GLY_drain': (0, 10000),
                       'EX_GLT_drain': (0, 10000), 'EX_THR_drain': (0, 10000), 'EX_LEU_drain': (0, 10000),
                       'EX_VAL_drain': (0, 10000), 'EX_ILE_drain': (0, 10000), 'EX_HIS_drain': (0, 10000),
                       'EX_MET_drain': (0, 10000), 'EX_LYS_drain': (0, 10000),
                       'EX_L-ALPHA-ALANINE_drain': (0, 10000)}

    mediums = {'leaf': med_photo, 'stem': med_respiration, 'berry_mature': med_respiration,
               'berry_green': med_respiration}

    # ################################################# Default ##################################################

    default = OmicsIntegration(dataset_id=dataset, model_id=modelid, model_file=modelfile, obj_reaction=obj,
                               protected_reacs=protec, integration_thresholds=[6, 8, 10], medium=mediums,
                               integration_strategy='default', samples=all_samples, macros=macro)

    default.reconstruction_pipeline()

    # ################################################# Global ##################################################

    glb_50 = OmicsIntegration(dataset_id=dataset, model_id=modelid, model_file=modelfile, obj_reaction=obj,
                              protected_reacs=protec, integration_thresholds=[0], medium=mediums,
                              integration_strategy='Global', global_threshold_upper=2, samples=all_samples,
                              macros=macro)

    glb_50.reconstruction_pipeline()

    glb_75 = OmicsIntegration(dataset_id=dataset, model_id=modelid, model_file=modelfile, obj_reaction=obj,
                              protected_reacs=protec, integration_thresholds=[0], medium=mediums,
                              integration_strategy='Global', global_threshold_upper=3, samples=all_samples,
                              macros=macro)

    glb_75.reconstruction_pipeline()

    # # ################################################## Local1 ##################################################

    thresholds = [[2, 2], [2, 3], [3, 2], [3, 3]]

    for ths in thresholds:

        local1 = OmicsIntegration(dataset_id=dataset, model_id=modelid, model_file=modelfile,
                                  obj_reaction=obj, protected_reacs=protec, integration_thresholds=[0], medium=mediums,
                                  integration_strategy='Local1', global_threshold_upper=ths[0], local_threshold=ths[1],
                                  samples=all_samples, macros=macro)

        local1.reconstruction_pipeline()

    # # ################################################## Local2 ##################################################

    thresholds = [[1, 3, 2], [1, 3, 3]]

    for ths in thresholds:

        local2 = OmicsIntegration(dataset_id=dataset, model_id=modelid, obj_reaction=obj, model_file=modelfile,
                                  protected_reacs=protec, integration_thresholds=[0], medium=mediums,
                                  integration_strategy='Local2', samples=all_samples, global_threshold_lower=ths[0],
                                  global_threshold_upper=ths[1], local_threshold=ths[2], macros=macro)

        local2.reconstruction_pipeline()
