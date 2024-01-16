import sys
import os
import warnings
import traceback

import pandas

# Warnings and prints are blocked to avoid massive outputs that appear after running COBAMP.
warnings.simplefilter("ignore")


# Disable
def block_print():
    sys.stdout = open(os.devnull, 'w')


# Restore
def enable_print():
    sys.stdout = sys.__stdout__


block_print()

import cobra
import re

from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ReconstructionWrapper
from cobamp.utilities.parallel import batch_run

enable_print()

NOMENCLATURE = 'ensemble_gene_id'
OMICS_TYPE = 'transcriptomics'
AND_OR_FUNC = (min, max)


def reconstruction_function(omics_container, parameters: dict):
    """
    This function is used to run the reconstruction algorithm.

    Parameters
    ----------
    omics_container : pandas.DataFrame
        The omics_integration data set.
    parameters : dict
        The parameters to be used for the reconstruction algorithm.

    Returns
    ----------
    rec_wrapper : Reconstruction Wrapper object with model and omics_integration data.
    """
    block_print()

    threshold, rec_wrapper, method, protected = [parameters[parameter] for parameter in
                                                 ['threshold', 'reconstruction_wrapper', 'algorithm',
                                                  'protected_reactions']]

    def integration_fx(data_map):
        return [[k for k, v in data_map.get_scores().items() if (v is not None and v > threshold) or k in protected]]

    def score_apply(data_map):
        dm = {k: 0 if v is None else (min(v, 10) - threshold) if k not in protected else 20
              for k, v in data_map.items()}
        return dm

    # noinspection PyBroadException
    try:
        if method == 'fastcore':
            return rec_wrapper.run_from_omics(omics_data=omics_container, algorithm=method, and_or_funcs=AND_OR_FUNC,
                                              integration_strategy=('custom', [integration_fx]), solver='CPLEX')

        elif method == 'tinit':
            return rec_wrapper.run_from_omics(omics_data=omics_container, algorithm=method, and_or_funcs=AND_OR_FUNC,
                                              integration_strategy=('continuous', score_apply), solver='CPLEX')

    except:
        traceback.print_exc()

        return {r: False for r in rec_wrapper.model_reader.r_ids}


def troppo_omics_integration(model: cobra.Model, algorithm: str, threshold: float, thread_number: int,
                             omics_dataset: pandas.DataFrame, dataset: str,
                             protected_reacs: list = None):
    """
    This function is used to run the Troppo's integration algorithms.

    Parameters
    ----------
    omics_dataset: pandas.DataFrame
        A dataframe containing the omics_integration dataset to integrate
    model: cobra.Model
        The COBRA model.
    algorithm: str
        The algorithm to be used.
    threshold: float
        The threshold to be used.
    thread_number: int
        The number of threads to be used.
    dataset: str
        name of the dataset to integrate
    protected_reacs: list, optional
            list of reactions to keep in the context-specific models

    Returns
    -------
    integration_results: dict
        Dataframe containing the results of the omics_integration integration.
        Each sample as a dictionary containing a boolean value for each reaction.

    """
    block_print()

    patt = re.compile('__COBAMPGPRDOT__[0-9]{1}')
    replace_alt_transcripts = lambda x: patt.sub('', x)

    details = [dataset, algorithm, threshold]

    template = model.copy()

    omics_data = TabularReader(path_or_df=omics_dataset, nomenclature=NOMENCLATURE,
                               omics_type=OMICS_TYPE).to_containers()

    print(omics_data)

    enable_print()
    print('Tabular Reader Finished.')
    block_print()

    reconstruction_wrapper = ReconstructionWrapper(model=template, ttg_ratio=9999,
                                                   gpr_gene_parse_function=replace_alt_transcripts)

    enable_print()
    print('Reconstruction Wrapper Finished.')
    block_print()

    parameters = {'threshold': threshold, 'reconstruction_wrapper': reconstruction_wrapper, 'algorithm': algorithm}

    if protected_reacs:
        parameters['protected_reactions'] = protected_reacs
    else:
        parameters['protected_reactions'] = []

    batch_fastcore_res = batch_run(reconstruction_function, omics_data, parameters, threads=thread_number)

    result_dict = dict(zip([sample.condition for sample in omics_data], batch_fastcore_res))

    enable_print()
    print(f'Omics Integration with {details[1]} (Threshold = {details[2]}) Finished.')

    return result_dict
