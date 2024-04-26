import os
from configparser import RawConfigParser

configs = RawConfigParser()
configs.read('/iplants_app/conf/iplantsapp.conf')

PROJECT_PATH = '/iplants_app'
OMICS_DATADIR = os.path.join(PROJECT_PATH, 'omics_data')


class TransytAPI:
    run_link = str(configs.get('iplantsapp-configurations', 'transyt_run_link'))
    status_link = str(configs.get('iplantsapp-configurations', 'transyt_status_link'))
    download_link = str(configs.get('iplantsapp-configurations', 'transyt_download_link'))


class DiamondAPI:
    diamondp_run_link = str(configs.get('iplantsapp-configurations', 'diamondp_run_link'))
    diamondx_run_link = str(configs.get('iplantsapp-configurations', 'diamondx_run_link'))
    status_link = str(configs.get('iplantsapp-configurations', 'diamond_status_link'))
    download_link = str(configs.get('iplantsapp-configurations', 'diamond_download_link'))


class API:
    """
    Change API base_url when updated
    """
    base_url = str(configs.get('iplantsapp-configurations', 'iplants_api'))
    create_node_model = str(configs.get('iplantsapp-configurations', 'create_node_model'))
    create_doc_model = str(configs.get('iplantsapp-configurations', 'create_doc_model'))
    reactions_from_enzyme_neo = str(configs.get('iplantsapp-configurations', 'reactions_from_enzyme_neo'))
    metabolites_from_reaction = str(configs.get('iplantsapp-configurations', 'metabolites_from_reaction'))
    pathways_from_reaction = str(configs.get('iplantsapp-configurations', 'pathways_from_reaction'))
    enzymes_from_reaction_mongo = str(configs.get('iplantsapp-configurations',
                                                  'enzymes_from_reaction_mongo'))
    components_from_enzyme = str(configs.get('iplantsapp-configurations', 'components_from_enzyme'))
    reactions_from_model = str(configs.get('iplantsapp-configurations', 'reactions_from_model'))
    model_detail = str(configs.get('iplantsapp-configurations', 'model_detail'))
    reaction_detail = str(configs.get('iplantsapp-configurations', 'reaction_detail'))
    pathway_detail = str(configs.get('iplantsapp-configurations', 'pathway_detail'))
    metabolite_detail = str(configs.get('iplantsapp-configurations', 'metabolite_detail'))
    organism_detail = str(configs.get('iplantsapp-configurations', 'organism_detail'))
    add_enzymes_model_neo = str(configs.get('iplantsapp-configurations', 'add_enzymes_model_neo'))
    add_reactions_model_neo = str(configs.get('iplantsapp-configurations', 'add_reactions_model_neo'))
    add_metabolites_model_neo = str(configs.get('iplantsapp-configurations', 'add_metabolites_model_neo'))
    add_enzymes_model_mongo = str(configs.get('iplantsapp-configurations', 'add_enzymes_model_mongo'))
    add_reactions_model_mongo = str(configs.get('iplantsapp-configurations', 'add_reactions_model_mongo'))
    add_metabolites_model_mongo = str(configs.get('iplantsapp-configurations',
                                                  'add_metabolites_model_mongo'))
    add_pathways_model_mongo = str(configs.get('iplantsapp-configurations', 'add_pathways_model_mongo'))
    add_genes_model_mongo = str(configs.get('iplantsapp-configurations', 'add_genes_model_mongo'))
    add_annotation_model = str(configs.get('iplantsapp-configurations', 'add_annotation_model'))
    add_gprs_model = str(configs.get('iplantsapp-configurations', 'add_gprs_model'))
    add_comparts_model = str(configs.get('iplantsapp-configurations', 'add_comparts_model'))
    add_new_reactions_mongo = str(configs.get('iplantsapp-configurations', 'add_new_reactions_mongo'))
    add_new_reactions_neo = str(configs.get('iplantsapp-configurations', 'add_new_reactions_neo'))
