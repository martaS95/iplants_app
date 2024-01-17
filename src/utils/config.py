import os

PROJECT_PATH = 'C:/Users/BiSBII/Documents/iplants_app'
OMICS_DATADIR = os.path.join(PROJECT_PATH, 'omics_data')


class API:
    """
    Change API base_url when updated
    """
    base_url = "http://127.0.0.1:8000"
    create_node_model = '/create_model/model_node/'
    create_doc_model = '/create_model/model_doc/'
    reactions_from_enzyme_neo = '/list/reactions/enzyme/'
    metabolites_from_reaction = '/list/metabolites/reaction/'
    pathways_from_reaction = '/list/pathways/reaction/'
    enzymes_from_reaction_mongo = '/list/enzymes/reaction/doc/'
    components_from_enzyme = '/list/components/enzyme/'
    add_enzymes_model_neo = '/create_model/add_rels_enzymes/'
    add_reactions_model_neo = '/create_model/add_rels_reactions/'
    add_metabolites_model_neo = '/create_model/add_rels_metabolites/'
    add_enzymes_model_mongo = '/create_model/add_enzymes/'
    add_reactions_model_mongo = '/create_model/add_reactions/'
    add_metabolites_model_mongo = '/create_model/add_metabolites//'
    add_pathways_model_mongo = '/create_model/add_pathways/'
    add_genes_model_mongo = '/create_model/add_genes/'
    add_annotation_model = '/create_model/add_annotation/'
    add_gprs_model = '/create_model/add_gprs/'
    add_comparts_model = '/create_model/add_compartments/'
    add_new_reactions_mongo = '/create_model/newreactions/'
    add_new_reactions_neo = '/create_model/newreactions_neo/'
    reactions_from_model = '/reactions/model/'
    model_detail = '/detail/metabolicmodel/'
    reaction_detail = '/detail/reaction/'
    pathway_detail = '/detail/pathway/'
    metabolite_detail = '/detail/metabolite/'
    organism_detail = '/detail/organism/'


class Mongo:
    host = 'palsson.di.uminho.pt'
    port = 1017
    database = 'plantcyc'


class Neo:
    host = 'palsson.di.uminho.pt'
    port = 1087
    username = 'neo4j'
    password = 'plant'
