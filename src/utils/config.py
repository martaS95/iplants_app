import os

PROJECT_PATH = 'C:/Users/BiSBII/Documents/plantdb'
OMICS_DATADIR = os.path.join(PROJECT_PATH, 'omics_data')


class API:
    base_url = "http://127.0.0.1:8000"
    create_node_model = 'modelneo'
    create_doc_model = 'modeldoc'
    reactions_from_enzyme_neo = '/reactions/enzyme/'
    metabolites_from_reaction = '/metabolites/reaction/'
    pathways_from_reaction = '/pathways/reaction/'
    enzymes_from_reaction_mongo = '/enzymes/reaction/'
    components_from_enzyme = '/components/enzyme/'
    add_enzymes_model_neo = '/modelenzymes_neo/'
    add_reactions_model_neo = '/modelreactions_neo/'
    add_metabolites_model_neo = '/modelmetabolites_neo/'
    add_enzymes_model_mongo = '/modelenzymes/'
    add_reactions_model_mongo = '/modelreactions/'
    add_metabolites_model_mongo = '/modelmetabolites/'
    add_pathways_model_mongo = '/modelpathways/'
    add_genes_model_mongo = '/modelgenes/'
    add_annotation_model = '/annotation/'
    add_gprs_model = '/gprs/'
    add_comparts_model = '/compartments/'
    add_new_reactions_mongo = '/newreactions/'
    add_new_reactions_neo = '/newreactions_neo/'
    reactions_from_model = '/reactions/model/'
    model_detail = '/metabolicmodel/'
    reaction_detail = '/reaction/'
    pathway_detail = '/pathway/'
    metabolite_detail = '/metabolite/'
    organism_detail = '/organism/'


class Mongo:
    host = 'palsson.di.uminho.pt'
    port = 1017
    database = 'plantcyc'


class Neo:
    host = 'palsson.di.uminho.pt'
    port = 1087
    username = 'neo4j'
    password = 'plant'
