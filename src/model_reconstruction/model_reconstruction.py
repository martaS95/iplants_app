import os
import logging
import re
import json
import datetime
from collections import defaultdict
from typing import List, Union, Dict, Any, Tuple, Optional
import requests
import pandas as pd
from Bio import SeqIO
from cobra import Model as CobraModel
from cobra import Reaction as CobraReaction
from cobra import Metabolite as CobraMetabolite
from cobra.core import Group as CobraGroup
from cobra.io import write_sbml_model
from ebiomass import BiomassFromGenome
from diamond import Diamond
from compartments import Compartments
from transyt import Transyt
from src.utils.config import PROJECT_PATH, API

logging.basicConfig(level=logging.DEBUG)


class ModelReconstruction:

    def __init__(self, model_id: str, organism: str, taxid: int, year: int, author: str, email: str,
                 genbank_file: str = None):
        """
        Class implementing the steps for model reconstruction
        Parameters
        ----------
        model_id: str
            identifier for the metabolic model under reconstruction
        organism: str
            name of the organism
        taxid: int
            taxonomy id for the organism
        year: int
            year of the reconstruction of the model
        author: str
            name of the author of the model reconstruction
        email: str
            email of the author to run diamond
        genbank_file: str
            name of the genbank file of the genome
        """
        self.model_id = model_id
        self.organism = organism
        self.taxid = taxid
        self.year = year
        self.author = author
        self.email = email
        self.genbank_file = genbank_file
        self.genes = []
        self.reactions = []
        self.enzymes = []
        self.metabolites = []
        self.pathways = []
        self.gprs = []
        self.genbank_annotation = {}
        self.diamond_annotation = {}
        self.biomass_reactions = {}
        self.drains = {}
        self.prot_compartments = {}
        self.output_diamond_file = None

        self.project_path = PROJECT_PATH
        self.genome_path = os.path.join(self.project_path, 'genome_files', self.model_id)

        self.api = API()

    def create_model_in_database(self):
        """
        Create the metabolic model instance node in the Neo4j database and the document in MongoDB
        """

        url_elements_mongo = [self.api.base_url, self.api.create_doc_model, self.model_id, self.organism,
                              self.taxid, self.year, self.author]

        url_elements_neo = [self.api.base_url, self.api.create_node_model, self.model_id, self.organism,
                            self.taxid, self.year, self.author]

        url_mongo = '/'.join(str(x) for x in url_elements_mongo)

        url_neo = '/'.join(str(x) for x in url_elements_neo)

        r_mongo = requests.get(url_mongo)

        if r_mongo.status_code == 200:
            logging.info('Metabolic model document for ' + str(self.model_id) + ' was created!')
        else:
            logging.info('An error has occurred when connecting to the mongodb API')

        r_neo = requests.get(url_neo)

        if r_neo.status_code == 200:
            logging.info('Metabolic model node for ' + str(self.model_id) + ' was created!')
        else:
            logging.info('An error has occurred when connecting to the neo4j API')

    def get_gene_locus(self) -> Optional[List[str]]:
        """
        Get the locus tag from the genbank annotation of the genome. It requires the genbank_file
        Returns
        -------
        list_of_genes: list
            locus_tag of the genes in the genome
        """
        if self.genbank_file:
            filepath = os.path.join(self.genome_path, self.genbank_file)
        else:
            logging.info('There is no genbank file specified. It is not possible to get the gene locus of the genome.')
            return None

        list_of_genes = []
        records = SeqIO.parse(filepath, format='gb')
        for record in records:
            for feat in record.features:
                if feat.type == 'CDS':
                    if 'gene' in feat.qualifiers:
                        gene = feat.qualifiers['gene'][0]
                        if 'gene' not in list_of_genes:
                            list_of_genes.append(gene)
        self.genes = list_of_genes
        return self.genes

    def get_gene_proteins(self) -> Optional[Dict[str, list]]:
        """
        Get the locus tag of genes and the associated product (protein_id) from the annotation of the genome
        (genbank file is needed)
        Returns
        -------
        gene_protein_dic: dic
            genes on the genome and associated products (protein identifiers)
        """
        if self.genbank_file:
            filepath = os.path.join(self.genome_path, self.genbank_file)
        else:
            logging.info('There is no genbank file specified. It is not possible to get the gene locus of the genome.')
            return None

        gene_protein_dic = {}
        records = SeqIO.parse(filepath, format='gb')
        for record in records:
            for feat in record.features:
                if feat.type == 'CDS':
                    if 'gene' in feat.qualifiers:
                        gene = feat.qualifiers['gene'][0]
                    else:
                        continue
                    if 'protein_id' in feat.qualifiers:
                        protein = feat.qualifiers['protein_id'][0]
                    else:
                        protein = None

                    if gene not in gene_protein_dic and protein:
                        gene_protein_dic[gene] = [protein]
                    elif protein and protein not in gene_protein_dic[gene]:
                        gene_protein_dic[gene].append(protein)

        self.genbank_annotation = gene_protein_dic
        return self.genbank_annotation

    def run_diamond(self, runtype: str = 'diamondp') -> Dict[str, str]:
        """
        Runs diamondP or DiamondX. Requires protein.faa or genomic.fna files and the SubmissionParameters.txt to run.
        It imports the class Diamond that runs the diamond in a container on the server and parses the result to a json
        Returns
        -------
        diamond_annotation: dict
            {genomic identifier: database identifier}
        """
        output_folder = os.path.join(self.project_path, 'reconstruction_results', self.model_id, 'results_diamond')

        if not os.path.isdir(os.path.join(output_folder)):
            os.makedirs(os.path.join(output_folder))

        self.output_diamond_file = os.path.join(output_folder, 'results_' +
                                                str(datetime.datetime.now().strftime("%Y%m%d_%I%M")))

        diamond_cls = Diamond(source_folder=self.genome_path, email=self.email, output_folder=self.output_diamond_file,
                              run_type=runtype)

        diamond_cls.run_diamond()

        self.diamond_annotation = diamond_cls.parse_results()

        return self.diamond_annotation

    def get_enzymes_in_model(self) -> Optional[List[str]]:
        """
        Get enzymes in model from diamond results
        Returns
        self.enzymes: list
            enzymes in model
        -------
        """
        if self.diamond_annotation:
            enzs = self.diamond_annotation.values()
            enzs = [e for e in enzs if e != 'no hits']
            self.enzymes = list(set(enzs))
            return self.enzymes

    def get_reacs_from_enzymes(self) -> Optional[List[str]]:
        """
        Get the list of reactions in model from the list of enzymes
        Returns
        -------
        self.reactions: list
            reactions in the model
        """

        url = self.api.base_url + self.api.reactions_from_enzyme_neo

        if self.enzymes:
            reactions_in_model = []

            for enz in self.enzymes:
                url_enz = url + enz

                r = requests.get(url_enz)
                if r.status_code == 200:
                    res = r.json()
                    if 'reactions' in res:
                        reactions = res['reactions']
                        reactions_in_model.extend(reactions)
                else:
                    logging.info('An error has occurred when connecting to neo4j API')

            self.reactions = list(set(reactions_in_model))

            return self.reactions

        else:
            logging.info('No enzymes in the model')

    def get_mets_from_reactions(self) -> Optional[List[str]]:
        """
        Get the list of metabolites in model from the list of reactions
        Returns
        -------
        self.metabolites: list
            metabolites in the model
        """
        if self.reactions:
            mets_in_model = []

            url = self.api.base_url + self.api.metabolites_from_reaction

            for reac in self.reactions:

                url_reac = url + reac

                r = requests.get(url_reac)
                if r.status_code == 200:
                    res = r.json()
                    if 'metabolites' in res:
                        metabolites = res['metabolites']
                        mets_in_model.extend(metabolites)
                else:
                    logging.info('An error has occurred when connecting to neo4j API')

            self.metabolites = list(set(mets_in_model))
            return self.metabolites

        else:
            logging.info('No reactions in the model')

    def get_paths_from_reactions(self) -> Optional[List[str]]:
        """
        Get the list of pathways in model from the list of reactions
        Returns
        -------
        self.pathways: list
            pathways in the model
        """

        if self.reactions:
            paths_in_model = []

            url = self.api.base_url + self.api.pathways_from_reaction

            for reac in self.reactions:

                url_reac = url + reac

                r = requests.get(url_reac)
                if r.status_code == 200:
                    res = r.json()
                    if 'pathways' in res:
                        pathways = res['pathways']
                        paths_in_model.extend(pathways)
                else:
                    logging.info('An error has occurred when connecting to neo4j API')

            self.pathways = list(set(paths_in_model))
            return self.pathways

        else:
            logging.info('No reactions in the model')

    def get_model_elems(self):
        """
        Get all information about the model, including enzymes, reactions, metabolites and pathways
        """
        self.get_enzymes_in_model()
        self.get_reacs_from_enzymes()
        self.get_mets_from_reactions()
        self.get_paths_from_reactions()

    def add_enzymes_to_model_db(self):
        """
        Create the relationships between the metabolic model node and the enzymes in Neo4j and update the
        attribute 'enzymes' in the metabolic model document in mongodb
        """
        if self.enzymes:
            enzymes = ','.join(self.enzymes)

            url_neo = self.api.base_url + self.api.add_enzymes_model_neo + self.model_id + '/' + enzymes
            r_neo = requests.get(url_neo)

            if r_neo.status_code == 200:
                logging.info('the relationships between the metabolic model and enzymes were created!')
            else:
                logging.info('an error has occurred connectiong to the neo4j API')

            url_mongo = self.api.base_url + self.api.add_enzymes_model_mongo + self.model_id + '/' + enzymes
            r_mongo = requests.get(url_mongo)

            if r_mongo.status_code == 200:
                logging.info('the enzyme list was added to the metabolic model document')
            else:
                logging.info('an error has occurred connectiong to the mongo API')

        else:
            logging.info('No enzymes in the model')

    def add_reacs_to_model_db(self):
        """
        Create the relationships between the metabolic model node and the reactions in Neo4j and update the
        attribute 'reactions' in the metabolic model document in mongodb
        """
        if self.reactions:
            reactions = ','.join(self.reactions)

            url_neo = self.api.base_url + self.api.add_reactions_model_neo + self.model_id + '/' + reactions
            r_neo = requests.get(url_neo)

            if r_neo.status_code == 200:
                logging.info('the relationships between the metabolic model and reactions were created!')
            else:
                logging.info('an error has occurred connectiong to the neo4j API')

            url_mongo = self.api.base_url + self.api.add_reactions_model_mongo + self.model_id + '/' + reactions
            r_mongo = requests.get(url_mongo)

            if r_mongo.status_code == 200:
                logging.info('the reaction list was added to the metabolic model document')
            else:
                logging.info('an error has occurred connectiong to the mongo API')

        else:
            logging.info('No reactions in the model')

    def add_mets_to_model_db(self):
        """
        Create the relationships between the metabolic model node and the metabolites in Neo4j and update the
        attribute 'metabolites' in the metabolic model document in mongodb
        """
        if self.metabolites:
            metabolites = ','.join(self.metabolites)

            url_neo = self.api.base_url + self.api.add_metabolites_model_neo + self.model_id + '/' + metabolites
            r_neo = requests.get(url_neo)

            if r_neo.status_code == 200:
                logging.info('the relationships between the metabolic model and metabolites were created!')
            else:
                logging.info('an error has occurred connectiong to the neo4j API')

            url_mongo = self.api.base_url + self.api.add_metabolites_model_mongo + self.model_id + '/' + metabolites
            r_mongo = requests.get(url_mongo)

            if r_mongo.status_code == 200:
                logging.info('the metabolite list was added to the metabolic model document')
            else:
                logging.info('an error has occurred connectiong to the mongo API')

        else:
            logging.info('No metabolites in the model')

    def add_paths_to_model_db(self):
        """
        Add list of pathways to the metabolic model object in monogdb
        """
        if self.pathways:
            pathways = ','.join(self.pathways)
            url_mongo = self.api.base_url + self.api.add_pathways_model_mongo + self.model_id + '/' + pathways
            r_mongo = requests.get(url_mongo)

            if r_mongo.status_code == 200:
                logging.info('the pathway list was added to the metabolic model document')
            else:
                logging.info('an error has occurred connectiong to the mongo API')

        else:
            logging.info('No pathways in the model')

    def add_genes_to_model_db(self):
        """
        Add list of genes to the metabolic model object in monogdb
        """
        if self.genes:
            genes = ','.join(self.genes)
            url_mongo = self.api.base_url + self.api.add_genes_model_mongo + self.model_id + '/' + genes
            r_mongo = requests.get(url_mongo)

            if r_mongo.status_code == 200:
                logging.info('the gene list was added to the metabolic model document')
            else:
                logging.info('an error has occurred connectiong to the mongo API')

        else:
            logging.info('No genes in the model')

    def add_annotation_to_model_db(self):
        """
        Add results from diamond annotation to the metabolic model object in monogdb
        """
        if self.output_diamond_file:
            try:
                annotation_file = open(os.path.join(self.output_diamond_file, 'DiamondResults_parsed.json'), 'r')
            except FileNotFoundError:
                logging.info('No annotation file was found')
                return None

            url = self.api.base_url + self.api.add_annotation_model + self.model_id

            r = requests.post(url, files={'file': annotation_file})

            if r.status_code == 200:
                logging.info('annotation data was added to the model document')
            else:
                logging.info('An error has occurred when connecting to mongo API')
        else:
            logging.info('No annotation file is defined')

    def add_info_to_model_db(self):
        """
        Create all the model relationships in the neo4j database, including to enzymes, reactions and metabolites
        """
        self.add_enzymes_to_model_db()
        self.add_reacs_to_model_db()
        self.add_mets_to_model_db()
        self.add_paths_to_model_db()
        self.add_genes_to_model_db()
        self.add_annotation_to_model_db()

    def get_unbalanced_reactions(self) -> List[List[Union[str, dict]]]:
        """
        Get the unbalanced reactions in the model
        Returns
        -------
        list_of_reacs: list
            list of lists with [reaction_id, {atom: amount}] for each unbalaced reaction
        """
        if self.reactions:
            list_of_reacs = []
            for reac in self.reactions:
                balance_dic = self.check_mass_balance(reacid=reac)
                if balance_dic:
                    list_of_reacs.append([reac, balance_dic])
            return list_of_reacs
        else:
            logging.info('No reactions in the model')

    def create_biomass_reaction(self, atp_coef: float = 53.26) -> Optional[Dict[str, Any]]:
        """
        Defines the biomass reaction by reading a csv file with the reactants of biomass and the respective coefficients
        ,mmol/gDW
            e-Protein,0.400
            e-Lipids,0.08
                ...
        Parameters
        ----------
        atp_coef: float
            coefficient for the ATP consumption in the biomass reaction

        Returns
        -------
        self.biomass_reactions: dict
            details abound the biomass reactions
        """
        biomass_path = os.path.join(self.project_path, 'reconstruction_results', self.model_id, 'biomass')

        if not os.path.isdir(biomass_path):
            logging.info('There is no biomass folder')
            return None

        biomass_reac_id = 'e-Biomass_' + self.model_id

        biomass_file = os.path.join(biomass_path, 'Biomass_reaction.csv')
        biomass_df = pd.read_csv(biomass_file, index_col=0)
        try:
            biomass_comps = biomass_df.to_dict()["mmol/gDW"]
        except KeyError:
            col = biomass_df.columns[0]
            biomass_comps = biomass_df.to_dict()[col]

        biomass_prods = {"e-Biomass": 1, "ADP": atp_coef, "Pi": atp_coef}

        self.biomass_reactions[biomass_reac_id] = {'reactants': biomass_comps, 'products': biomass_prods,
                                                   'direction': 'IRREVERSIBLE', 'lower_bound': 0,
                                                   'upper_bound': 100000, 'reaction_type': 'biomass',
                                                   'in_pathway': ['biomass_reactions'], 'compartment': ['cytosol'],
                                                   'models': {self.model_id: [[biomass_reac_id + '__cyto',
                                                                               'cytosol']]}}
        if biomass_reac_id not in self.biomass_reactions:
            self.reactions.append(biomass_reac_id)

        return self.biomass_reactions

    def calculate_biomass_from_genome(self, expr_file: str = None):
        """
        Estimate biomass composition in nucleotides and amino acids from DNA, RNA and protein sequences.
        If expression data is available, it can also be used to estimate aa percentage in the sequence.
        This function requires three input sequence files: protein.faa, dna.fna, rna.fna
        It creates the .csv files for protein, dna and rna reactions with the respective reactants and products with
        biocyc identifiers: Protein_reaction.csv, DNA_reaction.csv, RNA_reaction.csv
        Parameters
        ---------
        expr_file: str
            name of the transcriptomics file to use to calculate biomass composition
        """
        macros = ['Protein', 'DNA', 'RNA']

        protein_fasta = os.path.join(self.genome_path, 'protein.faa')
        dna_fasta = os.path.join(self.genome_path, 'dna.fna')
        rna_fasta = os.path.join(self.genome_path, 'rna.fna')

        biomass_path = os.path.join(self.project_path, 'reconstruction_results', self.model_id, 'biomass')
        if expr_file:
            expr_file_path = os.path.join(biomass_path, expr_file)
        else:
            expr_file_path = None

        biomass = BiomassFromGenome(protein_fasta=protein_fasta, dna_fasta=dna_fasta, rna_fasta=rna_fasta,
                                    expr_file=expr_file_path)

        for macro in macros:
            biomass.get_final_composition(molecule=macro, biomass_path=biomass_path)

    def create_macro_reaction(self, macro_csv: str, macro: str) -> Dict[str, Any]:
        """
        Defines a macro reaction (e.g.: e-Protein) by reading a csv file with the reactants and products and the
        respective coefficients
        Parameters
        ----------
        macro_csv: str
            name of the csv file with the metabolites of the macro reaction and the following structure:
            ,reactants,products
            MAL,0.2,
            UDP,,0.1
        macro: str
            name of the macro to create (e.g. protein, carbohydrates, DNA, etc...)
        Returns
        -------
        self.biomass_reactions: dict
            details abound the biomass reactions
        """
        macro_path = os.path.join(self.project_path, 'reconstruction_results', self.model_id, 'biomass')
        macro_file = os.path.join(macro_path, macro_csv)
        macro_df = pd.read_csv(macro_file, index_col=0)
        reacts = macro_df["reactants"][macro_df["reactants"].notnull()].to_dict()
        prods = macro_df["products"][macro_df["products"].notnull()].to_dict()

        macro_reac_id = 'e-' + macro + '_' + self.model_id

        self.biomass_reactions[macro_reac_id] = {'reactants': reacts, 'products': prods, 'direction': 'IRREVERSIBLE',
                                                 'lower_bound': 0, 'upper_bound': 100000, 'reaction_type': 'biomass',
                                                 'in_pathway': ['biomass_reactions'], 'compartment': ['cytosol'],
                                                 'models': {self.model_id: [[macro_reac_id + '__cyto',
                                                                             'cytosol']]}}
        if macro_reac_id not in self.reactions:
            self.reactions.append(macro_reac_id)

        return self.biomass_reactions

    def create_maintenance_reaction(self, bounds: tuple = (2, 2)) -> Dict[str, Any]:
        """
        Create the ATP maintenance reaction
        Parameters
        ----------
        bounds: tuple
            lower and upper bounds for the reaction
        Returns
        -------
        self.biomass_reactions: dict
            details abound the biomass reactions
        """
        reac_id = 'Maintenance_' + str(self.model_id)

        self.biomass_reactions[reac_id] = {'reactants': {"ATP": 1, "WATER": 1}, 'products': {"ADP": 1, "Pi": 1},
                                           'direction': 'IRREVERSIBLE', 'lower_bound': bounds[0],
                                           'upper_bound': bounds[1], 'reaction_type': 'biomass',
                                           'in_pathway': ['biomass_reactions'], 'compartment': ['cytosol'],
                                           'models': {self.model_id: [[reac_id + '__cyto', 'cytosol']]}}
        if reac_id not in self.reactions:
            self.reactions.append(reac_id)
        return self.biomass_reactions

    def create_all_biomass_reactions(self, atp_bounds: tuple = None):
        """
        Create the biomass reaction and all the associated reactions (for macromolecules and maintenance)
        Parameters
        ----------
        atp_bounds: tuple
            bounds for the atp maintenance reaction
        Returns
        -------
        """

        macros = ['DNA', 'RNA', 'Protein', 'Carbohydrates', 'Cellwall', 'Lipids', 'Cofactor', 'FattyAcid']

        list_of_csv = []
        for macro in macros:
            reac_csv = macro + '_reaction.csv'
            list_of_csv.append(reac_csv)

        self.create_biomass_reaction()

        for i in range(len(macros)):
            self.create_macro_reaction(macro_csv=list_of_csv[i], macro=macros[i])

        if atp_bounds:
            self.create_maintenance_reaction(bounds=atp_bounds)
        else:
            self.create_maintenance_reaction()

        filename = os.path.join(PROJECT_PATH, 'reconstruction_results', self.model_id, 'biomass.json')
        f = open(filename, 'w')

        json.dump(self.biomass_reactions, f)

    def create_drains(self, metabolites: Dict[str, tuple]) -> Dict[str, Any]:
        """
        Create drains for the metabolites in the list_of_metabolites

        Parameters
        ----------
        metabolites: dict
            metabolites in the model that need a drain and the respective bounds for the drain
            metabolite: (lower_bound, upper_bound)
        Returns
        -------
        self.drains: list
            drain in the model
        """
        drains = {}
        for met in metabolites:
            drain_id = 'EX_' + met + '_drain'
            if metabolites[met][0] == 0:
                direction = 'IRREVERSIBLE'
            else:
                direction = 'REVERSIBLE'
            drains[drain_id] = {'reactants': {met: 1}, 'products': {}, 'lower_bound': metabolites[met][0],
                                'upper_bound': metabolites[met][1], 'direction': direction, 'in_pathway': ['drains'],
                                'compartment': ['extracellular'], 'reaction_type': 'drain',
                                'models': {self.model_id: [[drain_id, 'extracellular']]}}
            if drain_id not in self.reactions:
                self.reactions.append(drain_id)
        self.drains = drains

        return self.drains

    def create_main_drains(self):
        """
        Create the main drains for metabolic models using metacyc identifiers for metabolites and generic bounds
        """
        reversible_drains = ['FE+2', 'MG+2', 'WATER', 'OXYGEN-MOLECULE', 'PROTON', 'Pi',
                             'CARBON-DIOXIDE', 'NITRATE', 'SULFATE']
        met_bounds = {}
        for met in reversible_drains:
            met_bounds[met] = (-100000, 100000)

        met_bounds['light'] = (-100, 100000)

        irreversible_drains = ['SUCROSE', 'ALPHA-GLUCOSE', 'BETA-D-FRUCTOSE', 'MALTOSE', 'L-ALPHA-ALANINE',
                               'L-ASPARTATE', 'ASN', 'ARG', 'CYS', 'HIS', 'GLY', 'GLN', 'GLT', 'LEU', 'ILE', 'LYS',
                               'MET', 'PRO', 'PHE', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

        for met in irreversible_drains:
            met_bounds[met] = (0, 100000)

        self.create_drains(metabolites=met_bounds)

        filename = os.path.join(PROJECT_PATH, 'reconstruction_results', self.model_id, 'drains.json')
        f = open(filename, 'w')

        json.dump(self.drains, f)

    def add_biomass_db(self):
        """
        Add biomass reactions to the databases
        """
        reaction_path = os.path.join(PROJECT_PATH, 'reconstruction_results', self.model_id, 'biomass.json')

        try:
            reaction_file = open(reaction_path, 'rb')
        except FileNotFoundError:
            logging.info('no file was found: ' + reaction_path)
            return None

        url_neo = self.api.base_url + self.api.add_new_reactions_neo + self.model_id

        url_mongo = self.api.base_url + self.api.add_new_reactions_mongo + self.model_id

        r_mongo = requests.post(url_mongo, files={'file': reaction_file})

        if r_mongo.status_code == 200:
            logging.info('new reactions were added to mongo database')
        else:
            logging.info('An error has occurred when connecting to the mongodb API')

        r_neo = requests.post(url_neo, files={'file': reaction_file})

        if r_neo.status_code == 200:
            logging.info('new reactions were added to neo4j database')
        else:
            logging.info('An error has occurred when connecting to the neo4j API')

    def add_drains_db(self):
        """
        Add drain reactions to the databases
        """
        reaction_path = os.path.join(PROJECT_PATH, 'reconstruction_results', self.model_id, 'drains.json')

        try:
            reaction_file = open(reaction_path, 'rb')
        except FileNotFoundError:
            logging.info('no file was found: ' + reaction_path)
            return None

        url_neo = self.api.base_url + self.api.add_new_reactions_neo + self.model_id

        url_mongo = self.api.base_url + self.api.add_new_reactions_mongo + self.model_id

        r_mongo = requests.post(url_mongo, files={'file': reaction_file})

        if r_mongo.status_code == 200:
            logging.info('new reactions were added to mongo database')
        else:
            logging.info('An error has occurred when connecting to the mongodb API')

        r_neo = requests.post(url_neo, files={'file': reaction_file})

        if r_neo.status_code == 200:
            logging.info('new reactions were added to neo4j database')
        else:
            logging.info('An error has occurred when connecting to the neo4j API')

    def create_gprs(self) -> List[Tuple[str, str]]:
        """
        Create gprs for the reactions in the model using database information about enzymatic complexes.
        Returns
        -------
        self.gprs: list of tuples
            ( reaction id , gpr rule )
        """
        gprs_dic = {}

        if not self.reactions:
            url = self.api.base_url + self.api.model_detail + self.model_id
            r = requests.get(url)

            if r.status_code == 200:
                res = r.json()
                self.reactions = res['reactions']
                self.diamond_annotation = res['annotation']
            else:
                logging.info('An error has occurred when connecting to the mongo API')

        for reaction in self.reactions:
            url = self.api.base_url + self.api.enzymes_from_reaction_mongo + reaction

            r = requests.get(url)

            if r.status_code != 200:
                logging.info('An error has occurred when connecting to the neo4j mongo for reaction ' + str(reaction))
                continue

            res = r.json()
            if res['enzymes']:
                all_prots = []

                for enz in res['enzymes']:

                    url_components = self.api.base_url + self.api.components_from_enzyme + enz

                    r = requests.get(url_components)

                    if r.status_code != 200:
                        logging.info('An error has occurred when connecting to the neo4j API for reaction ' + str(enz))
                        continue

                    res_enzyme = r.json()

                    if res_enzyme['components']:
                        if enz in self.diamond_annotation.values():
                            matches = [k for k, v in self.diamond_annotation.items() if v == enz]
                            if matches:
                                matches_str = ' or '.join(matches)
                                all_prots.append(matches_str)

                        else:
                            comps = []
                            for comp in res_enzyme['components']:
                                matches = [k for k, v in self.diamond_annotation.items() if
                                           v == comp]
                                if matches:
                                    if len(matches) == 1:
                                        matches_str = matches[0]
                                    else:
                                        matches_str = '(' + ' or '.join(matches) + ')'
                                    comps.append(matches_str)

                            if comps:
                                cplx_rule = '(' + ' and '.join(enz_comp for enz_comp in comps) + ')'
                                all_prots.append(cplx_rule)

                    else:
                        matches = [k for k, v in self.diamond_annotation.items() if v == enz]
                        if matches:
                            rule = ' or '.join(m for m in matches)
                            all_prots.append(rule)

                if all_prots:
                    final_rule = ' or '.join(prot for prot in all_prots)
                    if 'and' not in final_rule or 'or' not in final_rule:
                        final_rule = final_rule.replace('(', '').replace(')', '')
                    gprs_dic[reaction] = final_rule

        gprs_list = [(x, y) for x, y in gprs_dic.items()]

        self.gprs = gprs_list

        df = pd.DataFrame(gprs_list, columns=['reaction', 'rule'])
        df.to_csv(os.path.join(PROJECT_PATH, 'reconstruction_results', self.model_id, 'gprs.csv'), index=False)

        return self.gprs

    def add_gprs_to_model_db(self):
        """
        Add gprs to the metabolic model object in monogdb
        """
        gprs_path = os.path.join(PROJECT_PATH, 'reconstruction_results', self.model_id, 'gprs.csv')
        try:
            gprs_file = open(gprs_path, 'r')
        except FileNotFoundError:
            logging.info('No gprs file was found')
            return None

        url = self.api.base_url + self.api.add_gprs_model + self.model_id

        r = requests.post(url, files={'file': gprs_file})

        if r.status_code == 200:
            logging.info('gprs data was added to the model document')
        else:
            logging.info('An error has occurred when connecting to mongo API')

    def get_compartments(self, tool: str = 'loctree3') -> Dict[str, List[str]]:
        """
        Get the compartment prediction for each protein in the genome and map the compartments to the reactions
        Parameters
        ----------
        tool: str
            tool used for prediction protein subcellular localization. It only supports loctree3 or wolfpsort
        Returns
        -------
        reac_compartments: dict
            {reac_id: list of compartments}
        """
        res_folder = os.path.join(self.project_path, 'reconstruction_results', self.model_id, 'compartments')
        comps = Compartments(results_folder=res_folder)

        if tool == 'loctree3':
            pred_compartments = comps.parse_results_loctree3()
        else:
            pred_compartments = comps.parse_results_wolfpsort()

        self.prot_compartments = pred_compartments

        reac_compartments = {}

        if not self.gprs:
            url = self.api.base_url + self.api.model_detail + self.model_id
            r = requests.get(url)

            if r.status_code == 200:
                res = r.json()
                self.gprs = res['gprs']
                self.reactions = res['reactions']
            else:
                logging.info('An error has occurred when connecting to the mongo API')

        gprs_dic = {x: y for (x, y) in self.gprs}

        for reac in gprs_dic:
            reac_comparts = []
            rule = gprs_dic[reac].replace('(', '').replace(')', '')

            if 'or' in rule:
                or_objs = rule.split(' or ')
                and_objs = [x.split(' and ') if 'and' in x else x for x in or_objs]
                rule_list = []
                for elem in and_objs:
                    if isinstance(elem, list):
                        rule_list.extend(elem)
                    else:
                        rule_list.append(elem)

            elif 'and' in rule:
                rule_list = rule.split(' and ')
            else:
                rule_list = [rule]

            for protein in rule_list:

                try:
                    pred_comps = pred_compartments[protein]

                except KeyError:
                    pred_comps = 'cytosol'

                if isinstance(pred_comps, str):
                    if pred_comps not in reac_comparts:
                        reac_comparts.append(pred_comps)
                else:
                    for comp in pred_comps:
                        if comp not in reac_comparts:
                            reac_comparts.append(comp)

            reac_compartments[reac] = reac_comparts

        filename = os.path.join(PROJECT_PATH, 'reconstruction_results', self.model_id, 'compartments.json')
        f = open(filename, 'w')
        json.dump(reac_compartments, f)

        return reac_compartments

    def add_comparts_to_model_db(self):
        """
        Add compartment info to database
        """
        comparts_path = os.path.join(PROJECT_PATH, 'reconstruction_results', self.model_id, 'compartments.json')
        try:
            comparts_file = open(comparts_path, 'rb')
        except FileNotFoundError:
            logging.info('No comparts file was found')
            return None

        url = self.api.base_url + self.api.add_comparts_model + self.model_id

        r = requests.post(url, files={'file': comparts_file})

        if r.status_code == 200:
            logging.info('compartment data was added to the model document')
        else:
            logging.info('An error has occurred when connecting to mongo API')

    def run_transyt(self, **kwargs) -> Dict[str, Dict]:
        """
        Run transyt, parse the results and save them to the database
        Parameters
        ----------
        kwargs: dict
            keyword arguments representing the advanced options of transyt. The parameters not passed here will get
            the default value.
        Returns
        -------
        transporters: dict
            {transporter_id: dict of attributes}
        """
        prot_file = os.path.join(self.genome_path, "protein.faa")
        draft_model = os.path.join(self.project_path, 'reconstruction_results', self.model_id,
                                   self.model_id + '_draft.xml')
        output_path = os.path.join(self.project_path, 'reconstruction_results', self.model_id, 'transporters')

        transp = Transyt(modelid=self.model_id, taxid=self.taxid, protein_file=prot_file, draft_model=draft_model,
                         output_folder=output_path)

        transp.run_transyt(**kwargs)

        transporters = transp.parse_results(comparts_dic=self.prot_compartments)

        return transporters

    def add_transporters(self):
        """
        Add transport reactions to the databases
        """
        transporters_path = os.path.join(PROJECT_PATH, 'reconstruction_results', self.model_id, 'transporters.json')

        try:
            transporters_file = open(transporters_path, 'r')
        except FileNotFoundError:
            logging.info('No transporters file was found')
            return None

        url_mongo = self.api.base_url + self.api.add_new_reactions_mongo + self.model_id

        r_mongo = requests.post(url_mongo, files={'file': transporters_file})

        if r_mongo.status_code == 200:
            logging.info('transporters were added to the mongo database')
        else:
            logging.info('An error has occurred when connecting to mongo API')
            return None

        url_neo = self.api.base_url + self.api.add_new_reactions_neo + self.model_id

        r_neo = requests.post(url_neo, files={'file': transporters_file})

        if r_neo.status_code == 200:
            logging.info('transporters were added to the neo database')
        else:
            logging.info('An error has occurred when connecting to neo4j API')

        gprs_path = os.path.join(PROJECT_PATH, 'reconstruction_results', self.model_id, 'gprs_transporters.csv')

        try:
            gprs_file = open(gprs_path, 'r')
        except FileNotFoundError:
            logging.info('No gprs file was found')
            return None

        url = self.api.base_url + self.api.add_gprs_model + self.model_id

        r = requests.post(url, files={'file': gprs_file})

        if r.status_code == 200:
            logging.info('gprs data was added to the model document')
        else:
            logging.info('An error has occurred when connecting to mongo API')

    def export_sbml(self):
        """
        Export sbml model using cobra with data only from the mongo database.
        Returns
        -------
        cobra_model: Model
            instance of the cobra model
        """
        url = self.api.base_url + self.api.model_detail + self.model_id
        r = requests.get(url)

        if r.status_code == 200:
            res = r.json()
            self.metabolites = res['metabolites']
            self.reactions = set(res['reactions'])
            self.pathways = set(res['pathways'])
            self.gprs = res['gprs']
        else:
            logging.info('An error has occurred when connecting to the mongo API')
            return None

        gprs_dic = {x: y for (x, y) in self.gprs}

        cobra_model = CobraModel(self.model_id)

        compartments = {'C001': 'cytosol', 'C002': 'pmf', 'C003': 'chloroplast', 'C004': 'mitochondrion',
                        'C005': 'endoplasmatic reticulum', 'C006': 'peroxisome', 'C007': 'golgi',
                        'C008': 'vacuole', 'C009': 'extracellular', 'C010': 'lumen'}

        cobra_model.compartments = compartments

        path_objs = []
        for path in self.pathways:

            url_path = self.api.base_url + self.api.pathway_detail + path

            r_path = requests.get(url_path)

            if r_path.status_code == 200:
                res = r_path.json()
            else:
                continue

            path_obj = CobraGroup(id=path)
            if 'common_name' in res:
                path_obj.name = res['common_name']

            path_objs.append(path_obj)

        transp = CobraGroup(id='transporters', name='transport reactions')
        drains = CobraGroup(id='drains', name='drain reactions')
        biomass = CobraGroup(id='biomass', name='biomass reactions')

        path_objs.extend([transp, drains, biomass])

        cobra_model.add_groups(path_objs)

        comps_reversed = {v: k for k, v in compartments.items()}

        cobra_mets = []

        for met in self.metabolites:

            url_met = self.api.base_url + self.api.metabolite_detail + met

            r_met = requests.get(url_met)

            if r_met.status_code == 200:
                res_met = r_met.json()
            else:
                continue

            for met_comp in res_met['models'][self.model_id]:
                cobra_met = CobraMetabolite(met_comp)
                if 'common_name' in res_met:
                    cobra_met.name = res_met['common_name']
                else:
                    cobra_met.name = met
                if 'formula' in res_met:
                    cobra_met.formula = res_met['formula']
                compr_name = res_met['models'][self.model_id][met_comp]
                cobra_met.compartment = comps_reversed[compr_name]

                cobra_mets.append(cobra_met)

        cobra_model.add_metabolites(cobra_mets)

        cobra_reacs = []

        self.reactions = list(set(self.reactions))

        for reac in self.reactions:
            url_reac = self.api.base_url + self.api.reaction_detail + reac

            r_reac = requests.get(url_reac)

            if r_reac.status_code == 200:
                res_reac = r_reac.json()
            else:
                continue

            reac_comparts = {k: v for k, v in res_reac['models'][self.model_id]}

            for reac_comp in reac_comparts:
                if 'drain' in reac_comp:
                    suffix = 'extr'
                else:
                    suffix = reac_comp.split('__')[1]
                cobra_reac = CobraReaction(reac_comp)

                if 'common_name' in res_reac:
                    cobra_reac.name = res_reac['common_name']
                else:
                    cobra_reac.name = reac

                cobra_reac.lower_bound = res_reac['lower_bound']
                cobra_reac.upper_bound = res_reac['upper_bound']

                if 'ecnumber' in res_reac:
                    cobra_reac.annotation['ec-code'] = res_reac['ecnumber']

                reac_mets = {}

                for cmp in res_reac['reactants']:
                    if res_reac['reaction_type'] == 'transporter':
                        cmp_id = cmp
                    else:
                        cmp_id = cmp + '__' + suffix

                    try:
                        cmp_obj = cobra_model.metabolites.get_by_id(cmp_id)
                    except KeyError:
                        logging.info(cmp_id + ' not in the model')
                        continue

                    try:
                        reac_mets[cmp_obj] = float(res_reac['reactants'][cmp]) * -1
                    except ValueError:
                        reac_mets[cmp_obj] = -1

                for prd in res_reac['products']:
                    if res_reac['reaction_type'] == 'transporter':
                        prd_id = prd
                    else:
                        prd_id = prd + '__' + suffix

                    try:
                        prd_obj = cobra_model.metabolites.get_by_id(prd_id)
                    except KeyError:
                        logging.info(prd_id + ' not in the model')
                        continue

                    try:
                        reac_mets[prd_obj] = float(res_reac['products'][prd])
                    except ValueError:
                        reac_mets[prd_obj] = 1

                cobra_reac.add_metabolites(reac_mets)

                try:
                    rule = gprs_dic[reac]
                    cobra_reac.gene_reaction_rule = rule
                except KeyError:
                    pass

                cobra_reacs.append(cobra_reac)

                if 'in_pathway' in res_reac:
                    for path in res_reac['in_pathway']:
                        try:
                            path_obj = cobra_model.groups.get_by_id(path)
                            path_obj.add_members([cobra_reac])
                        except KeyError:
                            path_obj = CobraGroup(id=path)
                            path_obj.add_members([cobra_reac])
                            cobra_model.add_groups([path_obj])

        new_reacs = [r for r in cobra_reacs if r not in cobra_model.reactions]
        cobra_model.add_reactions(new_reacs)

        cobra_model.objective = 'e-Biomass_' + self.model_id + '__cyto'
        output_file = os.path.join(self.project_path, "reconstruction_results", self.model_id,
                                   self.model_id + "_draft.xml")
        write_sbml_model(cobra_model, output_file)

        return cobra_model

    def check_mass_balance(self, reacid: str) -> Optional[Dict[str, float]]:
        """
        Auxiliar function to calculate the mass balance of a reaction
        Parameters
        ----------
        reacid: str
            identifier of the reaction to evaluate
        Returns
        -------
        final_dic: dict
            dict of {atom: amount} for unbalaced atoms in a reaction.
        It is empty for balanced reactions
        """

        url_reac = self.api.base_url + self.api.reaction_detail + reacid

        r_reac = requests.get(url_reac)

        if r_reac.status_code == 200:
            res_reac = r_reac.json()
        else:
            logging.info('an error has occurred when connecting to mongo API')
            return None

        mets = {}
        if res_reac['reactants']:
            for m in res_reac['reactants']:
                try:
                    mets[m] = -1 * float(res_reac['reactants'][m])
                except ValueError:
                    mets[m] = -1

            for p in res_reac['products']:
                try:
                    mets[p] = float(res_reac['products'][p])
                except ValueError:
                    mets[p] = -1

        mass_balance = defaultdict(float)
        atom_re = re.compile("([A-Z][a-z]?)([0-9.]+[0-9.]?|(?=[A-Z])?)")

        for met, coeff in mets.items():

            url_met = self.api.base_url + self.api.metabolite_detail + met

            r_met = requests.get(url_met)

            if r_reac.status_code == 200:
                res_met = r_met.json()
            else:
                logging.info('an error has occurred when connecting to mongo API')
                return None

            if res_met['formula']:
                atom_composition = {}
                parsed = atom_re.findall(res_met['formula'])
                for (atom, count) in parsed:
                    if count == '':
                        count = 1
                    else:
                        try:
                            count = float(count)

                        except ValueError:
                            logging.info('Failed to parse ' + met + ' formula.')
                            return None

                    if atom in atom_composition:
                        atom_composition[atom] += count
                    else:
                        atom_composition[atom] = count

                for elem, amount in atom_composition.items():
                    mass_balance[elem] += coeff * amount

            else:
                logging.info('Error! Metabolite ' + met + ' has no formula defined in the database')

        final_dic = {k: v for k, v in mass_balance.items() if v != 0}

        return final_dic

