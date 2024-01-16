import os.path
from cobra.io import read_sbml_model
import requests
import zipfile
from io import BytesIO
import logging
import time
import datetime
import json
from typing import Dict
from src.utils.config import PROJECT_PATH
import pandas as pd
import itertools


logging.basicConfig(level=logging.DEBUG)


class Transyt:

    def __init__(self, modelid, taxid, protein_file, draft_model, output_folder):
        """
        Parse the results from transyt
        Parameters
        ----------
        modelid: str
            model identifier
        taxid: int
            taxonomy identifier of the organism
        protein_file: str
            fasta file with all the proteins in the genome
        draft_model: str
            draft model in xml for the organism
        output_folder: str
            folder to save transyt results
        """
        self.modelid = modelid
        self.taxid = taxid
        self.protein = protein_file
        self.model = draft_model
        self.output_folder = output_folder
        # self.res_folder = None
        self.res_folder = 'results_20230227_1057'

    def run_transyt(self, **kwargs):
        """
        Runs transyt in a container on the server. It requires the fasta file with the proteins in the genome and a
        draft of the model under reconstruction in xml.
        Parameters
        ----------
        **kwargs: dict
            keyword arguments representing the advanced options of transyt. The parameters not passed here will get
            the default value.
        """

        default_values = {'dbs': 'MetaCyc', 'ontologies': 'False', 'alpha': '0.75', 'beta': '0.3', 'minHits': '2',
                          'bitScore': '50', 'queryCov': '80', 'blastEvalThresh': '1e-20', 'scoreThresh': '0.75',
                          'simScore': '30', 'alphaFam': '0.4', 'autoAcceptEval': '1e-0', 'perc': '10',
                          'limitAcceptEval': '1e-50', 'ignoreM2': 'False'}

        args_list = [str(self.taxid)]

        for key in default_values:
            if key in kwargs:
                value = str(kwargs[key])
            else:
                value = default_values[key]
            args_list.append(value)

        url = '/'.join(args_list)
        url_connection = "https://transyt.bio.di.uminho.pt/submitMerlinPlugin/" + url

        genome = open(self.protein, "rb")
        model = open(self.model, "rb")

        r_submit = requests.post(url_connection, files=[('file', genome), ('file', model)])

        if r_submit.status_code == 202 or r_submit.status_code == 201:
            text = r_submit.json()
            subid = str(text['submissionID'])
            status_url = "https://transyt.bio.di.uminho.pt/status/" + str(subid) + "/True"
            status = requests.get(status_url)

            if status.status_code not in [200, 201, 202]:
                logging.info(status.text)
                return None

            while status.status_code != 200:
                status = requests.get(status_url)
                time.sleep(120)

            download_url = "https://transyt.bio.di.uminho.pt/download/" + str(subid)

            results = requests.get(download_url)

            zipf = zipfile.ZipFile(BytesIO(results.content))
            zipf.extractall(self.output_folder)
            new_name = 'results_' + str(datetime.datetime.now().strftime("%Y%m%d_%I%M"))
            self.res_folder = new_name
            os.rename(os.path.join(self.output_folder, 'results'), os.path.join(self.output_folder, new_name))

            logging.info('transyt results were downloaded')

        else:
            logging.info(r_submit.text)

    def parse_results(self, comparts_dic) -> Dict[str, Dict]:
        """
        Parse the results from transyt
        Parameters
        ----------
        comparts_dic: dict
            compartment predictions for each protein in the genome

        Returns
        -------
        transporters: dict
            dict with the info of each transport reaction to load to mongo database
        """
        transyt_file = os.path.join(self.output_folder, self.res_folder, 'transyt.xml')
        transyt_model = read_sbml_model(transyt_file)

        compart_suffix = {'cytosol': 'cyto', 'pmf': 'pmf', 'chloroplast': 'chlo', 'mitochondrion': 'mito',
                          'endoplasmatic reticulum': 'er', 'peroxisome': 'pero', 'golgi': 'golgi', 'vacuole': 'vacu',
                          'extracellular': 'extr', 'secreted': 'extr', 'lumen': 'lum', 'endoplasmic reticulum': 'er'}

        transporter_mets = {}
        transporters = {}
        transporter_gprs = []

        for reaction in transyt_model.reactions:
            if reaction.lower_bound == 0:
                direction = 'IRREVERSIBLE'
            else:
                direction = 'REVERSIBLE'

            genes = reaction.genes
            rule = reaction.gene_reaction_rule

            comparts = []
            for gene in genes:
                try:
                    compart = comparts_dic[gene.id]
                except KeyError:
                    compart = 'cytosol'
                if isinstance(compart, str):
                    if compart not in comparts:
                        comparts.append(compart)
                else:
                    for alt in compart:
                        if alt not in comparts:
                            comparts.append(alt)

            for c in comparts:
                new_reac = reaction.id + '__' + compart_suffix[c]

                reactants = {}
                products = {}

                for met in reaction.metabolites:
                    if met.id.endswith('_e0'):
                        if c == 'cytosol':
                            new_suffix = '__' + compart_suffix['extracellular']
                            comp = 'extracellular'
                        else:
                            new_suffix = '__' + compart_suffix['cytosol']
                            comp = 'cytosol'

                        new_met = met.id.replace('_e0', new_suffix)
                        if new_suffix not in new_reac:
                            new_reac += new_suffix

                    else:
                        new_suffix = '__' + compart_suffix[c]
                        new_met = met.id.replace('_c0', new_suffix)
                        comp = c

                    if met.id[:-3] in transporter_mets:
                        if new_met not in transporter_mets[met.id[:-3]]:
                            transporter_mets[met.id[:-3]][new_met] = comp
                    else:
                        transporter_mets[met.id[:-3]] = {new_met: comp}

                    if met in reaction.reactants:
                        reactants[new_met] = reaction.metabolites[met] * (-1)
                    else:
                        products[new_met] = reaction.metabolites[met]

                transporters[new_reac] = {'reactants': reactants, 'products': products, 'direction': direction,
                                          'lower_bound': int(reaction.lower_bound),
                                          'upper_bound': int(reaction.upper_bound), 'reaction_type': 'transporter',
                                          'in_pathway': ['transporters'], 'compartment': [],
                                          'models': {self.modelid: [(new_reac, c)]}}

                transporter_gprs.append([new_reac, rule])

        transporters['metabolites'] = transporter_mets

        filename = os.path.join(PROJECT_PATH, 'reconstruction_results', self.modelid, 'transporters.json')
        f = open(filename, 'w')
        json.dump(transporters, f)

        df = pd.DataFrame(transporter_gprs, columns=['reaction', 'rule'])
        df.to_csv(os.path.join(PROJECT_PATH, 'reconstruction_results', self.modelid, 'gprs_transporters.csv'),
                  index=False)

        return transporters

    def get_tc_distribution(self):
        """
        Parse the summary file from transyt to get the TC numbers in the model
        """
        transyt_file = os.path.join(self.output_folder, self.res_folder, 'scoresMethod1.txt')
        ecs = {'1': [], '2': [], '3': [], '4': [], '5': []}

        fic = open(transyt_file, 'r')
        line = fic.readline()
        i = 0
        while line:
            if not line.startswith('>') and not line.startswith('\n'):
                ec = line[0]
                transps = line.strip().split('\t')[-1][1:-1]
                transp_list = transps.split(', ')
                for transp in transp_list:
                    if transp not in ecs[ec]:
                        ecs[ec].append(transp)
            line = fic.readline()

        react_list = [ecs[x] for x in ecs]

        print(len(set(list(itertools.chain.from_iterable(react_list)))))

        # with open(os.path.join(self.output_folder, 'tcs.json'), 'w') as fp:
        #     json.dump(ecs, fp)



if __name__ == '__main__':
    path = "C:/Users/BiSBII/OneDrive - Universidade do Minho (1)/PycharmProjects/grapevine_model/" \
             "grapevine_genome_files/protein.faa"
    drmodel = "C:/Users/BiSBII/Documents/plantdb/reconstruction_results/vvinif2021_draft.xml"
    output = os.path.join(PROJECT_PATH, 'reconstruction_results', 'vvinif2023', 'transporters')
    tr = Transyt(modelid='vvinif2023', protein_file=path, draft_model=drmodel, taxid=29760, output_folder=output)
    # tr.run_transyt()
    tr.get_tc_distribution()
