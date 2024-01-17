import os
import logging
import time
import json
import zipfile
from io import BytesIO
from typing import Dict, Optional
import requests
from Bio import SeqIO
from Bio.Blast import NCBIXML

logging.basicConfig(level=logging.DEBUG)


class Diamond:

    def __init__(self, source_folder, email, output_folder, run_type='diamondp'):
        """
        Class to run diamond
        Parameters
        ----------
        source_folder: str
            folder that contains the genomic sequence files for running diamond: protein.faa or genomic.fna, and
            SubmissionParameters.txt
        email: str
            email to perform the submission for diamond
        output_folder: str
            where to save diamond results
        run_type: str
            diamond to run. if type = 'diamondp' it runs Diamond blastp. if type = 'diamondx' it runs Diamond blastx
        """
        self.genome_path = source_folder
        self.output_folder = output_folder
        self.email = email
        self.run_type = run_type

    def run_diamond(self):
        """
        Runs diamondP or diamondX. Requires the protein.faa and the SubmissionParameters.txt files to run.
        It runs the diamond in a container on the server.
        """

        if self.run_type == 'diamondp':
            url_connection = "http://rosalind.di.uminho.pt:6082/iPlantsSubmission/DiamondP-iPlants/" + self.email

            protein_file = os.path.join(self.genome_path, 'protein.faa')

            if not os.path.isfile(protein_file):
                logging.info('Error! protein.faa file not found!')
                return None

            fasta = open(protein_file, 'rb')

        else:
            url_connection = "http://rosalind.di.uminho.pt:6082/iPlantsSubmission/DiamondX-iPlants/" + self.email

            genomic_file = os.path.join(self.genome_path, 'genomic.fna')

            if not os.path.isfile(genomic_file):
                logging.info('Error! genomic.fna file not found!')
                return None

            fasta = open(genomic_file, 'rb')

        settings_file = os.path.join(self.genome_path, 'SubmissionParameters.txt')

        if not os.path.isfile(settings_file):
            logging.info('A Settings file named SubmissionParameters.txt is needed to run diamond')
            return None

        settings = open(settings_file, 'rb')

        r = requests.post(url_connection, files=[('files', settings), ('files', fasta)])

        if r.status_code in [200, 201]:
            text = r.json()

            subid = str(text['submissionID'])

            time.sleep(10)

            status_url = "http://rosalind.di.uminho.pt:6082/status/" + self.email + "/" + subid

            status = requests.get(status_url)

            if status.status_code not in [200, 201, 202]:
                logging.info(status.text)
                return None

            while status.status_code != 200:
                status = requests.get(status_url)
                time.sleep(20)

            download_url = "http://rosalind.di.uminho.pt:6082/download/" + self.email + "/" + subid

            res = requests.get(download_url)

            zipf = zipfile.ZipFile(BytesIO(res.content))
            zipf.extractall(self.output_folder)

    def parse_results(self) -> Optional[Dict[str, str]]:
        """
        Parsing of the diamond results xml file. Creates a json file that associates each query identifier
        from the genome annotation to its match in the database (enzymes identifier).
        Returns
        -------
        results: dict
            annotation results {genomic identifier: database identifier}
        """
        res_file = os.path.join(self.output_folder, 'DiamondResults.xml')

        if self.run_type == 'diamondp':
            input_file = os.path.join(self.genome_path, 'protein.faa')
        else:
            input_file = os.path.join(self.genome_path, 'genomic.fna')

        if not os.path.isfile(res_file):
            logging.info('An error has occurred! Diamond results not found')
            return None

        recs = SeqIO.parse(input_file, 'fasta')
        ids = [rec.id for rec in recs]

        result_handle = open(res_file)
        blast_records = NCBIXML.parse(result_handle)
        results = {}
        for record in blast_records:
            key = record.query_id
            index = int(key.replace('Query_', '')) - 1
            query = ids[index]
            if record.alignments:
                for align in record.alignments:
                    match = align.title[:-1]
                    results[query] = match
            else:
                results[query] = 'no hits'

        outfile = os.path.join(self.output_folder, 'DiamondResults_parsed.json')

        with open(outfile, 'w') as json_file:
            json.dump(results, json_file)

        return results
