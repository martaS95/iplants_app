import os
from typing import Dict


class Compartments:

    def __init__(self, results_folder):
        """
        Parse the results of the prediction of the subcellular location of proteins
        Parameters
        ----------
        results_folder: str
            folder with the results
        """
        self.results_folder = results_folder
        self.loctree_res = {}
        self.wolfpsort_res = {}
        self.curated_res = {}

    def parse_results_wolfpsort(self, wolfpsort_output='wolfpsort_output.txt') -> Dict[str, str]:
        """
        Parse the results of wolfpsort. It only saves the first compartment predicted. The resulting compartments ids
        are converted to compartments names in the dictionary "convert_names".
        Parameters
        ----------
        wolfpsort_output: str
            file with the output from wolfpsort. It must be in the results folder
        Returns
        -------
        self.wolpsort_res: dict
            {protein_id: compartment predicted}
        """
        res_file = os.path.join(self.results_folder, wolfpsort_output)
        protein_dic = {}

        convert_names = {'plas': 'cytosol', 'cyto': 'cytosol', 'chlo': 'chloroplast', 'nucl': 'cytosol',
                         'extr': 'extracellular', 'pero': 'peroxisome', 'E.R.': 'endoplasmic reticulum',
                         'vacu': 'vacuole', 'mito': 'mitochondrion', 'golg': 'golgi', 'cysk': 'cytosol'}

        wolfpsort = open(res_file, 'r')

        wolfpsort.readline()
        line = wolfpsort.readline()

        while line:
            res = line.strip().split(' ')
            protein_id = res[0]
            comp = res[1]
            if '_' in comp:
                new_comp = comp.split('_')
                new_comp = [convert_names[x] for x in new_comp]
                if new_comp[0] == new_comp[1]:
                    new_comp = new_comp[0]
            else:
                new_comp = convert_names[comp]
            protein_dic[protein_id] = new_comp
            line = wolfpsort.readline()

        self.wolfpsort_res = protein_dic
        return self.wolfpsort_res

    def parse_results_loctree3(self, loctree_output='loctree3_output.lc3') -> Dict[str, str]:
        """
        Parse the results of loctree3. It only saves the first compartment predicted. The resulting compartments ids
        are converted to compartments names in the dictionary "convert_names".
        Parameters
        ----------
        loctree_output: str
            file with the output from loctree3. It must be in the results folder
        Returns
        -------
        self.loctree_res: dict
            {protein_id: compartment predicted}
        """
        convert_names = {'plasma': 'cytosol', 'nucleus': 'cytosol', 'cytoplasm': 'cytosol', 'plastid': 'chloroplast',
                         'secreted': 'extracellular', 'golgi apparatus': 'golgi',
                         'endoplasmic reticulum': 'endoplasmatic reticulum'}

        res_file = os.path.join(self.results_folder, loctree_output)

        loctree = open(res_file, 'r')
        loctree_dic = {}

        line = loctree.readline().strip()
        while line:
            if not line.startswith('#'):
                res = line.split('\t')
                protein_id = res[0]
                try:
                    comp = res[2]
                    comp = comp.replace(' membrane', '')
                    if comp in convert_names:
                        comp = convert_names[comp]
                except IndexError:
                    continue
                loctree_dic[protein_id] = comp
            line = loctree.readline().strip()

        self.loctree_res = loctree_dic
        return self.loctree_res


if __name__ == '__main__':
    folder = 'C:/Users/BiSBII/OneDrive - Universidade do Minho (1)/PycharmProjects/grapevine_model_update'
    comps = Compartments(results_folder=folder)
    comps.parse_results_wolfpsort(wolfpsort_output='wolfres.txt')
