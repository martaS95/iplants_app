import pandas as pd
from Bio import SeqIO
import os
import logging

logging.basicConfig(level=logging.DEBUG)


class BiomassFromGenome:

    def __init__(self, protein_fasta, dna_fasta, rna_fasta, expr_file=None):

        """
        Class to obtain the composition from the genome of the biomass components: protein, dna e rna.
        Parameters
        ----------
        protein_fasta: str
            fasta file with the protein sequences of the genome
        dna_fasta: str
            fasta file with the dna sequence of the genome
        rna_fasta: str
            fasta file with the rna sequences of the genome
        expr_file: str
            csv file with transcriptomics data, following this structure:
            protein identifier (same identifier of the fasta file),expression
            XP_010650054.1,1.24
        """

        self.protein = protein_fasta
        self.dna = dna_fasta
        self.rna = rna_fasta
        self.expr = expr_file

        self.list_aa = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V",
                        "W", "Y"]

        self.list_dna = ["A", "C", "G", "T"]

        self.list_rna = ["A", "C", "G", "U"]

        self.mw_aa = {"A": 89.09, "R": 175.21, "N": 132.12, "D": 132.10, "C": 121.16, "Q": 146.15, "E": 146.12,
                      "G": 75.07, "H": 155.16, "I": 131.18, "L": 131.18, "K": 147.20, "M": 149.21, "F": 165.19,
                      "P": 115.13, "S": 105.09, "T": 119.12, "W": 204.23, "Y": 181.19, "V": 117.15}

        self.mw_dna = {"A": 491.2, "C": 467.2, "G": 507.2, "T": 482.2}

        self.mw_rna = {"A": 503.2, "C": 479.1, "G": 519.1, "U": 480.1}

        self.ids_aa = {"A": "Charged-ALA-tRNAs", "C": "Charged-CYS-tRNAs", "D": "Charged-ASP-tRNAs",
                       "E": "Charged-GLT-tRNAs", "F": "Charged-PHE-tRNAs", "G": "Charged-GLY-tRNAs",
                       "H": "Charged-HIS-tRNAs", "I": "Charged-ILE-tRNAs", "K": "Charged-LYS-tRNAs",
                       "L": "Charged-LEU-tRNAs", "M": "Charged-MET-tRNAs", "N": "Charged-ASN-tRNAs",
                       "P": "Charged-PRO-tRNAs", "Q": "Charged-GLN-tRNAs", "R": "Charged-ARG-tRNAs",
                       "S": "Charged-SER-tRNAs", "T": "Charged-THR-tRNAs", "V": "Charged-VAL-tRNAs",
                       "W": "Charged-TRP-tRNAs", "Y": "Charged-TYR-tRNAs"}

        self.ids_dna = {"A": "DATP", "T": "TTP", "C": "DCTP", "G": "DGTP"}

        self.ids_rna = {"A": "ATP", "U": "UTP", "C": "CTP", "G": "GTP"}

    def get_mmol_gDW(self, molecule: str, biomass_path: str) -> pd.Series:
        """
        Get mmol/gDW of each monomer in the sequence
        Parameters
        ----------
        molecule: str
            macromolecule that composes biomass: protein, DNA or RNA
        biomass_path: str
            folder where the biomass files are saved
        Returns
        -------
        mmol_gDW_series: Series
            mmol/gDW for each monomer
        """
        biomass_file = os.path.join(biomass_path, 'Biomass_reaction.csv')
        biomass_df = pd.read_csv(biomass_file, index_col=0)

        biomass_composition = {'Protein': biomass_df.loc['e-Protein', 'mmol/gDW'],
                               'DNA': biomass_df.loc['e-DNA', 'mmol/gDW'], 'RNA': biomass_df.loc['e-RNA', 'mmol/gDW']}

        if molecule == "Protein":
            list_monomers = self.list_aa
            mw = self.mw_aa
            sequence_file = self.protein
        elif molecule == "DNA":
            list_monomers = self.list_dna
            mw = self.mw_dna
            sequence_file = self.dna
        elif molecule == "RNA":
            list_monomers = self.list_rna
            mw = self.mw_rna
            sequence_file = self.rna
        else:
            raise ValueError('Molecule not valid! Choose one of these options: "Protein", "DNA" or "RNA"')

        if self.expr:
            content = get_aa_perc_from_genome_and_transc(protein_file=sequence_file, list_aa=self.list_aa,
                                                         expr_file=self.expr)
        else:
            content = get_perc_from_genome(sequence_file=sequence_file, monomers_list=list_monomers)

        mass_g = {}
        for elem in content:
            mass_g[elem] = content[elem] * mw[elem]
        total_mass_g = sum(mass_g.values())

        mass_g_perc = {}
        for elem in mass_g:
            mass_g_perc[elem] = mass_g[elem] / total_mass_g

        mmol_gDW = {}
        for elem in mass_g_perc:
            mmol_gDW[elem] = (mass_g_perc[elem] * biomass_composition[molecule] * 1000) / mw[elem]

        mmol_gDW_series = pd.Series(mmol_gDW)
        mmol_gDW_series.name = "mmol/gDW"
        mmol_gDW_series.to_csv(os.path.join(biomass_path, molecule + "_mmol/gDW_composition.csv"))

        return mmol_gDW_series

    def get_mmol_g(self, molecule: str, biomass_path: str) -> pd.Series:
        """
        Get mmol/g of each monomer in the sequence
        Parameters
        ----------
        molecule: str
            macromolecule that composes biomass: protein, DNA or RNA
        biomass_path: str
           folder where the biomass files are saved

        Returns
        -------
        mmol_g_series: Series
            mmol/g for each monomer
        """
        if molecule == "Protein":
            list_monomers = self.list_aa
            mw = self.mw_aa
            sequence_file = self.protein
        elif molecule == "DNA":
            list_monomers = self.list_dna
            mw = self.mw_dna
            sequence_file = self.dna
        elif molecule == "RNA":
            list_monomers = self.list_rna
            mw = self.mw_rna
            sequence_file = self.rna
        else:
            raise ValueError('Molecule not valid! Choose one of these options: "Protein", "DNA" or "RNA"')

        if self.expr:
            content = get_aa_perc_from_genome_and_transc(protein_file=sequence_file, list_aa=self.list_aa,
                                                         expr_file=self.expr)
        else:
            content = get_perc_from_genome(sequence_file=sequence_file, monomers_list=list_monomers)

        mass_g = {}
        for elem in content:
            mass_g[elem] = content[elem] * mw[elem]
        total_mass_g = sum(mass_g.values())

        mass_g_perc = {}
        for elem in mass_g:
            mass_g_perc[elem] = mass_g[elem] / total_mass_g

        mmol_g_gtotal = {}
        for elem in mass_g_perc:
            mmol_g_gtotal[elem] = (mass_g_perc[elem] * 1000) / mw[elem]

        mmol_g_series = pd.Series(mmol_g_gtotal)
        mmol_g_series.name = "mmol/g"
        mmol_g_series.to_csv(os.path.join(biomass_path, molecule + "_mmolg_composition.csv"))

        return mmol_g_series

    def get_final_composition(self, molecule: str, biomass_path: str) -> pd.DataFrame:
        """
        Get the final csv with the biocyc ids for all components of the molecule and the respective mmol/g-molecule.
        It also adds the products of the macro reaction
        Parameters
        ----------
        molecule: str
            macromolecule of biomass: Protein, DNA or RNA
        biomass_path: str
            folder where the biomass files are saved
        """
        comp_file = os.path.join(biomass_path, molecule + '_mmolg_composition.csv')

        if not os.path.isfile(comp_file):
            self.get_mmol_g(molecule=molecule, biomass_path=biomass_path)

        df_comp = pd.read_csv(comp_file, index_col=0)
        column_name = "mmol/g"

        new_inds = []
        if molecule == 'DNA':
            id_dic = self.ids_dna
            prods_ind = ['e-' + str(molecule), "PPI"]
            prod_values = [1, sum(df_comp[column_name])]
        elif molecule == "RNA":
            id_dic = self.ids_rna
            prods_ind = ['e-' + str(molecule), "PPI"]
            prod_values = [1, sum(df_comp[column_name])]
        elif molecule == 'Protein':
            id_dic = self.ids_aa
            prods_ind = []
            prod_values = []
            for aa in id_dic:
                prods_ind.append(id_dic[aa][8:])
                prod_values.append(df_comp.loc[aa, column_name])
            prods_ind.append('e-' + str(molecule))
            prods_ind.append('WATER')
            prod_values.append(1)
            prod_values.append(sum(df_comp[column_name]))
        else:
            raise ValueError('Molecule not valid! Choose one of these options: "Protein", "DNA" or "RNA"')

        for ind in df_comp.index:
            new_ind = id_dic[ind]
            new_inds.append(new_ind)

        new_df = df_comp.copy()
        new_df.index = new_inds

        prod_series = pd.Series(prod_values, index=prods_ind)

        final_df = pd.concat([new_df, prod_series])
        final_df.columns = ["reactants", "products"]

        final_df.to_csv(os.path.join(biomass_path, molecule + "_reaction.csv"))

        return final_df


def get_perc_from_genome(sequence_file: str, monomers_list: list) -> dict:
    """
    Auxiliar function to get the percentage of a monomer in a genome sequence (dna, rna or portein)
    Parameters
    ----------
    sequence_file: str
        fasta file with the genome sequence
    monomers_list: list
        monomers to count in the genome sequence
    Returns
    -------
    perc: dict
        percentage of each monomer in the sequence
    """
    counts = {}
    tamanho_total = 0

    records = SeqIO.parse(sequence_file, "fasta")
    for record in records:
        tamanho_total += len(record.seq)
        for elem in monomers_list:
            if elem != "U":
                elem_count = str(record.seq).upper().count(elem)
            else:
                elem_count = str(record.seq).upper().count("T")
            if elem not in counts:
                counts[elem] = elem_count
            else:
                counts[elem] += elem_count

    perc = {}
    for key in counts:
        perc[key] = counts[key] / tamanho_total

    return perc


def get_aa_perc_from_genome_and_transc(protein_file: str, list_aa: list, expr_file: str) -> dict:
    """
    Auxiliar function to get the percentage of an aa in a genome sequence using transcriptomics data
    Parameters
    ----------
    protein_file: str
        fasta file of the protein sequence
    list_aa: list
        amino acids in the sequence
    expr_file: str
        csv file with transcriptomics data with the following structure:
            protein identifier (same identifier of the fasta file),expression
            XP_010650054.1,1.24
    Returns
    -------
    perc: dict
        percentage of each amino acid in the sequence
    """
    aa_counts = {}

    gene_expr = pd.read_csv(expr_file, header=0, index_col=0)

    records = SeqIO.parse(protein_file, "fasta")
    for record in records:
        seq_id = record.id
        expression = gene_expr.loc[seq_id, 'expression']
        for aa in list_aa:
            aa_count_expr = str(record.seq).count(aa) * expression
            if aa in aa_counts:
                aa_counts[aa] += aa_count_expr
            else:
                aa_counts[aa] = aa_count_expr

    total = sum(aa_counts.values())

    aa_perc = {}
    for aa in aa_counts:
        aa_perc[aa] = aa_counts[aa] / total

    return aa_perc
