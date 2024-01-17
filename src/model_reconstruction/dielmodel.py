import os
import logging

import cobra
import pandas as pd
from cobra.io import read_sbml_model, write_sbml_model
from cobra import Reaction as CobraReaction
from cobra import Metabolite as CobraMetabolite
from cobra import Model as CobraModel
from cobra.core import Group as CobraGroup
from src.utils.config import PROJECT_PATH
from cobra.flux_analysis import pfba


class DielModel:

    def __init__(self, model_id: str, model_file: str, storage_metabolites: dict):
        """
        Class implementing the reconstruction of a diel model (light vs dark metabolism)
        Parameters
        ----------
        model_id: str
            identifier for the metabolic model under reconstruction
        model_file: str
            name of the multitissue model file to use for the creation of the diel model
        storage_metabolites: dict
            dictionary with the name of the metabolites to exchange between light and dark phases and the bounds for
            each transporter (some transporters can be unidirectional)
        """

        self.model_id = model_id
        self.storage_metabolites = storage_metabolites

        self.multitissue_model_file = os.path.join(PROJECT_PATH, 'reconstruction_results', self.model_id, model_file)

        self.base_model = read_sbml_model(self.multitissue_model_file)

        self.new_model = CobraModel(self.model_id)
        self.new_model.compartments = self.base_model.compartments

    def replicate_reactions(self) -> cobra.Model:
        """
        Replicates the reactions in the model for Light and Dark conditions
        Returns
        -------
        self.new_model: cobra.Model
            model with the replicate reactions
        """
        new_reactions = []

        states = ['light', 'dark']

        not_replicate = ['TotalBiomass', 'EX_TotalBiomass_drain', 'EX_e-Biomass__cyto_berry_drain',
                         'EX_e-Biomass__cyto_stem_drain', 'EX_e-Biomass__cyto_leaf_drain']

        reactions_to_replicate = [r for r in self.base_model.reactions if r.id not in not_replicate]

        for reaction in reactions_to_replicate:
            for state in states:
                new_reac_obj = reaction.copy()
                new_reac_obj.id = reaction.id + '_' + state

                for met in reaction.metabolites:
                    new_met_id = met.id + '_' + state
                    try:
                        new_met_obj = self.base_model.metabolites.get_by_id(new_met_id)
                    except KeyError:
                        new_met_obj = met.copy()
                        new_met_obj.id = new_met_id

                    new_reac_obj.add_metabolites({met: 0,
                                                  new_met_obj: reaction.metabolites[met]}, combine=False)
                new_reactions.append(new_reac_obj)

        self.new_model.add_reactions(new_reactions)

        new_groups = []

        groups = self.base_model.groups
        for group in groups:
            new_group = CobraGroup(id=group.id, name=group.name)
            new_members = []
            for reaction in group.members:
                try:
                    r_dark = self.new_model.reactions.get_by_id(reaction.id + '_dark')
                    r_light = self.new_model.reactions.get_by_id(reaction.id + '_light')
                except KeyError:
                    continue

                new_members.append(r_dark)
                new_members.append(r_light)

            new_group.add_members(new_members)

            new_groups.append(new_group)

        self.new_model.add_groups(new_groups)

        return self.new_model

    def update_totalbiomass(self) -> cobra.Model:
        """
        Creates the sum of biomass reactions for each tissue (light + dark)
        Update the totalbiomass reaction to include the created metabolites (sum biomass for each tissue)
        Returns
        -------
        self.new_model: cobra.Model
            model with the totalbiomass reactions
        """
        biomasses = {}
        reactions_bio = self.new_model.groups.get_by_id('biomass').members
        for bio in reactions_bio:
            if bio.id.startswith('e-Biomass'):
                products = bio.products
                met_bio = [k for k in products if k.id.startswith('e-Biomass')][0]
                tissue = bio.id.split('_')[-2]
                if tissue not in biomasses:
                    biomasses[tissue] = [met_bio]
                else:
                    biomasses[tissue].append(met_bio)

        new_reactions = []
        new_prods = []
        for tissue in biomasses:
            new_reaction = CobraReaction(id=tissue + '_biomass', name='Sum ' + tissue + ' Biomass')
            new_product = CobraMetabolite(id=tissue + '_biomass', compartment='C011')
            new_reaction.add_metabolites({met: -1 for met in biomasses[tissue]})
            new_reaction.add_metabolites({new_product: 1})

            new_reactions.append(new_reaction)
            new_prods.append(new_product)

        self.new_model.add_reactions(new_reactions)

        totalbio_reaction = self.base_model.reactions.get_by_id('TotalBiomass')
        totalbio_reaction.add_metabolites({met: 0 for met in totalbio_reaction.reactants}, combine=False)
        totalbio_reaction.add_metabolites({met: -1 for met in new_prods}, combine=False)

        drain = self.base_model.reactions.get_by_id('EX_TotalBiomass_drain')

        drains_bio = [drain]

        for p in new_prods:
            drain_name = 'EX_' + p.id + '_drain'
            drain_reac = CobraReaction(drain_name)
            drain_reac.bounds = (0, 10000)
            drain_reac.add_metabolites({p: -1})
            drains_bio.append(drain_reac)

        self.new_model.add_reactions([totalbio_reaction] + drains_bio)
        self.new_model.objective = 'TotalBiomass'

        group_bio = self.new_model.groups.get_by_id('biomass')
        group_bio.add_members(new_reactions + [totalbio_reaction])
        group_drain = self.new_model.groups.get_by_id('drains')
        group_drain.add_members(drains_bio)

        return self.new_model

    def add_storage_metabolites(self) -> cobra.Model:
        """
        Add the reactions for the storage of metabolites between light and dark phases
        Parameters
        ----------
        Returns
        -------
        new_model: cobra.Model
            model with the storage reactions
        """

        metabolites = {met.id: met for met in self.new_model.metabolites}

        met_bounds = {}

        for comp, mets in self.storage_metabolites.items():
            for met in mets:
                met_bounds[met] = self.storage_metabolites[comp][met]

        storage_dic = {}

        for compartment, mets in self.storage_metabolites.items():
            for met in mets:
                met_objs = [value for key, value in metabolites.items() if key.startswith(met + '__' + compartment)]

                tissues = {}
                for met_obj in met_objs:
                    tissue = met_obj.id.split('_')[-2]
                    state = met_obj.id.split('_')[-1]
                    if tissue not in tissues:
                        tissues[tissue] = {state: met_obj}
                    else:
                        tissues[tissue][state] = met_obj

                storage_dic[met] = tissues

        storage_reactions = []

        for metabolite in storage_dic:
            for tissue in storage_dic[metabolite]:
                met_light = storage_dic[metabolite][tissue]['light']
                met_dark = storage_dic[metabolite][tissue]['dark']
                reac_id = met_light.id + '_dark_storage'
                new_reaction = CobraReaction(id=reac_id, name='Storage of ' + metabolite + ' in Light')
                new_reaction.bounds = met_bounds[metabolite]
                new_reaction.add_metabolites({met_light: -1, met_dark: 1})

                storage_reactions.append(new_reaction)

        self.new_model.add_reactions(storage_reactions)

        new_group = CobraGroup(id='Storage', name='Storage between light and dark phases')
        new_group.add_members(storage_reactions)
        self.new_model.add_groups([new_group])

        return self.new_model

    def create_diel(self):
        """
        Creates and writes the final diel model
        """
        logging.info('duplicating the reactions for light and dark phases')
        self.replicate_reactions()

        logging.info('adding the storage transporters')
        self.add_storage_metabolites()

        logging.info('update biomass reactions')
        self.update_totalbiomass()

        write_sbml_model(self.new_model, os.path.join(PROJECT_PATH, 'reconstruction_results', self.model_id,
                                                      self.multitissue_model_file.replace('.xml',
                                                                                          '_diel.xml')))

    def diel_model_simulation(self):
        """
        Simulates the diel model under photorespiration conditions and saves the fluxes in a csv file

        Returns
        -------
        df_simulation: pd.Dataframe
            simulation results
        """

        if self.new_model.reactions:
            final_model = self.new_model
        else:
            diel_model_file = os.path.join(PROJECT_PATH, 'reconstruction_results', self.model_id,
                                           self.multitissue_model_file.replace('.xml', '_diel.xml'))
            final_model = read_sbml_model(diel_model_file)

        final_model.reactions.get_by_id('EX_Light_drain_dark').lower_bound = 0
        final_model.reactions.get_by_id('EX_Light_drain_light').lower_bound = -300

        final_model.reactions.get_by_id('EX_SUCROSE_drain_light').lower_bound = 0
        final_model.reactions.get_by_id('EX_SUCROSE_drain_dark').lower_bound = 0

        final_model.reactions.get_by_id('EX_CARBON-DIOXIDE_drain_light').lower_bound = -10000
        final_model.reactions.get_by_id('EX_CARBON-DIOXIDE_drain_dark').lower_bound = 0

        final_model.reactions.get_by_id('RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN__chlo_leaf_light').bounds = (0, 10000)
        final_model.reactions.get_by_id('RXN-961__chlo_leaf_light').bounds = (0, 10000)
        final_model.reactions.get_by_id('GLUTAMATE-SYNTHASE-FERREDOXIN-RXN__chlo_leaf_light').bounds = (0, 10000)
        final_model.reactions.get_by_id('GLUTAMATE-SYNTHASE-NADH-RXN__chlo_leaf_light').bounds = (0, 0)

        final_model.reactions.get_by_id('T_PROTON__vacu_stem_light').bounds = (0, 0)
        final_model.reactions.get_by_id('T_PROTON__vacu_stem_dark').bounds = (0, 0)
        final_model.reactions.get_by_id('T_PROTON__vacu_leaf_light').bounds = (0, 0)
        final_model.reactions.get_by_id('T_PROTON__vacu_leaf_dark').bounds = (0, 0)
        final_model.reactions.get_by_id('T_PROTON__vacu_berry_light').bounds = (0, 0)
        final_model.reactions.get_by_id('T_PROTON__vacu_berry_dark').bounds = (0, 0)

        final_model.reactions.get_by_id('EX_e-Biomass_drain_dark').bounds = (0, 0)
        final_model.reactions.get_by_id('EX_e-Biomass_drain_light').bounds = (0, 0)
        final_model.reactions.get_by_id('EX_leaf_biomass_drain').bounds = (0, 0)
        final_model.reactions.get_by_id('EX_stem_biomass_drain').bounds = (0, 0)
        final_model.reactions.get_by_id('EX_berry_biomass_drain').bounds = (0, 0)

        # RUBISCO CONSTRAINT
        q = 3
        same_flux = final_model.problem.Constraint(final_model.reactions.get_by_id(
            'RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN__chlo_leaf_light').flux_expression -
            final_model.reactions.get_by_id('RXN-961__chlo_leaf_light').flux_expression * q, lb=0, ub=0)

        final_model.add_cons_vars(same_flux)

        # NITRATE UPTATE CONTRAINT
        q = 3 / 2
        same_flux_nitrate = final_model.problem.Constraint(final_model.reactions.get_by_id(
            'EX_NITRATE_drain_light').flux_expression - final_model.reactions.get_by_id(
            'EX_NITRATE_drain_dark').flux_expression * q, lb=0, ub=0)

        final_model.add_cons_vars(same_flux_nitrate)

        res = pfba(final_model)

        print(final_model.summary(solution=res))

        data = []
        for ind, row in res.fluxes.iteritems():
            if row >= 1e-5 or row <= -1e-5:
                reac = final_model.reactions.get_by_id(ind)
                groups = [(c.id, c.name) for c in final_model.get_associated_groups(reac)]
                line = [reac.id, reac.reaction, row, '|'.join(str(x) for x in groups)]
                data.append(line)

        df_simulation = pd.DataFrame(data, columns=['reaction id', 'equation', 'flux', 'pathway'])

        model_name = self.multitissue_model_file.split('_')[-1].replace('.xml', '')

        df_simulation.to_csv(os.path.join(PROJECT_PATH, 'reconstruction_results', self.model_id,
                                          'diel_simulations_maxbio_' + model_name + '.csv'), index=False)

        return df_simulation


if __name__ == '__main__':
    modelid = 'vvinif2023'
    model = 'vvinif2023_multissue_mature.xml'
    storage = {'chlo': {'Starch': (-10000, 10000)},
               'vacu': {'BETA-D-FRUCTOSE': (-10000, 10000),
                        'SUCROSE': (-10000, 10000),
                        'ALPHA-GLUCOSE': (-10000, 10000),
                        'MAL': (-10000, 10000),
                        'CIT': (-10000, 10000),
                        'NITRATE': (-10000, 10000)},
               'cyto': {'L-ALPHA-ALANINE': (0, 10000),
                        'L-ASPARTATE': (0, 10000),
                        'ASN': (0, 10000),
                        'ARG': (0, 10000),
                        'HIS': (0, 10000),
                        'GLY': (0, 10000),
                        'GLN': (0, 10000),
                        'GLT': (0, 10000),
                        'LEU': (0, 10000),
                        'ILE': (0, 10000),
                        'LYS': (0, 10000),
                        'MET': (0, 10000),
                        'PHE': (0, 10000),
                        'SER': (0, 10000),
                        'THR': (0, 10000),
                        'TRP': (0, 10000),
                        'TYR': (0, 10000),
                        'VAL': (0, 10000),
                        'PRO': (0, 10000),
                        'CYS': (0, 10000)}}

    dm = DielModel(model_id=modelid, model_file=model, storage_metabolites=storage)
    # dm.create_diel()
    dm.diel_model_simulation()
