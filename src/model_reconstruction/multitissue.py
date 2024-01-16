import os
import logging
import cobra
from src.utils.config import PROJECT_PATH
from cobra.io import read_sbml_model, write_sbml_model
from cobra import Model as CobraModel
from cobra import Reaction as CobraReaction
from cobra import Metabolite as CobraMetabolite


logging.basicConfig(level=logging.DEBUG)


class Multitissue:

    def __init__(self, model_id: str, dataset_id: str, tissue_models: dict, base_tissue: str, tissue_pools: dict,
                 tissue_drains: dict, totalbiomass_composition: dict):
        """
        Class implementing the reconstruction of a multitissue model

        Parameters
        ----------
        model_id: str
            identifier for the metabolic model under reconstruction
        dataset_id: str
            identifier for the dataset used to build the tissue-specific models
        tissue_models: dict
            name of the model file for each tissue
        base_tissue: str
            tissue that performs the main exchanges with the environment
        tissue_pools: dict
            the common pool each tissue contacts with
        tissue_drains: dict
            drains to include in each tissue
        totalbiomass_composition: dict
            composition of the biomass of each tissue in the total biomass
        """
        self.model_id = model_id
        self.dataset_id = dataset_id

        path = os.path.join(PROJECT_PATH, 'reconstruction_results', self.model_id, 'results_troppo', self.dataset_id,
                            'reconstructed_models', 'TISSUE_MODELS_CHOSEN')

        self.tissue_models = {k: os.path.join(path, v) for k, v in tissue_models.items()}
        self.base_tissue = base_tissue
        self.tissue_pools = tissue_pools
        self.tissue_drains = tissue_drains
        self.totalbio_comp = totalbiomass_composition
        self.multitissue = CobraModel(self.model_id + '_multi')
        self.multitissue_reactions = []

    @staticmethod
    def update_mets_ids(model: cobra.Model, tissue: str, pools: list):
        """
        Method to update the identifiers of the metabolites.
        For internal metabolites, it adds '_tissue' to the id
        For external metabolites, it creates copies of the metabolite object for each pool
        Parameters
        ----------
        model: cobra.Model
            cobra model of a tissue
        tissue: str
            name of the tissue to add to the ids
        pools: list
            list of the pools in which the metabolites will be
        Returns
        -------
        model: cobra.Model
            tissue model with the modifications
        """
        extr_metabolites = [met for met in model.metabolites if '__extr' in met.id]
        other_metabolites = [met for met in model.metabolites if met not in extr_metabolites]

        for met in other_metabolites:
            new_metid = met.id + '_' + tissue
            met.id = new_metid

        new_metabolites = []
        for met in extr_metabolites:
            for pool in pools:
                new_met = met.copy()
                new_met.id = met.id.replace('extr', pool)
                new_met.compartment = pool
                new_metabolites.append(new_met)

        model.add_metabolites(new_metabolites)

        return model

    @staticmethod
    def update_internal_reactions(model: cobra.Model, internal_reactions: list, tissue: str):
        """
        Update the name of internal reactions by adding '_tissue' to the reaction identifier
        Parameters
        ----------
        model: cobra.Model
            cobra model of a tissue
        internal_reactions: list
            list of internal reactions in the model
        tissue: str
            name of the tissue
        Returns
        -------
        model: cobra.Model
            tissue model with the modifications
        """
        for reaction in internal_reactions:
            new_reacid = reaction.id + '_' + tissue
            reaction.id = new_reacid
        return model

    @staticmethod
    def update_drain_ids(model: cobra.Model, drains: list):
        """
        Update the metabolite identifiers of drains by replacing 'extr' with 'env'
        Parameters
        ----------
        model: cobra.Model
            cobra model of a tissue
        drains: list
            list of drains in the model
        Returns
        -------
        model: cobra.Model
            tissue model with the modifications
        """
        for drain in drains:
            for met in drain.metabolites:
                new_metid = met.id.replace('extr', 'env')
                try:
                    new_met_obj = model.metabolites.get_by_id(new_metid)
                except KeyError:
                    new_met_obj = met.copy()
                    new_met_obj.id = new_metid
                    new_met_obj.compartment = 'env'
                drain.add_metabolites({met: 0, new_met_obj: -1}, combine=False)
        return model

    @staticmethod
    def update_extr_reactions(model: cobra.Model, extr_reactions: list, pools_list: list):
        """
        Update the name of external reactions by adding them to all pools in pools_list
        Parameters
        ----------
        model: cobra.Model
            cobra model of a tissue
        extr_reactions: list
            list of external reactions in the model (reactions that occur in the external compartment)
        pools_list: list
            common pools where the reactions occur
        Returns
        -------
        model: cobra.Model
            tissue model with the modifications
        """
        new_reactions = []
        old_reactions = []
        for extra in extr_reactions:
            for cp in pools_list:
                new_reac_obj = extra.copy()
                new_reac_obj.id = new_reac_obj.id.replace('extr', cp)
                groups = model.get_associated_groups(extra)
                for group in groups:
                    group.add_members([new_reac_obj])
                for met in new_reac_obj.metabolites:
                    new_metid = met.id.replace('extr', cp)
                    new_met_obj = model.metabolites.get_by_id(new_metid)
                    new_reac_obj.add_metabolites({met: 0,
                                                  new_met_obj: new_reac_obj.metabolites[met]}, combine=False)
                new_reactions.append(new_reac_obj)
            old_reactions.append(extra)

        model.add_reactions(new_reactions)
        model.remove_reactions(old_reactions)

        return model

    def update_extr_transporters(self, model: CobraModel, extr_transporters: list, pools_list: list, tissue: str):
        """
        Update the transporters that transport external metabolites from/to the env and the common pools
        Parameters
        ----------
        model: cobra.Model
            cobra model of a tissue
        extr_transporters: list
            list of reactions that transport an external metabolite
        pools_list: list
            list of the pools to where the transport will occur
        tissue: str
            name of the tissue
        Returns
        -------
        model: cobra.Model
            tissue model with the modifications
        """
        new_transporters = []
        old_transporters = []

        for transporter in extr_transporters:

            if tissue != self.base_tissue and transporter.id in self.tissue_drains[tissue]:
                pools = pools_list + ['env']
            else:
                pools = pools_list
            for cp in pools:
                new_reac_obj = transporter.copy()
                if 'extr' in new_reac_obj.id:
                    new_reac_obj.id = new_reac_obj.id.replace('extr', cp) + '_' + tissue
                else:
                    new_reac_obj.id = new_reac_obj.id + '_' + cp + '_' + tissue

                for met in new_reac_obj.metabolites:
                    if 'extr' in met.id:
                        new_metid = met.id.replace('extr', cp)
                    else:
                        new_metid = met.id
                    new_met_obj = model.metabolites.get_by_id(new_metid)
                    new_reac_obj.add_metabolites({met: 0,
                                                  new_met_obj: new_reac_obj.metabolites[met]}, combine=False)
                new_transporters.append(new_reac_obj)
            old_transporters.append(transporter)

        model.add_reactions(new_transporters)
        model.remove_reactions(old_transporters)

        transp_group = model.groups.get_by_id('transporters')
        transp_group.add_members(new_transporters)

        return model

    @staticmethod
    def create_total_biomass_reaction(model: cobra.Model, metabolites: dict):
        """
        Create the total biomass reaction that adds the biomass reactions of all tissues. It also creates the drains
        for the total biomass metabolite as well as for the tissue's biomass metabolites.
        Parameters
        ----------
        model: cobra.Model:
            the multitissue model
        metabolites: dict
            the composition of each tissue's biomass in total biomass reaction
        Returns
        -------
        model: cobra.Model
            multitissue model with the modifications
        """
        new_reaction = CobraReaction(id='TotalBiomass')
        new_reaction.bounds = (0, 10000)
        mets_obj = {}
        reactions_to_add = []
        for m in metabolites:
            try:
                new_met = model.metabolites.get_by_id(m)
            except KeyError:
                new_met = CobraMetabolite(id=m, compartment='C011')

            mets_obj[new_met] = metabolites[m]

            drain_name = 'EX_' + m + '_drain'

            drain_reac = CobraReaction(drain_name)
            drain_reac.bounds = (0, 10000)
            drain_reac.add_metabolites({new_met: -1})
            reactions_to_add.append(drain_reac)

        new_reaction.add_metabolites(mets_obj)

        group = model.groups.get_by_id('drains')
        group.add_members(reactions_to_add)

        reactions_to_add.append(new_reaction)
        model.add_reactions(reactions_to_add)

        group = model.groups.get_by_id('biomass')
        group.add_members([new_reaction])

        return model

    def add_base_tissue(self):
        """
        First step in the reconstruction. Start with the model for the base tissue. The base tissue is the one that
        contacts more with the envirnoment and with the other tissues.
        Returns
        -------
        model_base: cobra.Model
            model for the base tissue
        """
        model_base = read_sbml_model(self.tissue_models[self.base_tissue])

        model_base = self.update_mets_ids(model=model_base, tissue=self.base_tissue,
                                          pools=self.tissue_pools[self.base_tissue] + ['env'])

        transporters = model_base.groups.get_by_id('transporters').members

        extr_transporters = [t for t in transporters if '__extr' in t.id or
                             (t.id.startswith('T_') and t.id.endswith('__cyto'))]

        drains = [d for d in model_base.groups.get_by_id('drains').members]

        extr_reacs = [r for r in model_base.reactions if '_extr' in r.id and r not in extr_transporters]

        all_excepts = extr_transporters + drains + extr_reacs

        internal = [reac for reac in model_base.reactions if reac not in all_excepts]

        model_base = self.update_internal_reactions(model=model_base, internal_reactions=internal,
                                                    tissue=self.base_tissue)

        model_base = self.update_drain_ids(model_base, drains=drains)

        model_base = self.update_extr_reactions(model_base, extr_reactions=extr_reacs,
                                                pools_list=self.tissue_pools[self.base_tissue])

        model_base = self.update_extr_transporters(model_base, extr_transporters=extr_transporters,
                                                   pools_list=self.tissue_pools[self.base_tissue] + ['env'],
                                                   tissue=self.base_tissue)

        reactions = [r.id for r in model_base.reactions]

        self.multitissue_reactions = reactions

        return model_base

    def add_other_tissues(self):
        """
        Second step in the reconstruction. Build the models for the other tissue with the modifications in the ids that
        are needed for the final integration
        Returns
        -------
        other_models: dict
            dictionary containing the model object for each tissue
        """

        other_tissues = {k: v for k, v in self.tissue_models.items() if k != self.base_tissue}

        other_models = {}

        for tissue in other_tissues:
            model_tissue = read_sbml_model(self.tissue_models[tissue])

            model_tissue = self.update_mets_ids(model=model_tissue, tissue=tissue,
                                                pools=self.tissue_pools[tissue])

            transporters = model_tissue.groups.get_by_id('transporters').members

            extr_transporters = [t for t in transporters if '__extr' in t.id or
                                 (t.id.startswith('T_') and t.id.endswith('__cyto'))]

            if tissue == 'leaf':
                reac_light = model_tissue.reactions.get_by_id('T_Light__chlo')
                extr_transporters.append(reac_light)

            drains = model_tissue.groups.get_by_id('drains').members

            extr_reacs = [r for r in model_tissue.reactions if '_extr' in r.id and r not in extr_transporters]

            all_excepts = extr_transporters + drains + extr_reacs

            internal = [reac for reac in model_tissue.reactions if reac not in all_excepts]

            model_tissue = self.update_internal_reactions(model=model_tissue, internal_reactions=internal,
                                                          tissue=tissue)

            model_tissue = self.update_drain_ids(model_tissue, drains=drains)

            model_tissue = self.update_extr_reactions(model_tissue, extr_reactions=extr_reacs,
                                                      pools_list=self.tissue_pools[tissue])

            model_tissue = self.update_extr_transporters(model_tissue, extr_transporters=extr_transporters,
                                                         pools_list=self.tissue_pools[tissue],
                                                         tissue=tissue)

            other_models[tissue] = model_tissue

            logging.info('model for ' + tissue + ' was created')

        return other_models

    def run_complete_pipeline(self):
        """
        Runs the complete pipeline to build the multitissue model.
        First, it creates the tissue base model, which is the tissue that receives most drains and contacts with the
        other tissues. Then, it creates the models for the other tissues.
        Second, it adds the reactions in the models of the other tissues and verifies if the base model has all the
        necessary transporters. If not, it creates them.
        Finally, it creates the totalbiomass reaction and updates the pathway groups in the final model.
        """

        model_base = self.add_base_tissue()

        logging.info('model for ' + self.base_tissue + ' was created')

        model_tissues = self.add_other_tissues()

        reactions_to_add = []

        groups_map = {}

        logging.info('integrating all tissues....')

        for tissue in model_tissues:
            model_tissue = model_tissues[tissue]
            pool = self.tissue_pools[tissue][0]
            for reaction in model_tissue.reactions:
                if reaction.id.endswith(pool + '_' + tissue) and reaction.id != 'T_Light__chlo_cp1_leaf':
                    reaction_base_id = reaction.id.replace(tissue, self.base_tissue)
                    if reaction_base_id not in self.multitissue_reactions:
                        reaction_base = reaction.copy()
                        reaction_base.id = reaction_base_id

                        groups = model_tissue.get_associated_groups(reaction)
                        for group in groups:
                            group.add_members([reaction_base])

                        for met in reaction_base.metabolites:
                            if met.id.endswith(tissue):
                                new_metid = met.id.replace(tissue, self.base_tissue)
                                try:
                                    new_met_obj = model_tissue.metabolites.get_by_id(new_metid)
                                except KeyError:
                                    new_met_obj = met.copy()
                                    new_met_obj.id = new_metid

                                reaction_base.add_metabolites({met: 0, new_met_obj: reaction_base.metabolites[met]},
                                                              combine=False)
                        reactions_to_add.append(reaction_base)
                        self.multitissue_reactions.append(reaction_base.id)

                        reaction_env_id = reaction_base_id.replace(pool, 'env')
                        if reaction_env_id not in self.multitissue_reactions:
                            reaction_env = reaction_base.copy()
                            reaction_env.id = reaction_env_id

                            groups = model_tissue.get_associated_groups(reaction)
                            for group in groups:
                                group.add_members([reaction_env])

                            for met in reaction_env.metabolites:
                                if met.id.endswith(pool):
                                    new_metid = met.id.replace(pool, 'env')
                                    try:
                                        new_met_obj = model_tissue.metabolites.get_by_id(new_metid)
                                    except KeyError:
                                        new_met_obj = met.copy()
                                        new_met_obj.id = new_metid

                                    reaction_env.add_metabolites({met: 0, new_met_obj: reaction_env.metabolites[met]},
                                                                 combine=False)
                            reactions_to_add.append(reaction_env)
                            self.multitissue_reactions.append(reaction_env.id)

                    reactions_to_add.append(reaction)
                    self.multitissue_reactions.append(reaction.id)

                elif reaction.id not in self.multitissue_reactions:
                    reactions_to_add.append(reaction)
                    self.multitissue_reactions.append(reaction.id)

            for group in model_tissue.groups:
                members = [r.id for r in group.members]
                if group.id not in groups_map:
                    groups_map[group.id] = members
                else:
                    for reac in members:
                        if reac not in groups_map[group.id]:
                            groups_map[group.id].append(reac)

        model_base.add_reactions(reactions_to_add)

        logging.info('create total biomass reaction....')

        model_base = self.create_total_biomass_reaction(model=model_base, metabolites=self.totalbio_comp)

        model_base.objective = 'TotalBiomass'

        for group in groups_map:
            group_obj = model_base.groups.get_by_id(group)
            for r in groups_map[group]:
                try:
                    r_obj = model_base.reactions.get_by_id(r)
                    group_obj.add_members([r_obj])
                except KeyError:
                    pass

        old_mets_to_remove = []
        for met in model_base.metabolites:
            if not met.reactions:
                old_mets_to_remove.append(met)

        model_base.remove_metabolites(old_mets_to_remove)

        write_sbml_model(model_base, os.path.join(PROJECT_PATH, 'reconstruction_results', self.model_id,
                                                  self.model_id + '_' + 'multissue_mature.xml'))

        logging.info('pipeline finished')


if __name__ == '__main__':
    modelid = 'vvinif2023'
    datasetid = 'RNAseq'

    tmodels = {'stem': 'stem_drains.xml',
               'leaf': 'leaf_drains.xml',
               'berry': 'berry_mature_drains.xml'}

    ttissues = {'leaf': ['T_CARBON-DIOXIDE__cyto', 'T_OXYGEN-MOLECULE__cyto', 'T_WATER__cyto', 'T_Light__chlo'],
                'berry': ['T_CARBON-DIOXIDE__cyto', 'T_OXYGEN-MOLECULE__cyto', 'T_WATER__cyto']}

    firsttissue = 'stem'

    tpools = {'stem': ['cp1', 'cp2'], 'leaf': ['cp1'], 'berry': ['cp2']}

    totalbiomass = {'e-Biomass__cyto_berry': -0.12, 'e-Biomass__cyto_leaf': -0.39, 'e-Biomass__cyto_stem': -0.49,
                    'TotalBiomass': 1}

    mt = Multitissue(model_id=modelid, dataset_id=datasetid, tissue_models=tmodels, base_tissue=firsttissue,
                     tissue_pools=tpools, tissue_drains=ttissues, totalbiomass_composition=totalbiomass)
    mt.add_base_tissue()
    mt.add_other_tissues()
    mt.run_complete_pipeline()
