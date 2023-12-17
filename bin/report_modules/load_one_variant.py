#!/usr/bin/env python3
#
# Project        : PSB Pipeline
# Filename       : load_one_variant.py
# Authors        : Chris Moth and R. Michael Sivley
# project Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : chris.moth@vanderbilt.edu
# Date           : 2023 Januar 27
# Description    : class CalculationResultsLoader is the "data" collection
#                : servant for the rep_one_variant.py code that creates
#                : the complete pipeline webpage for a single variant
# =============================================================================#

import os
import sys
import string
import json
from typing import Dict
from typing import List
import logging
import pandas as pd
import numpy as np

from lib import PDBMapGlobals
from psb_shared.ddg_monomer import DDG_monomer
from psb_shared.ddg_cartesian import DDG_cartesian
from psb_shared.ddg_repo import DDG_repo

LOGGER = logging.getLogger(__name__)
from lib import PDBMapTranscriptBase

def html_div_id_create(row_or_dict) -> str:
    """"
    It is convenient to be able to assemble an html div id from either a reporting dict or row
    """
    map_punct_to_None = dict.fromkeys(
        string.punctuation)  # for each punctuation character as key, the value is None

    return str("%s%s%s" % (row_or_dict.structure_id.strip(),
                           row_or_dict.chain_id.strip(),
                           '' if str(row_or_dict.mers.strip()) == 'monomer' else str(
                               row_or_dict.mers.strip()))).translate(
        str.maketrans(map_punct_to_None))

# This complex routine gathers PathProx, and ddG results from the pipeline run, where available
# It creates a single .html file representing the mutation, and returns a dictionary that is
# input to the overall case report in gathered_info

class CalculationResultsLoader:
    """
    Load all Pathprox, ddG, and any other calculation results for a complex into RAM
    Once loaded, other processes will inject these data into web pages.
    """

    def __init__(self, case_root_dir: str,
                 variant_directory_segment: str,
                 parent_report_row: Dict,
                 config_pathprox_dict: Dict):
        """
        @param variant_directory_segment, usually str of form "GENE_refseq_A123B"
        @param workstatus_filename:       Filename of status of all calculations, previously updated by psb_monitor.py,
        @param parent_report_row: In the event that no calculations were launched for the variant, then
           the unp from the parent row (on the main page report) populates variables so a minimal report can be
           completed
        @param config_dict: The dictionary obtained from the command line by the caller

        Load all the calculation results from the filesystem, making them available as dataframes and dictionaries
        indexed by (method_structureid_chain_mers)
        """
        if not os.path.exists(
                CalculationResultsLoader.case_root_dir_variant_directory_segment(case_root_dir,variant_directory_segment)):
            LOGGER.critical(
                "The  variant directory %s has not been created by the pipeline under %s.",
                variant_directory_segment, self._case_root_dir)
            sys.exit(1)

        self._config_pathprox_dict = config_pathprox_dict

        self._case_root_dir = case_root_dir
        self._variant_directory_segment = variant_directory_segment

        # This "full path" can be a spectacularly long entity - and
        # We probably need to force execution of psb_rep.py with UDN/UDN123456 as the current working directory
        # to stop other things from leaking out as issues.
        self._variant_directory_fullpath = os.path.join(self._case_root_dir, self._variant_directory_segment)

        # These dataframes are raw loads of .csv files created by the psb_plan.py and psb_monitor.py
        # over-arching pipeline control programs.  These are accessed via @properties of the same name
        # which lack the initial underscore.
        self._workplan_df = None  # Original list of jobs to be run, created by psb_plan.py
        self._workstatus_df = None  # pdb_monintor.py updates workplan rows with exit codes to create this
        self._structure_report_df = None  # Original list of structures on which calculations are planned run
        self._dropped_structures_df = None  # Original list of structures rejected by psb_plan.py

        # We build up dictionaries of dataframes, which link each 3D structure to the raw dataframes of computed
        # results for that structure
        self.ddG_monomer_results_dict_of_dfs = {}
        self.ddG_cartesian_results_dict_of_dfs = {}
        self.pathprox_disease1_results_dict_of_dfs = {}
        self.pathprox_disease2_results_dict_of_dfs = {}

        # self.rate4site_results_dict = {}

        # Once the above raw dataframes are loaded, we iterate through the outputs and gather
        # filenames for display images, PDBs to display in the NGL viewer, as well as variant lists.
        self.structure_graphics_dicts = []

        # it is possible that we only have a missense.csv file in the case directory so far,
        # In that case, we need to leave the constructor with all
        if not os.path.exists(self._variant_directory_fullpath):
            LOGGER.warning(
                "The  variant directory %s has not been created by psb_plan.py under %s.",
                variant_directory, self._case_root_dir)
            return

        # If psb_plan.py has been run, then we should have a workplan file with all the jobs laid planned out
        # However, we won't yet have a workstatus file until psb_launch.py program
        # In any case, make sure the data collection simply returns all empty dictionaries if no planning has been done
        self._workplan_filename_exists = False
        self._workstatus_filename_exists = False
        if os.path.exists(self.workplan_filename):
            self._workplan_filename_exists = True
            # We don't check for a workstatus file until we know there is a plan file.
            if os.path.exists(self.workstatus_filename):
                self._workstatus_filename_exists = True

        # Now that the *filename_exists booleans have been set above, we can access the properties
        # and get the dataframes back
        if len(self.workplan_df) > 0:
            if len(self.workstatus_df) == 0:
                LOGGER.warning("No workstatus file (%s) found.  However, %d jobs are to be launched.",
                               self.workstatus_filename,
                               len(self.workplan_df))
        else:
            self._info_dict = parent_report_row
            LOGGER.info("No work was planned for this variant")

        # Throughout all rows of the workstatus file, these columns are unchanged.  So just grab them as a reference
        if self._workstatus_filename_exists:
            self._info_dict = self.workplan_df.iloc[0]

            # Assert that we have a 'unp' in the workstatus dataframe, else
            # really the pipelien must halt altogether... at this juncture
            if 'unp' not in self._info_dict:
                msg = "No 'unp' column found in %s.  The report generator cannot continue" % (
                    self.workplan_filename)
                LOGGER.critical(msg)
                sys.exit(msg)
        else:  # Since we don't have any rows in the workstatus file, just use the information from the main page report
            self._info_dict = parent_report_row

#        self._unp_transcript = PDBMapTranscriptUniprot(self.unp)

    @property
    def project(self):
        if 'project' in self._info_dict:
            return self._info_dict['project']
        return args.projectORstructures

    # What fun - define read-only @property class members with one-line anonymous functions
    unp = property(lambda self: self._info_dict['unp'])
    gene = property(lambda self: self._info_dict['gene'])
    mutation = property(lambda self: self._info_dict['mutation'])
    # Try to move to "Variant" to describe that which was previously a mutation
    variant = property(lambda self: self._info_dict['mutation'])
    # unp_transcript = property(lambda self: self._unp_transcript)

    @property
    def notes(self):
        """
        @return: A string that might contain an overall condition, like no successful calculations
        """
        return ""

    def _makefilename(self, filename_segment: str) -> str:
        """
        @param filename_segment: "workstatus", "structure_report", or "dropped_structures"
        @return: os.path.join(self._case_root_dir,filename_segment)
        """

        # Example 'workstatus' becomes $UDN/variant_directory/gene_refseq_variant_workstatus.csv
        return os.path.join(self._variant_directory_fullpath,  # #UDN/variant_directory
                            "%s_%s.csv" % (self._variant_directory_segment, filename_segment))

    @property  # The name of the file listing all creations, created by psb_plan.py
    def workplan_filename(self) -> str:
        return self._makefilename("workplan")

    @property  # The filename created by psbb_launch.py, updated by psb_monitor.py:= workplan rows + job statuses
    def workstatus_filename(self) -> str:
        return self._makefilename("workstatus")

    @property  # The original listing of covering pdbs and models created by psb_plan.py
    def structure_report_filename(self) -> str:
        return self._makefilename("structure_report")

    @property  # The original listing of sub-optimal pdbs and models created by psb_plan.py
    def dropped_structures_filename(self) -> str:
        return self._makefilename("dropped_structures")

    @property
    def workplan_df(self) -> pd.DataFrame:
        if self._workplan_df is not None:
            return self._workplan_df

        if self._workplan_filename_exists:
            self._workplan_df = pd.read_csv(self.workplan_filename, delimiter='\t')
        else:
            msg = "Work plan file: %s not found.  Should have been created by psb_plan.py" % self.workplan_filename
            LOGGER.warning(msg)
            self._workplan_df = pd.DataFrame(columns=['uniquekey'])

        return self._workplan_df

    @property
    def workstatus_df(self) -> pd.DataFrame:
        # If nothing was done, return a known empty dataframe
        if self._workstatus_df is not None:
            return self._workstatus_df

        if self._workstatus_filename_exists:
            self._workstatus_df = pd.read_csv(self.workstatus_filename, delimiter='\t', keep_default_na=False,
                                              na_filter=False, dtype=str)
        else:
            self._workstatus_df = pd.DataFrame(columns=['uniquekey'])

        return self._workstatus_df

    @property
    def dropped_structures_df(self) -> pd.DataFrame:
        if self._dropped_structures_df is not None:
            return self._dropped_structures_df

        if not os.path.exists(self.dropped_structures_filename):
            self._dropped_structures_df = pd.DataFrame()
        else:
            # msg = "Work status file: %s not found.  Should have been created by pipeline" % \
            # self.dropped_structures_filename
            # LOGGER.critical(msg)
            # sys.exit(msg)
            self._dropped_structures_df = pd.read_csv(self.dropped_structures_filename, delimiter='\t')

        return self._dropped_structures_df

    @property
    def structure_report_df(self) -> pd.DataFrame:
        if self._structure_report_df is not None:
            return self._structure_report_df

        if not os.path.exists(self.structure_report_filename):
            LOGGER.warning("Structure Report file: %s not found.  Should have been created by pipeline psb_plan.py",
                           self.structure_report_filename)
            # Create an empty dataframe
            self._structure_report_df = pd.DataFrame(columns=['method', 'structure_id', 'chain_id', 'mers'])
        else:
            self._structure_report_df = pd.read_csv(self.structure_report_filename, delimiter='\t')

        return self._structure_report_df

    def load_Pathprox_results_dataframe(self, workplan_df_row: pd.Series, disease_1_or_2: str) -> (pd.Series, str):
        disease_variant_set_sql_label = self._config_pathprox_dict['%s_variant_sql_label' % disease_1_or_2]
        # import pdb; pdb.set_trace()
        # variant_short_description = self._config_pathprox_dict['%s_variant_short_description'%disease_1_or_2]

        # Define directory prefixes
        # This is where we'll search for those final pathprox .pdb files, etc
        output_flavor_directory = os.path.join(os.path.basename(workplan_df_row['outdir']), workplan_df_row['flavor'])

        pathprox_prefix = "%s_%s_%s_%s_D" % (
            workplan_df_row['pdbid'],
            workplan_df_row['chain'],
            disease_variant_set_sql_label,
            self._config_pathprox_dict['neutral_variant_sql_label'])

        # Load the Pathprox summary of results for this mutation
        pathprox_summary_filename = os.path.join(self._variant_directory_fullpath, output_flavor_directory,
                                                 "%s_summary.csv" % pathprox_prefix)
        try:
            pathprox_summary_df = pd.read_csv(pathprox_summary_filename, sep='\t')
        except FileNotFoundError:
            return None, "PathProx summary file not found: %s" % pathprox_summary_filename

        # Stuff the filesystem locations in the dataframe, so we don't have to recompute them later
        # when we grab graphics files, and variant json files, from these areas
        pathprox_summary_df['output_flavor_directory'] = output_flavor_directory
        pathprox_summary_df['pathprox_prefix'] = pathprox_prefix
        pathprox_summary_df['pdbmut'] = workplan_df_row.pdbmut
        pathprox_summary_df['mutation'] = workplan_df_row.mutation

        # Map dictionary keys to Pathprox-output .png filename ends that will be looked for in the filesystem
        # We only include neutral variant pictures for disease 1, because they are identical for disease 2
        # and it saves some space.
        pathprox_neutral_png_dict = {
            # Pathprox creates some legacy static png files using chimera_headless which we want on our
            # final report.  These *_png keynames will be rolled up at the end into the website .tar file.
            "variant_unknown_significance_png": "structure.png",
            "neutral_variants_png": "neutral.png",
            # Pathprox creates a variety of plots that describe important details about the calculation

            "neutral_K_plot_png": "neutral_K_plot.png",  # Most discriminating K plot, was _K before
        }

        pathprox_disease_png_dict = {
            "pathogenic_variants_png": "pathogenic.png",
            # Ripley's D results for the pathogenic
            "D_plot_png": "D_plot.png",

            # Load the PathProx Mapping and Performance for the pathogenic calculations
            "pp_plot_png": "pathprox.png",
            "roc_plot_png": "pathprox_roc.png",
            "pr_plot_png": "pathprox_pr.png",

            "pathogenic_K_plot_png": "pathogenic_K_plot.png"  # Most discriminating K plot, was _K before

        }

        pathprox_png_dict = pathprox_disease_png_dict.copy()
        if disease_1_or_2 == 'disease1':
            pathprox_png_dict.update(pathprox_neutral_png_dict)

        LOGGER.debug("Retrieving pathprox image pngs from %s" % os.path.join(self._variant_directory_segment,
                                                                             output_flavor_directory))

        pathprox_png_filename_dict = dict.fromkeys(pathprox_png_dict.keys(), '')

        for pathprox_png_key, pathprox_png_filename_ending in pathprox_png_dict.items():
            # We can save some space by NOT loading the neutral png files for the second disease.
            # These are the same as the first one, and it saves space to skip them
            png_filename = os.path.join(self._variant_directory_segment, output_flavor_directory,
                                        "%s_%s" % (pathprox_prefix, pathprox_png_filename_ending))
            if os.path.exists(png_filename):
                LOGGER.debug("Pathprox putput png file %s found" % png_filename)
                # The actual filename in the final web html will be relative to the isoform_variant directory
                pathprox_png_filename_dict[pathprox_png_key] = os.path.join(output_flavor_directory,
                                                                            os.path.basename(png_filename))
            else:
                # This is not a problem - just means pathprox did not do this particular calculation
                LOGGER.debug("Pathprox putput png file %s not found" % png_filename)

        pathprox_png_filename_pd = pd.DataFrame(pathprox_png_filename_dict, index=[0])

        pathprox_summary_df = pathprox_summary_df.join(pathprox_png_filename_pd)

        # An odd thing about the pathprox output .csv file is that some columns are named with the mutation
        # So, Let's create alias columns with standard names
        # I.e. pathprox_summary_df['pathprox'] = pathprox_summary_df['M123V_pathprox']
        #      pathprox_summary_df['neutcon'] = pathprox_summary_df['M123V_neutcont']
        #      pathprox_summary_df['pathcon'] = pathprox_summary_df['M123V_pathcon']
        mutation = workplan_df_row.mutation
        for score_string in ['pathprox', 'neutcon', 'pathcon']:
            pathprox_summary_df[score_string] = np.nan
            if "%s_%s" % (mutation, score_string) in pathprox_summary_df.columns:
                pathprox_summary_df[score_string] = pathprox_summary_df["%s_%s" % (mutation, score_string)]
        return pathprox_summary_df, None

    def load_ddG_results_dataframe(self, workplan_df_row: pd.Series, calculation_flavor: str) -> (pd.DataFrame, str):
        ddg_repo = DDG_repo(PDBMapGlobals.config['ddg_config'],
                            calculation_flavor=calculation_flavor)

        structure_id = workplan_df_row['pdbid']
        chain_id = workplan_df_row['chain']
        if workplan_df_row['method'] == 'swiss':
            ddg_structure_dir = ddg_repo.set_swiss(structure_id, chain_id)
        elif workplan_df_row['method'] == 'modbase':
            ddg_structure_dir = ddg_repo.set_modbase(structure_id, chain_id)
        elif workplan_df_row['method'] == 'alphafold':
            ddg_structure_dir = ddg_repo.set_alphafold(structure_id, chain_id)
        elif workplan_df_row['method'] == 'usermodel':
            ddg_structure_dir = ddg_repo.set_usermodel(structure_id, chain_id)
        else:  # Has to be a pdb
            ddg_structure_dir = ddg_repo.set_pdb(structure_id.lower(), chain_id)

        variant = workplan_df_row['pdbmut']
        ddg_repo.set_variant(variant)
        ddg_monomer_or_cartesian = DDG_monomer(ddg_repo,
                                               variant) if calculation_flavor == 'ddG_monomer' else DDG_cartesian(
            ddg_repo, variant)

        ddg_results_df = ddg_monomer_or_cartesian.retrieve_result()
        if ddg_results_df is None:
            msg = 'No %s results found for %s.%s:%s' % (calculation_flavor, structure_id, chain_id, variant)
            LOGGER.warning(msg)
            return None, msg

        return ddg_results_df, None

    def load_dataframes(self):
        """
        Iterate through the workplan file, xfer to workstatus, and load numerical outputs for all
        calculations, for all structures representing the isoform.

        After this call, the following Dataframes will be filled, or left at None if no calculation was performed
            self.ddg_results_dict_of_dfs
            self.pathprox_disease1_results_dict_of_dfs
            self.pathprox_disease2_results_dict_of_dfs

        @return:
        """

        structure_report_df_indexed = self.structure_report_df.set_index(['method', 'structure_id', 'chain_id', 'mers'],
                                                                         drop=False)
        workstatus_df_indexed = self.workstatus_df.set_index(['uniquekey'], drop=False)
        for workplan_uniquekey, workplan_row in self.workplan_df.set_index(['uniquekey'], drop=False).iterrows():
            try:
                workstatus_row = workstatus_df_indexed.loc[workplan_uniquekey]
            except KeyError:
                LOGGER.critical("psb_launch.py did not create status record for Calculation %s" % workplan_uniquekey)
                workstatus_row = workplan_row.copy()
                workstatus_row['Notes'] = "Not Launched"
                workstatus_row['ExitCode'] = ''

            if workstatus_row['ExitCode'] == '0':
                # Excellent go harvest raw outputs to later combine into reports
                method_pdbid_chain_mers_tuple = (
                    workstatus_row['method'], workstatus_row['pdbid'], workstatus_row['chain'], workstatus_row['mers'])

                # Today, Sequence Annotation results are not included in web outputs
                if workstatus_row['flavor'] == 'SequenceAnnotation':
                    continue

                try:
                    structure_row = structure_report_df_indexed.loc[method_pdbid_chain_mers_tuple]
                except KeyError:
                    msg = "Structure %s referenced from job is not in the structure file" % str(
                        method_pdbid_chain_mers_tuple)
                    LOGGER.critical(msg)
                    sys.exit(msg)

                if workstatus_row['flavor'].endswith(self._config_pathprox_dict['disease1_variant_sql_label']):
                    pathprox_disease1_results_df, msg = self.load_Pathprox_results_dataframe(
                        workstatus_row, 'disease1')
                    if pathprox_disease1_results_df is None:
                        workstatus_row['Notes'] = msg
                    else:
                        self.pathprox_disease1_results_dict_of_dfs[
                            method_pdbid_chain_mers_tuple] = pathprox_disease1_results_df

                elif workstatus_row['flavor'].endswith(self._config_pathprox_dict['disease2_variant_sql_label']):
                    pathprox_disease2_results_df, msg = self.load_Pathprox_results_dataframe(
                        workstatus_row, 'disease2')
                    if pathprox_disease2_results_df is None:
                        workstatus_row['Notes'] = msg
                    else:
                        # structure_row['disease2_pathprox'] = pathprox_disease2_results_df['pathprox']
                        # structure_row['disease2_pr_aoc'] = pathprox_disease2_results_df['pr_auc']
                        self.pathprox_disease2_results_dict_of_dfs[
                            method_pdbid_chain_mers_tuple] = pathprox_disease2_results_df

                elif workstatus_row['flavor'] in ['ddG_monomer', 'ddG_cartesian']:
                    ddg_results_df, msg = self.load_ddG_results_dataframe(workstatus_row, workstatus_row['flavor'])
                    if ddg_results_df is None:
                        workstatus_row['Notes'] = msg
                    else:
                        if workstatus_row['flavor'] == 'ddG_monomer':
                            # structure_row['ddG_monomer'] = ddg_results_df['ddG']
                            self.ddG_monomer_results_dict_of_dfs[method_pdbid_chain_mers_tuple] = ddg_results_df
                        else:
                            self.ddG_cartesian_results_dict_of_dfs[method_pdbid_chain_mers_tuple] = ddg_results_df
                else:
                    die_msg = "flavor %s in workstatus row is neither %s nor %s nor ddG_monomer not ddG_cartesian - cannot continue:\n%s" % (
                        workstatus_row['flavor'],
                        self._config_pathprox_dict['disease1_variant_sql_label'],
                        self._config_pathprox_dict['disease2_variant_sql_label'],
                        str(workstatus_row))
                    LOGGER.critical(die_msg)
                    sys.exit(die_msg)
            else:
                workstatus_row['Notes'] = workstatus_row['JobState']

    def _load_pathprox_residues_of_interest_json(self, pathprox_result: pd.Series) -> Dict:
        """
        Supports load_structure_graphics_dicts by loading json information that was output in a specific pathprox run
        Replaces the pdbSSfilename entry with a relative path suitable for web processing.
        @return: pathprox_output_json dictionary
        """
        pathprox_output_json_filename = os.path.join(self._variant_directory_segment,
                                                     pathprox_result['output_flavor_directory'],
                                                     pathprox_result['pathprox_prefix'] + "_ResiduesOfInterest.json")
        if not os.path.exists(pathprox_output_json_filename):
            LOGGER.warning(
                "Pathprox left no %s file\nwhich would have contained variants and a PDB filename",
                pathprox_output_json_filename)
            return {}
        else:
            with open(pathprox_output_json_filename) as json_f:
                pathprox_output_json = json.load(json_f)

            # if 'cifSSfilename' in pathprox_output_json:
            #    cifSSbasename = os.path.basename(pathprox_output_json['cifSSfilename'])
            #    cifSSdirname = os.path.realpath(os.path.dirname(pathprox_output_json['cifSSfilename']))
            #    pathprox_output_json['cifSSfilename'] = \
            #    os.path.join(os.path.relpath(cifSSdirname,variant_report_directory),cifSSbasename)

            def _variant_count(variant_type: str):
                return sum([len(chain_variants[1]) for chain_variants in pathprox_output_json[variant_type].items()])

            LOGGER.info('Loaded %d variants, %d neutrals, %d pathogenic from %s' % (
                _variant_count('variants'),
                _variant_count('neutrals'),
                _variant_count('pathogenics'),
                pathprox_output_json_filename))

            return pathprox_output_json

    def _load_pathprox_scores_all_residues_json(self, pathprox_result: pd.Series) -> Dict:
        """
        Supports load_structure_graphics_dicts by loading json information that was output in a specific pathprox run
        Replaces the pdbSSfilename entry with a relative path suitable for web processing.
        @return: pathprox_output_json dictionary
        """
        pathprox_output_json_filename = os.path.join(self._variant_directory_segment,
                                                     pathprox_result['output_flavor_directory'],
                                                     pathprox_result[
                                                         'pathprox_prefix'] + "_renum_pathprox.json")

        if not os.path.exists(pathprox_output_json_filename):
            LOGGER.warning(
                "Pathprox left no all-residues pathprox.json output file: %s",
                pathprox_output_json_filename)
            return {}
        else:
            with open(pathprox_output_json_filename) as json_f:
                pathprox_output_json = json.load(json_f)
        return pathprox_output_json

    def _load_alpha_fold_metrics_all_residues_json(self, pathprox_result: pd.Series) -> Dict:
        """
        Supports load_structure_graphics_dicts by loading json information that was output in a specific pathprox run

        @return: alpha_fold_metrics list
        """
        alpha_fold_metrics_json = os.path.join(self._variant_directory_segment,
                                               pathprox_result['output_flavor_directory'],
                                               pathprox_result['pathprox_prefix'] + "_alphafold_metrics.json")

        if not os.path.exists(alpha_fold_metrics_json):
            LOGGER.warning(
                "Pathprox left no alphafold model metrics json output file: %s",
                alpha_fold_metrics_json)
            return {}
        else:
            with open(alpha_fold_metrics_json) as json_f:
                alpha_fold_metrics = json.load(json_f)
        return alpha_fold_metrics

    def _load_cosmis_json(self, pathprox_result: pd.Series) -> Dict:
        """
        Supports load_structure_graphics_dicts by loading cosmis scores for each chain left from the pathprox runs

        @return: Rate4Site json file.  Dictionary per chain, with then values per residue.
        """
        cosmis_scores_json_filename = os.path.join(self._variant_directory_segment,
                                                   pathprox_result['output_flavor_directory'],
                                                   pathprox_result['pathprox_prefix'] + "_cosmis.json")

        if not os.path.exists(cosmis_scores_json_filename):
            LOGGER.warning(
                "Pathprox left no Cosmis json output file: %s",
                cosmis_scores_json_filename)
            return {}
        else:
            with open(cosmis_scores_json_filename) as json_f:
                cosmis_scores = json.load(json_f)
        return cosmis_scores

    def _load_rate4site_json(self, pathprox_result: pd.Series) -> Dict:
        """
        Supports load_structure_graphics_dicts by loading rate4site scores for each chain left from the pathprox runs

        @return: Rate4Site json file.  Dictionary per chain, with then values per residue.
        """
        rate4site_scores_json_filename = os.path.join(self._variant_directory_segment,
                                                      pathprox_result['output_flavor_directory'],
                                                      pathprox_result['pathprox_prefix'] + "_rate4site.json")

        if not os.path.exists(rate4site_scores_json_filename):
            LOGGER.warning(
                "Pathprox left no Rate4Site json output file: %s",
                rate4site_scores_json_filename)
            return {}
        else:
            with open(rate4site_scores_json_filename) as json_f:
                rate4site_scores = json.load(json_f)
        return rate4site_scores








    def load_structure_graphics_dicts(self) -> List[str]:
        """
            self.structure_graphics_dicts will be populated in same structure sequence as self.structure_report_df

        @return: additional_website_filelist to roll into the final website
        """
        additional_website_filelist = []

        # For each structure, gather the supplementary program outputs, png files, PDB file
        # references created by the pathprox jobs.

        structure_report_df_appended = self.structure_report_df.copy(). \
            set_index(['method', 'structure_id', 'chain_id', 'mers'], drop=False)

        structure_report_df_appended['ddG_monomer'] = np.nan
        for ddg_results_key in self.ddG_monomer_results_dict_of_dfs:
            structure_report_df_appended.at[ddg_results_key, 'ddG_monomer'] = \
                self.ddG_monomer_results_dict_of_dfs[ddg_results_key].ddG

        structure_report_df_appended['ddG_cartesian'] = np.nan
        for ddg_results_key in self.ddG_cartesian_results_dict_of_dfs:
            structure_report_df_appended.at[ddg_results_key, 'ddG_cartesian'] = \
                self.ddG_cartesian_results_dict_of_dfs[ddg_results_key].total

        structure_report_df_appended['disease1_pathprox'] = np.nan
        structure_report_df_appended['disease1_pr_auc'] = np.nan
        for pathprox_results_key in self.pathprox_disease1_results_dict_of_dfs:
            structure_report_df_appended.at[pathprox_results_key, 'disease1_pathprox'] = \
                self.pathprox_disease1_results_dict_of_dfs[pathprox_results_key].pathprox
            structure_report_df_appended.at[pathprox_results_key, 'disease1_pr_auc'] = \
                self.pathprox_disease1_results_dict_of_dfs[pathprox_results_key].pr_auc

        structure_report_df_appended['disease2_pathprox'] = np.nan
        structure_report_df_appended['disease2_pr_auc'] = np.nan
        for pathprox_results_key in self.pathprox_disease2_results_dict_of_dfs:
            structure_report_df_appended.at[pathprox_results_key, 'disease2_pathprox'] = \
                self.pathprox_disease2_results_dict_of_dfs[pathprox_results_key].pathprox
            structure_report_df_appended.at[pathprox_results_key, 'disease2_pr_auc'] = \
                self.pathprox_disease2_results_dict_of_dfs[pathprox_results_key].pr_auc

        structure_report_df_appended['html_div_id'] = ''
        # If the dataframe is empty, (i.e. no structural coverage) then do not run through the apply process,
        # because it returns a bad thing....



        if not structure_report_df_appended.empty:
            structure_report_df_appended['html_div_id'] = structure_report_df_appended.apply(
                lambda df_row: html_div_id_create(df_row), axis=1)

        # Let's round some columns for sane output
        structure_report_df_appended = structure_report_df_appended.round(
            {
                'ddG_monomer': 2,
                'ddG_cartesian': 2,
                'disease1_pathprox': 2,
                'disease1_pr_auc': 3,
                'disease2_pathprox': 2,
                'disease2_pr_auc': 3
            }
        )

        # Combine all the dataframes into a structure_graphics_dict with "root" elements
        # that dereference specific dictionary sets.
        for index, structure in structure_report_df_appended.iterrows():
            structure_graphics_dict = structure.fillna('').to_dict()
            structure_graphics_dict['html_div_id'] = html_div_id_create(structure)
            structure_key = (structure.method, structure.structure_id, structure.chain_id, structure.mers)
            alpha_fold_metrics = None
            rate4site_scores = None
            cosmis_scores = None

            for disease1_or_2, pathprox_results_dict_of_dfs in zip(
                    ['disease1', 'disease2'],
                    [self.pathprox_disease1_results_dict_of_dfs, self.pathprox_disease2_results_dict_of_dfs]):
                structure_graphics_dict["pathprox_%s_results" % disease1_or_2] = {}
                structure_graphics_dict['%s_pathogenics' % disease1_or_2] = {}
                if structure_key in pathprox_results_dict_of_dfs:
                    pathprox_result_series_for_structure = pathprox_results_dict_of_dfs[structure_key].iloc[0]
                    pathprox_output_json = self._load_pathprox_residues_of_interest_json(
                        pathprox_result_series_for_structure)

                    if not alpha_fold_metrics and structure.method == "alphafold":
                        alpha_fold_metrics = self._load_alpha_fold_metrics_all_residues_json(
                            pathprox_result_series_for_structure)
                        if alpha_fold_metrics:
                            structure_graphics_dict['ngl_alpha_fold_metrics'] = \
                                '[' + ', '.join(
                                    str(alpha_fold_metric) for alpha_fold_metric in alpha_fold_metrics[1]) + ']'

                    if not rate4site_scores:
                        rate4site_scores = self._load_rate4site_json(pathprox_result_series_for_structure)
                        if rate4site_scores:
                            ngl_formatted_residue_rate4site_pairs = []
                            for chain in rate4site_scores:
                                for residue_no in rate4site_scores[chain]:
                                    ngl_formatted_residue_rate4site_pairs.append(
                                        ("%s:%s" %
                                         (residue_no, chain),
                                         rate4site_scores[chain][residue_no]['score'])
                                    )

                            # Now create the final javascript-compatible version of this...
                            javascript_dict_residues_rate4site_scores = \
                                "{" + ", ".join(["'%s': %s" % dict_key \
                                                 for dict_key in ngl_formatted_residue_rate4site_pairs]) + \
                                "}"

                            # Add the Rate4site score json format to the dictionary seen in the django template.
                            structure_graphics_dict['ngl_rate4site_scores'] = \
                                javascript_dict_residues_rate4site_scores

                    if not cosmis_scores:
                        cosmis_scores = self._load_cosmis_json(pathprox_result_series_for_structure)
                        if cosmis_scores:
                            ngl_formatted_residue_cosmis_pairs = []
                            for chain in cosmis_scores:
                                for residue_no in cosmis_scores[chain]:
                                    ngl_formatted_residue_cosmis_pairs.append(
                                        ("%s:%s" %
                                         (residue_no, chain),
                                         cosmis_scores[chain][residue_no]['cosmis'])
                                    )

                            # Now create the final javascript-compatible version of this...
                            javascript_dict_residues_cosmis_scores = \
                                "{" + ", ".join(["'%s': %s" % dict_key \
                                                 for dict_key in ngl_formatted_residue_cosmis_pairs]) + \
                                "}"

                            # Add the Cosmis score json format to the dictionary seen in the django template.
                            structure_graphics_dict['ngl_cosmis_scores'] = \
                                javascript_dict_residues_cosmis_scores

                    if 'pdbSSfilename' in pathprox_output_json:
                        pdbSSbasename = os.path.basename(pathprox_output_json['pdbSSfilename'])

                        structure_graphics_dict['pdbSSfilename'] = os.path.join(
                            # self._variant_directory_segment,
                            pathprox_result_series_for_structure['output_flavor_directory'],
                            pdbSSbasename)

                        additional_website_filelist.append(
                            os.path.join(self._variant_directory_fullpath, structure_graphics_dict['pdbSSfilename']))

                        # This only works if we are running the report in the same directory system in
                        # which we ran the report - which is lousy for development
                        # pdbSSdirname = os.path.realpath(os.path.dirname(pathprox_output_json['pdbSSfilename']))
                        # structure_graphics_dict['pdbSSfilename'] = os.path.join(
                        #     os.path.relpath(pdbSSdirname, self._variant_directory_fullpath), pdbSSbasename)

                    # Integrate pathprox-output png files to the website tar build.
                    for pathprox_results_key in pathprox_result_series_for_structure.index:
                        if pathprox_results_key.endswith('_png'):
                            # IF we have a filename in this, and the filename does not contain neutral and
                            if pathprox_result_series_for_structure[pathprox_results_key]:
                                additional_website_filelist.append(os.path.join(
                                    self._variant_directory_fullpath,
                                    pathprox_result_series_for_structure[pathprox_results_key]))

                    structure_graphics_dict[
                        "pathprox_%s_results" % disease1_or_2] = pathprox_result_series_for_structure.fillna(
                        "").to_dict()
                    # For final .html - we need to be specific regarding type of pathogenic.  We weren't when json file was written
                    structure_graphics_dict['%s_pathogenics' % disease1_or_2] = pathprox_output_json['pathogenics']

                    if pathprox_output_json:
                        pathprox_scores_all_residues = self._load_pathprox_scores_all_residues_json(
                            pathprox_result_series_for_structure)
                        if pathprox_scores_all_residues:
                            ngl_formatted_residue_pathprox_pairs = []
                            for chain in pathprox_scores_all_residues:
                                for residue_no in pathprox_scores_all_residues[chain]:
                                    ngl_formatted_residue_pathprox_pairs.append(
                                        ("%s:%s" % (residue_no, chain),
                                         pathprox_scores_all_residues[chain][residue_no])
                                    )

                            javascript_dict_residues_pathprox_scores = \
                                "{" + ", ".join(["'%s': %s" % dict_key \
                                                 for dict_key in ngl_formatted_residue_pathprox_pairs]) + \
                                "}"
                            structure_graphics_dict['ngl_%s_pathprox_scores' % disease1_or_2] = \
                                javascript_dict_residues_pathprox_scores

                    def pathprox_residue_list_to_ngl_format(residues_onechain, chain_id):
                        """
                        Space out residues from pathprox into a single string for ngl viewer display
                        Append a chain_id to each residue, if supplied
                        @param residues_onechain:  pathprox output residues
                        @param chain_id:  pathprox output chain ID
                        @return: string ready for ngl display
                        """
                        if len(residues_onechain) == 0:
                            return 0, None
                        ngl_residue_string = ''
                        for transcript_resno in residues_onechain:
                            if ngl_residue_string:
                                ngl_residue_string += ' '
                            ngl_residue_string += "%d" % transcript_resno
                            if chain_id:
                                ngl_residue_string += ":%s" % chain_id

                        return len(residues_onechain), ngl_residue_string
                        # def residue_list_to_ngl_CAs(residues_onechain):
                        #    ngl_res_count,ngl_res_string = residue_list_to_ngl(residues_onechain)
                        #    return ngl_res_count,"(" + ngl_res_string + ") and .CA" if ngl_res_string else ngl_res_string

                    def residues_to_ngl(residues):
                        """
                        Given raw pathprox chain/residue oiutput, convert each chain to ngl format.
                        @param residues:
                        @return:
                        """
                        ngl_res_count = 0
                        ngl_res_string = None
                        for chain_id in residues:
                            ngl_res_chain_count, ngl_res_chain_string = pathprox_residue_list_to_ngl_format(
                                residues[chain_id], chain_id)
                            if ngl_res_chain_count and ngl_res_chain_string:
                                if ngl_res_string:
                                    ngl_res_string += ' '
                                else:
                                    ngl_res_string = ''

                                ngl_res_count += ngl_res_chain_count
                                ngl_res_string += ngl_res_chain_string
                        return ngl_res_count, ngl_res_string

                    def residues_to_ngl_CAs(residues):
                        ngl_res_count, ngl_res_string = residues_to_ngl(residues)
                        return ngl_res_count, "(" + ngl_res_string + ") and .CA" if ngl_res_string else ngl_res_string

                    structure_graphics_dict['unique_div_name'] = pathprox_output_json['unique_div_name']
                    structure_graphics_dict['ngl_variant_residue_count'], structure_graphics_dict[
                        'ngl_variant_residues'] = residues_to_ngl(
                        pathprox_output_json['variants'])
                    structure_graphics_dict['ngl_neutral_residue_count'], structure_graphics_dict[
                        'ngl_neutral_residues'] = residues_to_ngl_CAs(
                        pathprox_output_json['neutrals'])
                    structure_graphics_dict['ngl_%s_residue_count' % disease1_or_2], \
                    structure_graphics_dict['ngl_%s_residues' % disease1_or_2] = residues_to_ngl_CAs(
                        pathprox_output_json['pathogenics'])
                    # structure_graphics_dict['ngl_%s_pathprox_scores' % disease1_or_2] = pathprox_scores_all_residues

            self.structure_graphics_dicts.append(structure_graphics_dict)
        return additional_website_filelist

    @staticmethod
    def case_root_dir_variant_directory_segment(case_root_dir: str, variant_directory_segment: str):
        return os.path.join(case_root_dir, variant_directory_segment)

    def _load_ddG_details(self):
        self.workstatus_df = pd.read_csv(self.workstatus_filename, '\t', keep_default_na=False, na_filter=False,
                                         dtype=str)
        msg = "%d rows read from work status file %s" % (len(self.workstatus_df), self.workstatus_filename)
        if len(self.workstatus_df) < 1:
            LOGGER.critical(msg)
            return None
        else:
            if not infoLogging:
                print(msg)
            LOGGER.info(msg)

    """ OBSOLETEdef ddG_html_vars(structure, thestruct):
        if not 'ddG_monomer' in thestruct:
            return structure, None

        df_row = thestruct['ddG_monomer']

        ddg_repo = DDG_repo(config_dict['ddg_config'],
                            calculation_flavor='ddg_monomer')

        output_flavor_directory = os.path.join(df_row['outdir'], 'ddG', 'ddg')

        structure_id = df_row['Structure ID']
        chain_id = df_row['Chain']
        if df_row['method'] == 'swiss':
            ddg_structure_dir = ddg_repo.set_swiss(structure_id, chain_id)
        elif df_row['method'] == 'modbase':
            ddg_structure_dir = ddg_repo.set_modbase(structure_id, chain_id)
        elif df_row['method'] == 'usermodel':
            ddg_structure_dir = ddg_repo.set_usermodel(structure_id, chain_id)
        else:  # Has to be a pdb
            ddg_structure_dir = ddg_repo.set_pdb(structure_id.lower(), chain_id)

        variant = df_row['pdbmut']
        ddg_repo.set_variant(variant)
        ddg_monomer = DDG_monomer(ddg_repo, variant)

        ddg_results_df = ddg_monomer.retrieve_result()
        if ddg_results_df is None:
            LOGGER.info('No results found for %s.%s:%s' % (structure_id, chain_id, variant))
            return structure, None

        return structure, ddg_results_df
    """

    def load_calculations(self) -> None:
        """
        We iterate over all the jobs in the workstatus file.  For each one, we load
        a small dataframe from the file system, which is suitable for direct conversion
        to an .html table.

        We build up a variant wide dataframe of all the calculation results,

        Simultaneously, we create an dictory so that "All calculation results for structure abcd.pdf"
        are readily loadable
        """

        for index, row in self.workstatus_df.iterrows():
            # We now set about loading pathprox, ddG, etc calculation
            # results from the file system.
            # Each calculation is loaded into its own dataframe

            # The final .html (.pdf) output is organized structure by structure.
            # To stay sane, we somewhat need a dictionary with one element per structure
            # Here, we'll have struct_dict[Method/pdbid/chain] with sub hash elements:
            # mutation/pdbid/chain/
            # COSMIC dict/Clinvar dict/ddg Dict are dictionaries which have the entire rows
            # From the completed job status

            struct_dict = {}
            for index, row in df_complete_jobs.iterrows():
                method_pdbid_chain_mers_tuple = (row['Method'], row['pdbid'], row['chain'], row['mers'])
                thestruct = struct_dict.get(method_pdbid_chain_mers_tuple, {})

                thestruct['mutation'] = row['mutation']
                thestruct['pdbid'] = row['pdbid']
                thestruct['chain'] = row['chain']
                thestruct['mers'] = row['mers']
                if row['flavor'].endswith(self._config_pathprox_dict['disease1_variant_sql_label']):
                    thestruct['disease1'] = row.to_dict()
                elif row['flavor'].endswith(self._config_pathprox_dict['disease2_variant_sql_label']):
                    thestruct['disease2'] = row.to_dict()
                elif 'SequenceAnnotation' == row['flavor']:
                    continue  # UDN Sequence annotations are NOT a part of generated reports
                elif 'ddG_monomer' == row['flavor']:
                    thestruct['ddG_monomer'] = row.to_dict()
                elif 'ddG_cartesian' == row['flavor']:
                    thestruct['ddG_cartesian'] = row.to_dict()
                else:
                    LOGGER.critical("flavor in row is neither %s nor %s nor ddG_monomer - cannot continue:\n%s",
                                    self._config_pathprox_dict['disease1_variant_sql_label'],
                                    self._config_pathprox_dict['disease2_variant_sql_label'], str(row))
                    sys.exit(1)
                # This will often reassign over prior assignments - That's OK
                struct_dict[method_pdbid_chain_mers_tuple] = thestruct
