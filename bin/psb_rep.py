#!/usr/bin/env python3
#
# Project        : PSB Pipeline
# Filename       : psb_monitor.py
# Authors        : Chris Moth and R. Michael Sivley
# project Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : chris.moth@vanderbilt.edu
# Date           : 2018-01-22
# Description    : Output final html and pdf reports from completed jobs
#                : (output from psb_launch.py)
# =============================================================================#


"""\
Create final html and pdf reports for completed pipeline jobs from either
a project id (ex. UDN123456) or a single workstatus.csv file output from psb_monitor.py

Configuration options must be provided with
   -c global.config file
   -u user.config overrides
   structure_report.csv originally output by psb_plan.py
   workstatus.csv as previously updated with psb_monitor.py

HTML is created by populating template .html through jinja2 template library
"""

import datetime
import grp
import logging
import os
# import argparse, configparser
import pprint
import string
import sys
import time
from array import array
from logging.handlers import RotatingFileHandler

from psb_shared.ddg_monomer import DDG_monomer
from psb_shared.ddg_cartesian import DDG_cartesian
from psb_shared.ddg_repo import DDG_repo
import shutil
import numpy as np
import pandas as pd
from jinja2 import Environment, FileSystemLoader
from io import StringIO, BytesIO
import pycurl
import json
from xml2json import Cxml2json

import subprocess
from lib import PDBMapSQLdb
from lib import PDBMapTranscriptBase
from lib import PDBMapTranscriptUniprot
from lib import PDBMapTranscriptEnsembl
from lib import PDBMapComplex

from collections import defaultdict

from typing import Dict
# from subprocess import Popen, PIPE
# from weasyprint import HTML

from psb_shared import psb_config

print("4 of 4  %s Pipeline report generator.  -h for detailed help" % __file__)

sh = logging.StreamHandler()
LOGGER = logging.getLogger()
LOGGER.addHandler(sh)

# Now that we've added streamHandler, basicConfig will not add another handler (important!)
log_format_string = '%(asctime)s %(levelname)-4s [%(filename)16s:%(lineno)d] %(message)s'
date_format_string = '%H:%M:%S'
log_formatter = logging.Formatter('%(asctime)s %(levelname)-4s [%(filename)20s:%(lineno)d] %(message)s',
                                  datefmt="%H:%M:%S")

LOGGER.setLevel(logging.DEBUG)
sh.setLevel(logging.INFO)
sh.setFormatter(log_formatter)

root_dir_log_filename = "psb_rep.log"
need_roll = os.path.isfile(root_dir_log_filename)
rootdir_fh = RotatingFileHandler(root_dir_log_filename, backupCount=7)

rootdir_fh.setFormatter(log_formatter)
rootdir_fh.setLevel(logging.INFO)
LOGGER.addHandler(rootdir_fh)

if need_roll:
    rootdir_fh.doRollover()

sys.stderr.write("Log file for entire case is %s\n" % root_dir_log_filename)

# As report is generated, build up a list of just the files needed to make a portable website
# This list omits the intermediate files left behind in the calculation
website_filelist = ['html/']

cmdline_parser = psb_config.create_default_argument_parser(__doc__, os.path.dirname(os.path.dirname(__file__)))

cmdline_parser.add_argument("-s", "--slurm",
                            help="Create a slurm file to launch all the reports", action="store_true")
cmdline_parser.add_argument("projectORstructures", type=str,
                            help="The project ID UDN123456 to report on all mutations.  Else, a specific structures file from psb_plan.py  Example: ....../UDN/UDN123456/GeneName_NM_12345.1_S123A_structure_report.csv",
                            default=os.path.basename(os.getcwd()), nargs='?')
cmdline_parser.add_argument("workstatus", nargs='?', type=str,
                            help="Blank for all mutations, or specific output file from psb_monitor.py\n" +
                                 "Example: ....../UDN/UDN123456/GeneName_NM_12345.1_S123A_workstatus.csv")
args, remaining_argv = cmdline_parser.parse_known_args()

# A period in the argument means the user wants to monitor one mutation only,
# directly from a single mutation output file of psb_launch.py
oneMutationOnly = (
        '.' in args.projectORstructures and os.path.isfile(args.projectORstructures) and args.workstatus is not None)
infoLogging = False

if args.debug:
    infoLogging = True
    sh.setLevel(logging.DEBUG)
    rootdir_fh.setLevel(logging.DEBUG)
elif args.verbose:
    infoLogging = True
    sh.setLevel(logging.INFO)

# Load the structure details of all jobs that were not dropped

# print "Command Line Arguments"
LOGGER.info("Command Line Arguments:\n%s" % pprint.pformat(vars(args)))

if not os.path.exists(args.config):
    LOGGER.critical("Global config file not found: " + args.config)
    sys, exit(1)

required_config_items = ['output_rootdir', 'collaboration', 'ddg_config', 'rate4site_dir']

config, config_dict = psb_config.read_config_files(args, required_config_items)
config_dict_shroud_password = {x: config_dict[x] for x in required_config_items}
dbpass = config_dict.get('dbpass', '?')
config_dict_shroud_password['dbpass'] = '*' * len(dbpass)
LOGGER.info("Configuration File parameters:\n%s" % pprint.pformat(config_dict_shroud_password))

try:
    slurmParametersAll = dict(config.items("SlurmParametersAll"))
except Exception as ex:
    slurmParametersAll = {}

import copy

slurmParametersReport = copy.deepcopy(slurmParametersAll)

try:
    slurmParametersReport.update(dict(config.items("SlurmParametersReport")))
except Exception as ex:
    pass

config_pathprox_dict = dict(config.items("PathProx"))
LOGGER.info("Pathprox config entries:\n%s" % pprint.pformat(config_pathprox_dict))

# The collaboration_dir is the master directory for the case, i.e. for one patient
# Example: /dors/capra_lab/projects/psb_collab/UDN/UDN532183
udn_root_directory = os.path.join(config_dict['output_rootdir'], config_dict['collaboration'])

if oneMutationOnly:
    collaboration_absolute_dir = os.path.dirname(os.path.dirname(args.projectORstructures))
else:
    collaboration_absolute_dir = os.path.join(udn_root_directory, args.projectORstructures)

collaboration_absolute_dir = os.path.realpath(collaboration_absolute_dir)
LOGGER.info("Collaboration_absolute_dir = %s", collaboration_absolute_dir)

collaboration_dir = os.path.relpath(os.getcwd(), collaboration_absolute_dir)
LOGGER.info("Relative to cwd=%s, collaboration_directory= %s", os.getcwd(), collaboration_dir)


def copy_html_css_javascript():
    # The source html directory, and it's large set of supporting files, is located beneath the directory
    # where this file, psb_rep.py, is located.  All mutations share this one large html/css/javascript
    # sourcee directory
    src = os.path.join(os.path.dirname(os.path.realpath(__file__)), "html")
    dest = "%s/html" % collaboration_dir
    LOGGER.info("Copying supporting javascript and html from %s to %s" % (src, dest))

    try:
        shutil.rmtree(dest, ignore_errors=True)
        # Save some space by not coping the data from nglviewer
        shutil.copytree(src, dest, ignore=shutil.ignore_patterns('*/nglview/data/*'))
    except Exception as e:
        LOGGER.exception(e)
        sys.exit(1)

    # Make sure
    src_file_list = [os.path.join(root, name) for root, dirs, files in os.walk(src) for name in files]

    dest_file_list = [os.path.join(root, name) for root, dirs, files in os.walk(dest) for name in files]
    if len(src_file_list) == len(dest_file_list):
        infostr = "Copy of %d html support files succeeded" % len(src_file_list)
        print(infostr)
        LOGGER.info(infostr)
        # print '\n'.join(src_file_list)
    else:
        LOGGER.critical("FAILURE to cp html support files (%d source files,%d dest files)" % (
            len(src_file_list), len(dest_file_list)))
        sys.exit(1)

    # Now give read access to everyone, read-write to group, execute on directories
    try:
        capra_group = grp.getgrnam('capra_lab').gr_gid
    except KeyError:
        capra_group = os.getegid()

    os.chmod(dest, 0o775)
    os.chown(dest, -1, capra_group)

    for root, dirs, files in os.walk(dest):
        for dir in dirs:
            os.chmod(os.path.join(root, dir), 0o775)
            os.chown(os.path.join(root, dir), -1, capra_group)
        for file in files:
            fname = os.path.join(root, file)
            os.chmod(fname, 0o664)
            os.chown(os.path.join(fname), -1, capra_group)


# Return a dictionary of "generic Interactions" that can flow into the summary report
# Return an html table that can be placed at end of summary report
# Return html long report
# Return text that can be dumped to the log file
def gene_interaction_report(case_root, case, CheckInheritance):
    genepair_dict_file = os.path.join(case_root, 'casewide', 'DiGePred', '%s_all_gene_pairs_summary.json' % case)

    try:
        with open(genepair_dict_file, 'r') as fd:
            genepair_dict = json.load(fd)
    except FileNotFoundError:
        LOGGER.exception("Cannot open %s.  Did you run Souhrid's analysis scripts?" % genepair_dict_file)
        return None, None, None, None

    # genes_filename = os.path.join(case_root,"%s_genes.txt"%case)
    # try:
    #    with open(genes_filename) as fd:
    #       unstripped_genes = fd.readlines()
    # except:
    #    LOGGER.exception("Cannot open %s.  Did you run Souhrid's analysis scripts?"%genes_filename)
    #     return None,None,None,None

    geneInteractions = {}
    html_table = StringIO()
    html_report = StringIO()
    text_report = StringIO()

    structure_positive_genes = [gene.strip() for gene in list(sorted(genepair_dict.keys()))]

    pairs = {}
    pairs['direct'] = []
    pairs['pathways'] = []
    pairs['phenotypes'] = []
    pairs['other'] = []

    # We let the template html open the <table> with formatting as it prefers
    # print >>html_table, "<table>" 
    print("<thead><tr>", file=html_table)
    print("<th></th>", file=html_table)  # Top left is blank cell in header
    for gene1 in structure_positive_genes:
        print("<th>" + gene1 + "</th>", file=html_table)
    print("</tr></thead>", file=html_table)

    for gene1 in structure_positive_genes:
        print("%s:" % gene1, file=text_report)
        print('%s:' % gene1, file=html_report)
        print("<tr><td>%s</td>" % gene1, file=html_table)
        if gene1 in genepair_dict:
            print('<ul style="list-style: none;">', file=html_report)
        for gene2 in structure_positive_genes:
            print("<td>", file=html_table)
            if gene1 in genepair_dict and gene2 in genepair_dict[gene1] and (
                    (not CheckInheritance) or ('nherit' in genepair_dict[gene1][gene2])):
                pairs_hits = []
                if 'direct' in genepair_dict[gene1][gene2]:
                    pairs_hits.append('direct')
                if 'pwy' in genepair_dict[gene1][gene2]:
                    pairs_hits.append('pathways')
                if 'pheno' in genepair_dict[gene1][
                    gene2]:  # and genepair_dict[gene1][gene2]['Patient_Phenotype_Overlap'] > 0:
                    pairs_hits.append('phenotypes')
                # if genepair_dict[gene1][gene2]['Distance'] < 3 and genepair_dict[gene1][gene2]['Coexpression'] < 50:
                #    pairs_hits.append('other')
                if pairs_hits:
                    for pairs_hit in pairs_hits:
                        pairs[pairs_hit].append((gene1, gene2))
                        if pairs_hit != 'phenotypes':  # We don't get too excited for phenotypes at this point
                            if not gene1 in geneInteractions:
                                geneInteractions[gene1] = []
                            if not gene2 in geneInteractions[gene1]:
                                geneInteractions[gene1].append(gene2)
                    logging.debug("Gene Interactions: %s %s %s" % (gene1, gene2, str(pairs_hits)))
                    print(gene2, str(pairs_hits), file=text_report)
                    print("<li>", gene2 + ':', ", ".join(pairs_hits), "</li>", file=html_report)
                    print('<br>'.join(pairs_hits), file=html_table)
            print("</td>", file=html_table)
        print('</ul>', file=html_report)
        print("</tr>", file=html_table)
    # print("</table>", file=html_table)
    return geneInteractions, html_table.getvalue(), html_report.getvalue(), text_report.getvalue()


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
# It creates a sigle .html file representing the mutation, and returns a dictionary that is
# input to the overall case report in gathered_info

class CalculationResultsLoader:
    """
    Load all Pathprox, ddG, and any other calculation results for a complex into RAM
    Once loaded, other processes will inject these data into web pages.
    """

    def __init__(self, variant_directory_segment: str, parent_report_row: Dict = None):
        """
        @param variant_directory_segment, usually str of form "GENE_refseq_A123B"
        @param workstatus_filename:       Filename of status of all calculations, previously updated by psb_monitor.py,
        @param parent_report_row: In the event that no calculations were launched for the variant, then
           the unp from the parent row (on the main page report) populates variables so a minimal report can be
           completed

        Load all the calculation results from the filesystem, making them available as dataframes and dictionaries
        indexed by (method_structureid_chain_mers)
        """
        if not os.path.exists(
                CalculationResultsLoader.collaboration_dir_variant_directory_segment(variant_directory_segment)):
            LOGGER.critical(
                "The  variant directory %s has not been created by the pipeline under %s.",
                variant_directory_segment, collaboration_dir)
            sys.exit(1)

        self._variant_directory_segment = variant_directory

        # This "full path" can be a spectacularly long entity - and
        # We probably need to force execution of psb_rep.py with UDN/UDN123456 as the current working directory
        # to stop other things from leaking out as issues.
        self._variant_directory_fullpath = os.path.join(collaboration_dir, self._variant_directory_segment)

        if not os.path.exists(os.path.join(collaboration_dir, variant_directory_segment)):
            LOGGER.critical(
                "The  variant directory %s has not been created by the pipeline under %s.",
                variant_directory, collaboration_dir)
            sys.exit(1)

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

        self._workstatus_filename_exists = False
        if os.path.exists(self.workstatus_filename):
            self._workstatus_filename_exists = True

        if not self._workstatus_filename_exists:
            if len(self.workplan_df) > 0:
                LOGGER.critical("No workstatus file (%s) found.  However, %d jobs should have been launched.", \
                                len(self.workplan_df))
            else:
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

        self._unp_transcript = PDBMapTranscriptUniprot(self.unp)

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
    unp_transcript = property(lambda self: self._unp_transcript)

    @property
    def notes(self):
        """
        @return: A string that might contain an overall condition, like no successful calculations
        """
        return ""

    def _makefilename(self, filename_segment: str) -> str:
        """
        @param filename_segment: "workstatus", "structure_report", or "dropped_structures"
        @return: os.path.join(collaboration_dir,filename_segment)
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

        if not os.path.exists(self.workplan_filename):
            msg = "Work plan file: %s not found.  Should have been created by psb_plan.py" % self.workplan_filename
            LOGGER.critical(msg)
            sys.exit(msg)

        self._workplan_df = pd.read_csv(self.workplan_filename, delimiter='\t')
        return self._workplan_df

    @property
    def workstatus_df(self) -> pd.DataFrame:
        # If nothing was done, return a known empty dataframe
        if not self._workstatus_filename_exists:
            return self._workplan_df

        if self._workstatus_df is not None:
            return self._workstatus_df

        if not os.path.exists(self.workstatus_filename):
            msg = "Work status file: %s not found.  Should have been created by psb_monitor.py" % \
                  self.workstatus_filename
            LOGGER.critical(msg)
            sys.exit(msg)

        self._workstatus_df = pd.read_csv(self.workstatus_filename, delimiter='\t', keep_default_na=False,
                                          na_filter=False, dtype=str)
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
            msg = "Structure Report file: %s not found.  Should have been created by pipeline psb_plan.py" % self.structure_report_filename
            LOGGER.critical(msg)
            sys.exit(msg)

        self._structure_report_df = pd.read_csv(self.structure_report_filename, delimiter='\t')
        return self._structure_report_df

    def load_Pathprox_results_dataframe(self, workplan_df_row: pd.Series, disease_1_or_2: str) -> (pd.Series, str):
        disease_variant_set_sql_label = config_pathprox_dict['%s_variant_sql_label' % disease_1_or_2]
        # variant_short_description = config_pathprox_dict['%s_variant_short_description'%disease_1_or_2]

        # Define directory prefixes
        # This is where we'll search for those final pathprox .pdb files, etc
        output_flavor_directory = os.path.join(os.path.basename(workplan_df_row['outdir']), workplan_df_row['flavor'])

        pathprox_prefix = "%s_%s_%s_%s_D" % (
            workplan_df_row['pdbid'],
            workplan_df_row['chain'],
            disease_variant_set_sql_label,
            config_pathprox_dict['neutral_variant_sql_label'])

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
        ddg_repo = DDG_repo(config_dict['ddg_config'],
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

                if workstatus_row['flavor'].endswith(config_pathprox_dict['disease1_variant_sql_label']):
                    pathprox_disease1_results_df, msg = self.load_Pathprox_results_dataframe(
                        workstatus_row, 'disease1')
                    if pathprox_disease1_results_df is None:
                        workstatus_row['Notes'] = msg
                    else:
                        self.pathprox_disease1_results_dict_of_dfs[
                            method_pdbid_chain_mers_tuple] = pathprox_disease1_results_df

                elif workstatus_row['flavor'].endswith(config_pathprox_dict['disease2_variant_sql_label']):
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
                        config_pathprox_dict['disease1_variant_sql_label'],
                        config_pathprox_dict['disease2_variant_sql_label'],
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

    def load_structure_graphics_dicts(self):
        """
            self.structure_graphics_dicts will be populated in same structure sequence as self.structure_report_df

        @return:
        """

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
            structure_report_df_appended['html_div_id'] = structure_report_df_appended.apply(lambda df_row: \
                                                                                                 html_div_id_create(
                                                                                                     df_row), axis=1)

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

            for disease1_or_2, pathprox_results_dict_of_dfs in zip(
                    ['disease1', 'disease2'],
                    [self.pathprox_disease1_results_dict_of_dfs, self.pathprox_disease2_results_dict_of_dfs]):
                structure_graphics_dict["pathprox_%s_results" % disease1_or_2] = {}
                structure_graphics_dict['%s_pathogenics' % disease1_or_2] = {}
                if structure_key in pathprox_results_dict_of_dfs:
                    pathprox_result_series_for_structure = pathprox_results_dict_of_dfs[structure_key].iloc[0]
                    pathprox_output_json = self._load_pathprox_residues_of_interest_json(
                        pathprox_result_series_for_structure)
                    vus_chain = next(iter(pathprox_output_json['variants']))
                    vus_residue = pathprox_output_json['variants'][vus_chain][0]

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

                        # self.rate4site_results_dict = rate4site_scores[vus_chain][str(vus_residue)]

                    if 'pdbSSfilename' in pathprox_output_json:
                        pdbSSbasename = os.path.basename(pathprox_output_json['pdbSSfilename'])

                        structure_graphics_dict['pdbSSfilename'] = os.path.join(
                            # self._variant_directory_segment,
                            pathprox_result_series_for_structure['output_flavor_directory'],
                            pdbSSbasename)

                        website_filelist.append(
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
                                website_filelist.append(os.path.join(
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

    @staticmethod
    def collaboration_dir_variant_directory_segment(variant_directory_segment: str):
        return os.path.join(collaboration_dir, variant_directory_segment)

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
                if row['flavor'].endswith(config_pathprox_dict['disease1_variant_sql_label']):
                    thestruct['disease1'] = row.to_dict()
                elif row['flavor'].endswith(config_pathprox_dict['disease2_variant_sql_label']):
                    thestruct['disease2'] = row.to_dict()
                elif 'SequenceAnnotation' == row['flavor']:
                    continue  # UDN Sequence annotations are NOT a part of generated reports
                elif 'ddG_monomer' == row['flavor']:
                    thestruct['ddG_monomer'] = row.to_dict()
                elif 'ddG_cartesian' == row['flavor']:
                    thestruct['ddG_cartesian'] = row.to_dict()
                else:
                    LOGGER.critical("flavor in row is neither %s nor %s nor ddG_monomer - cannot continue:\n%s",
                                    config_pathprox_dict['disease1_variant_sql_label'],
                                    config_pathprox_dict['disease2_variant_sql_label'], str(row))
                    sys.exit(1)
                # This will often reassign over prior assignments - That's OK
                struct_dict[method_pdbid_chain_mers_tuple] = thestruct


def _initialize_local_logging_in_variant_directory(variant_directory: str):
    # pw_name = pwd.getpwuid( os.getuid() ).pw_name # example jsheehaj or mothcw
    log_filename = os.path.join(CalculationResultsLoader.collaboration_dir_variant_directory_segment(variant_directory),
                                "psb_rep.log")

    if oneMutationOnly:
        sys.stderr.write("psb_rep log file is %s\n" % log_filename)
    else:
        LOGGER.info("additionally logging to file %s" % log_filename)
    _need_roll = os.path.isfile(log_filename)
    local_fh = RotatingFileHandler(log_filename, backupCount=7)
    formatter = logging.Formatter('%(asctime)s %(levelname)-4s [%(filename)20s:%(lineno)d] %(message)s',
                                  datefmt="%H:%M:%S")
    local_fh.setFormatter(formatter)
    if args.debug:
        infoLogging = True
        local_fh.setLevel(logging.DEBUG)
    elif args.verbose:
        local_fh.setLevel(logging.INFO)

    LOGGER.addHandler(local_fh)
    if _need_roll:
        local_fh.doRollover()

    return local_fh  # Be sure to _end_local_logging(local_fh) on this later


def _end_local_logging(local_fh: RotatingFileHandler):
    local_fh.flush()
    local_fh.close()
    LOGGER.removeHandler(local_fh)


class PfamDomainGraphics:
    """
    Manages loading of graphics in form of received json string from pfam.xfam.org.
    Important: In cases where pfam.xfam.org lacks JSON graphics, or fails to return
    a graphic in timely fashion, we build the JSON ourselves.

    I've not documented the format because you can see it online:
    http://pfam.xfam.org/protein/P45381/graphic

    with examples of how to display the graphic and the table format, not far from that.  Ex:
    http://pfam.xfam.org/protein/P45381/
    """

    def __init__(self, uniprot_id, variant):
        """
        Download (or build) a PFAM-graphics-library-compatiable dictionary of predicted protein domain
        annotations, so that pipeline case users can see a linear graphic layout of each protein

        @param uniprot_id: The canonical, or non-canonical uniprot identifier for the protein to diagram
        @param variant:    Our variant of interest, to add to the graphics and the table beneath.
        """
        self.unp = uniprot_id
        self.variant = variant

        assert self.unp is not None

    def _fetch_canonical_pfam_graphic_json_from_xfam(self, timeout_seconds=60) -> str:
        """
        Fetch a canonical PFAM graphic JSON string from the internet, using curl, for the
        unp passed to the class constructor.

        @return: PFAM graphic JSON string or None if pfam.xfam server down or bad data returned

        Requires config_dict['interpro_dir'] to be defined
        Inside that directory must be 'match_humanonly.xml' file
        """

        # xfam lacks graphics URLs for non-canonical isoforms.  We assemble those ourselves, elsewhere.
        canonical_pfam_xfam_url = 'http://pfam.xfam.org/protein/{0}/graphic'.format(
            self.unp.split(' ')[0].split('-')[0])

        LOGGER.debug("Attempting curl connection to %s to fetch canonical pfam graphics.  Timeout=%d\n",
                     canonical_pfam_xfam_url, timeout_seconds)
        buffer = BytesIO()

        try:
            c = pycurl.Curl()
            c.setopt(pycurl.CONNECTTIMEOUT, timeout_seconds)
            c.setopt(pycurl.TIMEOUT, timeout_seconds)
            # c.setopt(pycurl.WRITEFUNCTION, lambda x: None)
            c.setopt(pycurl.WRITEFUNCTION, buffer.write)
            c.setopt(pycurl.HEADERFUNCTION, lambda x: None)
            c.setopt(pycurl.NOSIGNAL, 1)
            c.setopt(pycurl.URL, canonical_pfam_xfam_url)
            c.setopt(pycurl.HTTPGET, 1)

            c.perform()
            c.close()

        except Exception as ex:
            LOGGER.exception('Failed to get JSON string from %s: %s' % (canonical_pfam_xfam_url, str(ex)))
            return None

        domain_graphics_json = buffer.getvalue().decode('latin')

        if len(domain_graphics_json) > 3 and domain_graphics_json[0] == '[':
            domain_graphics_json = domain_graphics_json[1:-1]
        elif domain_graphics_json[0] == '<':
            # Sometimes we get back html in the response, to tell us that the PFAM website is down
            LOGGER.warning("XML, not json, was returned from pfam:\n%s", domain_graphics_json)
            return None
        LOGGER.debug("pfam sent us:\n%s", domain_graphics_json)
        return domain_graphics_json

    def create_graphics_legend_and_graphics_dict(self) -> (str, Dict):
        """
        From self.unp, create a graphics Legend, initialize self._domain_graphics_dict.
        This involves opening a large xml file, and either getting the JSON from the internet, or constructing it manually
        
        @return: domain_graphics_legend,domain_graphics_dict
        """
        domain_graphics_legend = "Pfam Domain Graphic unavailable"
        # Get Uniprot Graphics string or put in a blank
        domain_graphics_json = None
        if self.unp:  # Fetch the JSON - wait up to 30 seconds
            interpro_xml_handle = Cxml2json(os.path.join(config_dict['interpro_dir'], 'match_humanonly.xml'))
            if interpro_xml_handle.isCanonical(self.unp):  # Then we need to manually create the xml ourselves
                # Go for the web link because this is a canonical UNP - however that _can_ fail.
                domain_graphics_legend = "Downloaded Pfam Domain Graphic for Canonical Isoform %s" % self.unp
                LOGGER.info("Attempting download of Domain Graphics for %s from xfam.pfam" % self.unp)
                domain_graphics_json = self._fetch_canonical_pfam_graphic_json_from_xfam(30)
                if not domain_graphics_json:
                    LOGGER.info("Download returned nothing")
                    domain_graphics_legend = "Pfam Domain Graphic for Canonical Isoform %s" % self.unp
            else:
                domain_graphics_legend = "Pfam Domain Graphic for Non-canonical Isoform %s" % self.unp

            # _Either_ communications failure OR non-canonical isoform
            # So we create our own graphic from our local xml database
            if not domain_graphics_json:  # _Either_ communications failure OR non-canonical (no communication)
                LOGGER.info("Creating Domain Graphic for %s from xml", self.unp)
                nonCanonicalGraphicsJSON = interpro_xml_handle.PFAMgraphicsJSONfromUnp(self.unp)
                domain_graphics_json = json.dumps(nonCanonicalGraphicsJSON)

        domain_graphics_dict = json.loads(domain_graphics_json)

        return domain_graphics_legend, domain_graphics_dict

    def add_our_variant_to_graphics_dict(self, domain_graphics_dict: Dict) -> Dict:
        """
        Specifically add _our_ variant of interest to the pfam graphics, as a table entry, and as a diamond

        @return: reference to same domain_graphics_dict, with the ['markups'] altered.
        """
        # Add our variant point of interest to the PFAM-like domain graphics.
        # as a simple diamond in the graphics, and an entry in the table.
        variant_site_markup_dict = {'colour': '#e469fe',
                                    'display': True,
                                    'headStyle': 'diamond',
                                    'lineColour': '#333333',
                                    'metadata': {'database': 'UDN Case',
                                                 'description': '%s' % self.variant,
                                                 'start': self.variant[1:-1],
                                                 'type': 'Mutation'},
                                    'residue': 'X',
                                    'start': self.variant[1:-1],
                                    'type': 'UDN Mutation site',
                                    'v_align': 'top'}

        # Grab any _existing_ markups in our domain graphics.
        # Else default to an empty list.  Add _our_ variant markup to that list of markups.
        domain_graphics_dict['markups'] = domain_graphics_dict.get('markups', []) + [variant_site_markup_dict]

        # Now repackage in a list, so the "domain_graphics_json" is back to original form....
        return domain_graphics_dict

    def update_graphics_dict_hrefs_to_point_to_pfam(selfself, domain_graphics_dict: Dict) -> Dict:
        """
        Downloaded href's embedded in the graphics dictionary assume that the domains are described locally.
        We must convert them to point to pages on the UK webserver.

        @param domain_graphics_dict:
        @return:
        """

        for region in domain_graphics_dict.get('regions', []):
            if 'href' in region:
                region['href'] = 'http://pfam.xfam.org' + region['href']

        return domain_graphics_dict

    def create_pfam_html(self, doman_graphics_legend: str, domain_graphics_dict: Dict) -> str:
        """
        @param doman_graphics_legend: 
        @param domain_graphics_json: 
        @return html string:
        """

        # Build a table to explaine the pfam graphics
        PfamResultTableHTML = '''<table class="resultTable details"
                 id="imageKey"
                 summary="Key for the Pfam domain image">
            <thead>
              <tr>
                <th class="dh" rowspan="2">Source</th>
                <th class="dh" rowspan="2">Domain</th>
                <th class="dh" rowspan="2">Start</th>
                <th class="dh" rowspan="2">End</th>
                <th class="sh" style="display: none" colspan="2">Gathering threshold (bits)</th>
                <th class="sh" style="display: none" colspan="2">Score (bits)</th>
                <th class="sh" style="display: none" colspan="2">E-value</th>
              </tr>
              <tr>
                <th class="sh" style="display: none">Sequence</th>
                <th class="sh" style="display: none">Domain</th>
                <th class="sh" style="display: none">Sequence</th>
                <th class="sh" style="display: none">Domain</th>
                <th class="sh" style="display: none">Sequence</th>
                <th class="sh" style="display: none">Domain</th>
              </tr>
            </thead>'''

        PfamResultTableHTML += '<tbody>'

        table_rows = []
        for motif in domain_graphics_dict.get('motifs', []):
            table_rows.append({
                'type': motif['metadata']['type'],
                'domain': 'n/a',
                'start': motif['start'],
                'end': motif['end']})

        for region in domain_graphics_dict.get('regions', []):
            # It's nice to say something about a domain when we can
            # However not all subsequences of a protein is a domain of course
            domain_table_entry = ''
            if 'text' in region:
                region_text = region['text']
            else:
                region_text = 'Text Missing'

            if 'href' in region:
                domain_table_entry = "<a href=" + region['href'] + ">" + region_text + "</a>"
            else:
                domain_table_entry = region_text
            table_rows.append({
                'type': region['metadata']['type'],
                'domain': domain_table_entry,
                'start': region['start'],
                'end': region['end']})

        table_rows.append({
            'type': 'Variant',
            'domain': 'n/a',
            'start': self.variant[1:-1],
            'end': self.variant[1:-1]})

        sorted_table_rows = sorted(table_rows, key=lambda d: int(d['start']))

        for sr in sorted_table_rows:
            table_row = '<tr class="odd">'
            table_row += '<td class="domain">%s</td>' % sr['type']
            table_row += '<td><span class="inactive">%s</span></td>' % sr['domain']
            table_row += '<td>%s</td>' % sr['start']
            table_row += '<td>%s</td>' % sr['end']
            table_row += 6 * '<td class="sh" style="display: none"><span class="inactive">n/a</span></td>'
            table_row += '</tr>\n'
            PfamResultTableHTML += table_row

        PfamResultTableHTML += '</tbody></table>'

        return PfamResultTableHTML


def report_one_variant_one_isoform(variant_directory_segment: str, parent_report_row: Dict) -> Dict:
    """
    Generate extensive structure analysis report .html for a single variant of a single transcript by
    reading all PathPRox/ddG/etc calculations into dataframes, transforming the data a bit, and then
    writing out into .html via jinja2 templates.

    Variant_isoform summary data are returned to the caller in an extensive dictionary, which the caller
    uses to fill in tables on the main case-wide page.

    @param variant_directory: Typically a string like "Gene_transcriptId_aaVariant" that is the root subdirectory
            of the collaboration directory that houses all calculation results.  The string is also a prefix to
             key files that record calculation status

    @param parent_report_row: In the event that no calculations were launched for the variant, then
           the unp from the parent row (on the main page report) populates variables so a minimal report can be
           completed
    @return: An variant_isoform_summary dictionary that can help the caller assemble the main report page.
    """

    # Start logging - after which logging entries are mirrored to variant_directory/psb_rep.log
    variant_directory_fullpath = CalculationResultsLoader.collaboration_dir_variant_directory_segment(
        variant_directory_segment)
    local_logger_fh = _initialize_local_logging_in_variant_directory(variant_directory_fullpath)

    # The rate4site score for the variant is shown in a table of transcripts near top, if rate4site data is found
    rate4site_dict = {}
    if 'transcript' in parent_report_row and parent_report_row['transcript']:
        ensembl_transcript_ids =  parent_report_row['transcript'].split(';')
        for ensembl_transcript_id in ensembl_transcript_ids:
            rate4site_norm_rates_filename = os.path.join(config_dict['rate4site_dir'],
                                                     "%s_norm_rates.txt" % ensembl_transcript_id.split('.')[0])
            if os.path.exists(rate4site_norm_rates_filename):
                rate4site_df, alpha_parameter, average, std = PDBMapComplex._load_one_rate4site_file(
                    rate4site_norm_rates_filename)
                variant_pos = int(parent_report_row['mutation'][1:-1])
                # Now grab just th dataframe row of all the rate4site entries for this position
                rate4site_dict[ensembl_transcript_id] = rate4site_df[rate4site_df['pos'] == variant_pos].iloc[0].to_dict()
            else:
                LOGGER.warning("No rate4site results were found for %s: %s", parent_report_row['unp'], ensembl_transcript_id)

    calculation_results_loader = CalculationResultsLoader(variant_directory_segment, parent_report_row)
    calculation_results_loader.load_dataframes()
    calculation_results_loader.load_structure_graphics_dicts()

    # Various data max/min/etc computations from the gathered data, which are great to return to the caller
    # and also useful to compose the report.
    variant_isoform_summary = {}


    # Load the template psb_report.html that is hear
    env = Environment(loader=FileSystemLoader(os.path.join(os.path.dirname(os.path.realpath(__file__)))))
    variant_isoform_template = env.get_template("psb_report.html")
    structure_iframe_template = env.get_template("html/structureIframeTemplate.html")

    calculation_results_loader_as_dict = {key: value for key, value in vars(calculation_results_loader).items() \
                                          if not (
            key.startswith('_'))}  # -> don't exclude callable, you lose properties callable(value)))},

    for attribute in ['project', 'gene', 'unp', 'variant', 'unp_transcript']:
        calculation_results_loader_as_dict[attribute] = getattr(calculation_results_loader, attribute)

    # Generate the htrml for the doman graphics
    pfamGraphicsIframe_fname = "%s_%s_PfamGraphicsIframe.html" % (
        calculation_results_loader.unp,
        calculation_results_loader.variant)

    template_render_dict = {
        # Convert calculation_results to a dictionary because it is infinitely more compatable with
        # jinja2 template
        'calculation_results': calculation_results_loader_as_dict,
        # 'structure_report_df_appended': structure_report_df_appended,
        "today_date": time.strftime("%Y-%m-%d"),
        "pfamGraphicsIframe_fname": pfamGraphicsIframe_fname,
        # html template references disease variant strings like clinvar/cosmic...
        'config_pathprox_dict': config_pathprox_dict,
        'rate4site_dict': rate4site_dict

    }

    # Write out the individual structure.html files
    # website_filelist.append('color_chains_consecutively.js')
    for structure_graphics_dict in calculation_results_loader_as_dict['structure_graphics_dicts']:
        LOGGER.info("Creating NGL Viewer for %s", structure_graphics_dict['html_div_id'])
        structureIframe_out = structure_iframe_template.render(
            {'structure': structure_graphics_dict,
             'config_pathprox_dict': config_pathprox_dict}
        )
        structureIframe_base_filename = os.path.join(variant_directory_segment,
                                                     "%s_viewer.html" % structure_graphics_dict['html_div_id'])
        website_filelist.append(structureIframe_base_filename)
        with open(structureIframe_base_filename, "w") as html_f:
            html_f.write(structureIframe_out)

    if LOGGER.isEnabledFor(logging.DEBUG):
        pp = pprint.PrettyPrinter(indent=1)
        LOGGER.debug("Dictionary to render:\n%s" % pp.pformat(template_render_dict))

    variant_report_html_out = variant_isoform_template.render(template_render_dict)

    variant_report_base_filename = os.path.join(variant_directory_segment,
                                                "%s_report" % (variant_directory_segment))

    variant_isoform_summary['html_filename'] = variant_report_base_filename + ".html"
    website_filelist.append(variant_isoform_summary['html_filename'])
    with open(variant_isoform_summary['html_filename'], "w") as html_f:
        html_f.write(variant_report_html_out)

    pfg = PfamDomainGraphics(calculation_results_loader.unp, calculation_results_loader.variant)
    pfam_graphics_legend, pfam_domain_graphics_dict = pfg.create_graphics_legend_and_graphics_dict()
    pfam_domain_graphics_dict = pfg.add_our_variant_to_graphics_dict(pfam_domain_graphics_dict)
    pfam_domain_graphics_dict = pfg.update_graphics_dict_hrefs_to_point_to_pfam(pfam_domain_graphics_dict)
    template_render_dict["PfamResultTableHTML"] = pfg.create_pfam_html(pfam_graphics_legend, pfam_domain_graphics_dict)

    # Encapsulate in a list
    template_render_dict["DomainGraphicsJSON"] = '[' + json.dumps(pfam_domain_graphics_dict) + ']'

    template_pfam_graphics = env.get_template("html/pfamGraphicsIframeTemplate.html")
    pfam_graphics_html_full_iframe = template_pfam_graphics.render(template_render_dict)

    # The variant_report refers to the Pfam graphics without a directory indicator.
    # So we prefix the filename, here, with the report directory to write it.
    pfamGraphicsIframe_fname_with_directory = os.path.join(variant_directory_segment, pfamGraphicsIframe_fname)
    website_filelist.append(pfamGraphicsIframe_fname_with_directory)
    with open(pfamGraphicsIframe_fname_with_directory, "w") as html_f:
        html_f.write(pfam_graphics_html_full_iframe)

    # With the report fully written to disk, now populate summary details for return to caller.
    variant_isoform_summary['AA_seq'] = calculation_results_loader.unp_transcript.aa_seq
    variant_isoform_summary['AA_len'] = len(variant_isoform_summary['AA_seq'])

    # gathered_info['gene_id']  PDBMapProtein.unp2gene_id(gethered_info['unp'])

    # template_vars = gathered_info.copy()

    # mutation = row0['mutation']  # Sorry for alias - but this is all over the code

    variant_isoform_summary['Error'] = ''  # < This is False for truth but looks fine on the report
    variant_isoform_summary['ddG Monomer Max'] = variant_isoform_summary['ddG Monomer Min'] = None
    variant_isoform_summary['ddG Cartesian Max'] = variant_isoform_summary['ddG Cartesian Min'] = None
    variant_isoform_summary['disease1_pp Min'] = variant_isoform_summary['disease1_pp Max'] = None
    variant_isoform_summary['disease2_pp Min'] = variant_isoform_summary['disease2_pp Max'] = None

    if calculation_results_loader.ddG_monomer_results_dict_of_dfs:
        ddG_list = [calculation_results_loader.ddG_monomer_results_dict_of_dfs[ddg_results_key].ddG \
                    for ddg_results_key in calculation_results_loader.ddG_monomer_results_dict_of_dfs]
        variant_isoform_summary['ddG Monomer Max'] = max(ddG_list)
        variant_isoform_summary['ddG Monomer Min'] = min(ddG_list)

    if calculation_results_loader.ddG_cartesian_results_dict_of_dfs:
        ddG_list = [calculation_results_loader.ddG_cartesian_results_dict_of_dfs[ddg_results_key].total \
                    for ddg_results_key in calculation_results_loader.ddG_cartesian_results_dict_of_dfs]
        variant_isoform_summary['ddG Cartesian Max'] = max(ddG_list)
        variant_isoform_summary['ddG Cartesian Min'] = min(ddG_list)

    def nan_to_None(x: np.float64):
        return None if np.isnan(x) else x

    if calculation_results_loader.pathprox_disease1_results_dict_of_dfs:
        pathprox_list = [
            calculation_results_loader.pathprox_disease1_results_dict_of_dfs[pathprox_results_key].iloc[0].pathprox \
            for pathprox_results_key in calculation_results_loader.pathprox_disease1_results_dict_of_dfs]
        # Get rid of Nans
        variant_isoform_summary['disease1_pp Max'] = nan_to_None(np.nanmax(np.array(pathprox_list)))
        variant_isoform_summary['disease1_pp Min'] = nan_to_None(np.nanmin(np.array(pathprox_list)))

    if calculation_results_loader.pathprox_disease2_results_dict_of_dfs:
        pathprox_list = [
            calculation_results_loader.pathprox_disease2_results_dict_of_dfs[pathprox_results_key].iloc[0].pathprox \
            for pathprox_results_key in calculation_results_loader.pathprox_disease2_results_dict_of_dfs]

        variant_isoform_summary['disease2_pp Max'] = nan_to_None(np.nanmax(np.array(pathprox_list)))
        variant_isoform_summary['disease2_pp Min'] = nan_to_None(np.nanmin(np.array(pathprox_list)))

    return variant_isoform_summary

    env = Environment(loader=FileSystemLoader(os.path.join(os.path.dirname(os.path.realpath(__file__)))))
    isoform_variant_template = env.get_template("psb_report.html")

    if LOGGER.isEnabledFor(logging.DEBUG):
        pp = pprint.PrettyPrinter(indent=1)
        LOGGER.debug("Dictionary template_vars to render:\n%s" % pp.pformat(template_vars))

    html_out = isoform_variant_template.render(template_vars)

    templatePfamGraphics = env.get_template("html/pfamGraphicsIframeTemplate.html")
    htmlPfamGraphics = templatePfamGraphics.render(template_vars)

    save_cwd = os.getcwd()
    os.chdir(variant_report_directory)

    # WE ARE NOW OPERATING FROM ..../UDN/CaseName target directory
    LOGGER.info("Now working in %s", variant_report_directory)

    html_fname = base_fname % "html"
    website_filelist.append(os.path.join(web_dir, html_fname))
    with open(html_fname, "w") as html_f:
        html_f.write(html_out)

    with open(pfamGraphicsIframe_fname, "w") as html_f:
        html_f.write(htmlPfamGraphics)

    """write_pdfs = None
    if write_pdfs:
        pdf_fname = base_fname % "pdf"
        wkhtmltopdf_fname = base_fname % "wkhtml.pdf"
        print("\nWriting final reports to %s directory:\n  %s\n  %s\n  %s\n  %s\n" % (
            variant_report_directory, pdf_fname, wkhtmltopdf_fname, html_fname, pfamGraphicsIframe_fname))

        LOGGER.warning("Temporarily disabling all logging prior to calling weasyprint write_pdf()")

        logging.disable(sys.maxsize)
        HTML(string=html_out).write_pdf(pdf_fname, stylesheets=[
            "../html/css/typography.css"])  # ,CSS(string="@page { size: A3 landscape; }")])
        logging.disable(logging.NOTSET)
        LOGGER.warning("weasyprint write_pdf() has returned.  Restoring logging")

        # Write out another .pdf using the wkhtmltopdf tool

        wkhtml_command_list = ["wkhtmltopdf", "--print-media-type", "--no-outline", "--minimum-font-size", "12",
                               html_fname, wkhtmltopdf_fname]
        process = subprocess.Popen(wkhtml_command_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (output, err) = process.communicate()

        # So unfortunately, the stderr from wkhtmltopdf is full of ==================
        # things, which are a mess for the LOGGER.
        # Get rid of stuff that ends in \r  Keep all other lines and ending in \r
        err_nl_array = err.split('\n')
        err_legit_list = []
        for ena in err_nl_array:
            crs = ena.split('\r')
            if len(crs) > 0:  ## Pick off the last string, which means all ending with Carriage Return are ignored
                cr_last = crs[-1].strip()
                if cr_last:
                    err_legit_list.append(cr_last)

        err_legit = '\n'.join(err_legit_list)

        # ContentNotFound really just info - not warning
        if (process.returncode != 0) and ("Exit with code 1 due to network error: ContentNotFoundError" not in err):
            print("Unable to complete %s exitstatus=%d due to error %s\n  output: %s\n" % (
                str(wkhtml_command_list), process.returncode, err_legit, output))
            LOGGER.warning("wkhtmltopdf stderr:\n%s", err_legit)
        else:
            LOGGER.info("wkhtmltopdf stderr:\n%s", err_legit)
    """

    # WE HAVE NOW RETURNED TO the psb_pipeline/bin directory
    os.chdir(save_cwd)
    # Close out the local log file for this mutation
    _end_local_logging(local_logger_fh)
    gathered_info['variant_report_directory'] = os.path.basename(variant_report_directory)
    gathered_info['html_fname'] = html_fname
    return gathered_info  # End of function report_one_variant_one_isoform()


# =============================================================================
# Main logic here.  A period in the argument means the user wants to launch one mutation only,
# Whether one mutation, or the more typical set, we must have html/css/javascript statuc
# resources property installed in the destination tree
# =============================================================================

PDBMapSQLdb.set_access_dictionary(config_dict)

copy_html_css_javascript()

# directly from a single mutation output file of psb_plan.py
if oneMutationOnly:
    if args.slurm:
        LOGGER.warning("--slurm option is ignored when only one report is requested")
    template_vars = report_one_variant_one_isoform(args.projectORstructures, args.workstatus)
    if template_vars:  # The usual case, we had some Pathprox outputs
        print("Single mutation report saved to %s/.pdf/.wkhtml.pdf" % template_vars['html_fname'])
    else:
        print("Due to lack of pathprox outputs, no html (or pdf) reports were created from %s" % args.workstatus)
else:
    udn_csv_filename = os.path.join(collaboration_dir,
                                    "%s_missense.csv" % args.projectORstructures)  # The argument is an entire project UDN124356
    msg = "Retrieving project mutations from %s" % udn_csv_filename
    if not infoLogging:
        print(msg)
    LOGGER.info(msg)
    df_all_mutations = pd.read_csv(udn_csv_filename, sep=',', index_col=None, keep_default_na=False, encoding='utf8',
                                   comment='#', skipinitialspace=True)
    df_all_mutations.fillna('NA', inplace=True)

    if args.slurm:
        slurm_directory = os.path.join(collaboration_dir, "slurm")
        slurm_file = os.path.join(slurm_directory, "psb_reps.slurm")
        slurm_stdout_directory = os.path.join(slurm_directory, "stdout")
        if not os.path.exists(slurm_directory):  # make the containing directory if needed
            os.makedirs(slurm_directory)
        if not os.path.exists(slurm_stdout_directory):  # make the stdout directory if needed
            os.makedirs(slurm_stdout_directory)
        jobCount = len(df_all_mutations)
        msg = "Slurm script to run %d reports for %s is %s" % (jobCount, args.projectORstructures, slurm_file)
    else:  # not a slurm run - so just do one report after another- normal case
        msg = "Reporting on all %d project %s mutations" % (len(df_all_mutations), args.projectORstructures)
    if not infoLogging:
        print(msg)
    LOGGER.info(msg)

    slurm_array = []

    # Determine the formatting style of the main case page.
    # If CHROMs and POSs are availabine in the source data, THEN we create a heirarchical grouping above the
    # individual transcript processing.
    genome_variants = defaultdict(list)
    # mutation_summaries = [] # Used to populate old format report, one line per gene
    genome_headers = []  # Used for new format report, with one line per gene, expandable to transcript isoforms


    def fetch_gene_id(uniprot_id: str):

        query = """\
SELECT ID FROM Idmapping 
where Id_Type = 'GeneID' and unp = %(unp)s"""

        with PDBMapSQLdb(config_dict) as sql_connection:
            unp_base = uniprot_id.split('-')[0]
            sql_connection.execute(query, {'unp': unp_base})
            row = sql_connection.fetchone()
            if not row or len(row) != 1:
                LOGGER.warning("unp %s is invalid, or lacks a GeneID entry in the idmapping file" % unp_base)
                return ""
            else:
                return row[0].strip()


    if {'chrom', 'pos', 'change'}.issubset(df_all_mutations.columns) and not args.slurm:
        for index, row in df_all_mutations.iterrows():
            chrom_pos_change = (row['chrom'], row['pos'], row['change'])
            LOGGER.info("%s", chrom_pos_change)
            genome_variants[chrom_pos_change].append((index, row))

        for chrom_pos_change in genome_variants:
            # Crude - but get at the original missense row
            # under which our protein transcript isoform are grouped
            first_row = genome_variants[chrom_pos_change][0][1]
            msg = "Reporting on: %s %s" % (first_row['gene'], chrom_pos_change)
            if not infoLogging:
                print(msg)
            LOGGER.info(msg)
            genome_header = {
                'Gene': first_row['gene'],
                'Chrom': first_row['chrom'],
                'Pos': first_row['pos'],
                'Change': first_row['change'],
                'disease1_pp Min': None,
                'disease1_pp Max': None,
                'disease2_pp Min': None,
                'disease2_pp Max': None,
                'ddG Min': None,
                'ddG Max': None,
                'variant_isoform_summaries': {}
            }

            variant_isoform_summaries = []
            for index, genome_variant_row in genome_variants[chrom_pos_change]:
                msg = "%s %-10s %-10s %-6s" % (
                    "      ", genome_variant_row['refseq'], genome_variant_row['mutation'], genome_variant_row['unp'])
                if not infoLogging:
                    print(msg)
                LOGGER.info(msg)
                variant_directory = "%s_%s_%s" % (
                    genome_variant_row['gene'], genome_variant_row['refseq'], genome_variant_row['mutation'])
                variant_isoform_summary = report_one_variant_one_isoform(variant_directory, genome_variant_row)
                if variant_isoform_summary:
                    variant_isoform_summary['#'] = index
                    if not variant_isoform_summary['Error']:
                        msg = "%s %s %s report saved to %s" % (
                            genome_variant_row['gene'], genome_variant_row['refseq'], genome_variant_row['mutation'],
                            variant_isoform_summary['html_filename'])
                        if not infoLogging:
                            print(msg)
                        LOGGER.info(msg)
                    variant_isoform_summary['Gene'] = genome_variant_row['gene']

                    variant_isoform_summary['Unp'] = genome_variant_row['unp']
                    variant_isoform_summary['gene_id'] = fetch_gene_id(variant_isoform_summary['Unp'])
                    transcript_list = genome_variant_row['transcript'].split(';')
                    variant_isoform_summary['Transcript'] = transcript_list[0]
                    if len(transcript_list) > 1:
                        variant_isoform_summary['Transcript'] += '...'
                    variant_isoform_summary['Refseq'] = genome_variant_row['refseq']
                    variant_isoform_summary['Mutation'] = genome_variant_row['mutation']
                    variant_isoform_summaries.append(variant_isoform_summary)
                else:
                    msg = "Due to lack of pathprox or ddG outputs, %s %s %s has no html (or pdf) report" % (
                        genome_variant_row['gene'], genome_variant_row['refseq'], genome_variant_row['mutation'])
                    if not infoLogging:
                        print(msg)
                    LOGGER.info(msg)
            genome_header['variant_isoform_summaries'] = variant_isoform_summaries

            # Go to SQL to get the Gene ID for this gene.....
            genome_header['gene_id'] = variant_isoform_summary['gene_id']

            genome_header['ddG Monomer Min'] = min(
                [variant_isoform_summary['ddG Monomer Min'] for variant_isoform_summary in variant_isoform_summaries \
                 if 'ddG Monomer Min' in variant_isoform_summary and variant_isoform_summary[
                     'ddG Monomer Min'] is not None],
                default=None)
            genome_header['ddG Monomer Max'] = max(
                [variant_isoform_summary['ddG Monomer Max'] for variant_isoform_summary in variant_isoform_summaries \
                 if 'ddG Monomer Max' in variant_isoform_summary and variant_isoform_summary[
                     'ddG Monomer Max'] is not None],
                default=None)

            genome_header['ddG Cartesian Min'] = min(
                [variant_isoform_summary['ddG Cartesian Min'] for variant_isoform_summary in variant_isoform_summaries \
                 if 'ddG Cartesian Min' in variant_isoform_summary and variant_isoform_summary[
                     'ddG Cartesian Min'] is not None],
                default=None)
            genome_header['ddG Cartesian Max'] = max(
                [variant_isoform_summary['ddG Cartesian Max'] for variant_isoform_summary in variant_isoform_summaries \
                 if 'ddG Cartesian Max' in variant_isoform_summary and variant_isoform_summary[
                     'ddG Cartesian Max'] is not None],
                default=None)

            genome_header['AA_Len Min'] = min(
                [variant_isoform_summary['AA_len'] for variant_isoform_summary in variant_isoform_summaries \
                 if 'AA_len' in variant_isoform_summary and variant_isoform_summary['AA_len'] is not None],
                default=None)

            genome_header['AA_Len Max'] = max(
                [variant_isoform_summary['AA_len'] for variant_isoform_summary in variant_isoform_summaries \
                 if 'AA_len' in variant_isoform_summary and variant_isoform_summary['AA_len'] is not None],
                default=None)

            genome_header['disease1_pp Min'] = min(
                [variant_isoform_summary['disease1_pp Min'] for variant_isoform_summary in variant_isoform_summaries \
                 if 'disease1_pp Min' in variant_isoform_summary and variant_isoform_summary[
                     'disease1_pp Min'] is not None],
                default=None)
            genome_header['disease1_pp Max'] = max(
                [variant_isoform_summary['disease1_pp Max'] for variant_isoform_summary in variant_isoform_summaries \
                 if 'disease1_pp Max' in variant_isoform_summary and variant_isoform_summary[
                     'disease1_pp Max'] is not None],
                default=None)

            genome_header['disease2_pp Min'] = min(
                [variant_isoform_summary['disease2_pp Min'] for variant_isoform_summary in variant_isoform_summaries \
                 if 'disease2_pp Min' in variant_isoform_summary and variant_isoform_summary[
                     'disease2_pp Min'] is not None],
                default=None)
            genome_header['disease2_pp Max'] = max(
                [variant_isoform_summary['disease2_pp Max'] for variant_isoform_summary in variant_isoform_summaries \
                 if 'disease2_pp Max' in variant_isoform_summary and variant_isoform_summary[
                     'disease2_pp Max'] is not None],
                default=None)

            genome_header['Error'] = max(
                [variant_isoform_summary['Error'] for variant_isoform_summary in variant_isoform_summaries \
                 if 'Error' in variant_isoform_summary and variant_isoform_summary['Error'] is not None],
                default=None)

            genome_headers.append(genome_header)

    else:  # There are no chrom positions in this case.  Use the old format
        variant_isoform_summaries = []  # One line per variant
        for index, variant_row in df_all_mutations.iterrows():
            msg = "%s %-10s %-10s %-6s" % (
                "Generating slurm entry for" if args.slurm else "Reporting on", variant_row['gene'],
                variant_row['refseq'], variant_row['mutation'])
            if not infoLogging:
                print(msg)
            LOGGER.info(msg)

            if 'RefSeqNotFound_UsingGeneOnly' in variant_row['refseq']:
                variant_row['refseq'] = 'NA'
            # gene_refseq_mutation,structure_report_filename,workstatus_filename,dropped_structure_filenames = paths_for_variant_row(row)

            if args.slurm:
                slurm_array.append("./psb_rep.py -c %s -u %s %s %s" % (
                    args.config, args.userconfig, structure_report_filename, workstatus_filename))
            else:
                variant_directory = "%s_%s_%s" % (variant_row['gene'], variant_row['refseq'], variant_row['mutation'])
                variant_isoform_summary = report_one_variant_one_isoform(variant_directory, variant_row)
                if variant_isoform_summary:
                    variant_isoform_summary['#'] = index
                    if not variant_isoform_summary['Error']:
                        msg = "%s %s %s report saved to %s" % (
                            variant_row['gene'], variant_row['refseq'], variant_row['mutation'],
                            variant_isoform_summary['html_filename'])
                        if not infoLogging:
                            print(msg)
                        LOGGER.info(msg)
                    variant_isoform_summary['Gene'] = variant_row['gene']
                    variant_isoform_summary['Unp'] = variant_row['unp']
                    variant_isoform_summary['gene_id'] = fetch_gene_id(variant_isoform_summary['Unp'])
                    variant_isoform_summary['Refseq'] = variant_row['refseq']
                    variant_isoform_summary['Mutation'] = variant_row['mutation']
                    variant_isoform_summaries.append(variant_isoform_summary)
                else:
                    msg = "Due to lack of pathprox or ddG outputs, %s %s %s has no html (or pdf) report" % (
                        variant_row['gene'], variant_row['refseq'], variant_row['mutation'])
                    if not infoLogging:
                        print(msg)
                    LOGGER.info(msg)

        if args.slurm:
            with open(slurm_file, "w") as slurmf:
                slurmf.write("""\
    #!/bin/sh
    #
    # Project        : PSB Pipeline
    # Filename       : %s
    # Generated By   : %s
    # For case       : %s
    # Organization   : Vanderbilt Genetics Institute,
    #                : Program in Personalized Structural Biology,
    #                : Department of Biomedical Informatics,
    #                : Vanderbilt University
    # Generated on   : %s
    # Description    : Runs a python script on a cluster to launch psb_rep.py for 
    #                : each and all mutations in parallel
    #===============================================================================
    # Slurm Parameters
    """ % (slurm_file, __file__, args.projectORstructures, str(datetime.datetime.now())))

                slurmDict = slurmParametersReport
                slurmDict['output'] = os.path.join(slurm_stdout_directory,
                                                   "%s_%%A_%%a_psb_rep.out" % args.projectORstructures)
                slurmDict['job-name'] = "%s_psb_reps" % args.projectORstructures

                jobCount = len(slurm_array)
                if jobCount > 1:
                    slurmDict['array'] = "0-%d" % (jobCount - 1)

                # We can set a max on the number of simultaneously running jobs in the array.  For now, don't do this
                # if jobCount > 10:
                #   slurmDict['array'] += "%%10"

                for slurm_field in ['job-name', 'mail-user', 'mail-type', 'ntasks', 'time', 'mem', 'account', 'output',
                                    'array']:
                    if slurm_field in slurmDict:
                        slurmf.write("#SBATCH --%s=%s\n" % (slurm_field, slurmDict[slurm_field]))

                if jobCount > 1:
                    slurmf.write("""
    echo "SLURM_ARRAY_TASKID="$SLURM_ARRAY_TASKID
    """)
                slurmf.write("""
    echo "SLURM_JOBID="$SLURM_JOBID
    echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
    echo "SLURM_NNODES"=$SLURM_NNODES
    # echo "SLURMTMPDIR="$SLURMTMPDIR
    echo "SLURM_SUBMIT_DIR = "$SLURM_SUBMIT_DIR
    
    """)

                cwd = os.getcwd()

                slurmf.write("""
    cd %s
    if [ $? != 0 ]; then
    echo Failure at script launch: Unable to change to directory %s
    exit 1
    fi
    """ % (cwd, cwd))

                if jobCount == 1:
                    # No need to fiddle with the slurm case statement
                    slurmf.write(slurm_array[0])
                else:
                    slurmf.write("\ncase $SLURM_ARRAY_TASK_ID in\n")
                    slurm_array_id = 0
                    for launchstring in slurm_array:
                        # Write out the "case" tag
                        slurmf.write("%d)\n" % slurm_array_id)
                        # Save this to put back later
                        slurm_array_id += 1
                        slurmf.write("echo Command: '%s'\n" % launchstring)
                        slurmf.write(launchstring)
                        # end the case tag
                        slurmf.write("\n;;\n\n")
                    slurmf.write("esac\n")
            msg = "Created slurm script to launch all psb_rep.py processes for this case: %s" % slurm_file
            if args.verbose:
                LOGGER.info(msg)
            else:
                print(msg)

    if genome_headers or variant_isoform_summaries:  # Excellent, we have a home page report for many variants
        # If there are DigenicInteraction .svg files in the right place, integrate them
        # digenic_graphics_types = [
        #     'digenic_score', 
        #     # 'digenic_details', 
        #     'sys_bio_details']
        # digenic_graphics_score_only_filename = None
        # found_graphics_files = {}
        # for digenic_graphic_type in digenic_graphics_types:
        #     digenic_graphic_filename = os.path.join('casewide', 'DiGePred',
        #                                             '%s_all_gene_pairs_%s.svg' % (
        #                                             args.projectORstructures, digenic_graphic_type))
        #
        #      if os.path.exists(os.path.join(collaboration_dir, digenic_graphic_filename)):
        #         found_graphics_files[digenic_graphic_type] = digenic_graphic_filename
        #         website_filelist.append(os.path.join('.', digenic_graphic_filename))

        # Was there a generated .html file from DiGePred
        digepred_html_filename = os.path.join('casewide', 'DiGePred',
                                              '%s_all_gene_pairs_summary.html' % (
                                                  args.projectORstructures,))

        if os.path.exists(os.path.join(collaboration_dir, digepred_html_filename)):
            website_filelist.append(digepred_html_filename)
            LOGGER.info("Integrating DiGePred html: %s", digepred_html_filename)
        else:
            LOGGER.warning("DiGePred %s missing.  Did you run DiGePred?", digepred_html_filename)
            digepred_html_filename = None

        digepred_gene_pairs = []
        if digepred_html_filename:
            digepred_csv_filename = os.path.join('casewide', 'DiGePred',
                                                 '%s_all_gene_pairs_digenic_metrics.csv' % (
                                                     args.projectORstructures,))

            digepred_metrics_df = pd.DataFrame()
            if os.path.exists(os.path.join(collaboration_dir, digepred_csv_filename)):
                LOGGER.info("Integrating DiGePred csv: %s", digepred_csv_filename)
                digepred_metrics_df = pd.read_csv(digepred_csv_filename, sep=',')
                LOGGER.info("%d rows read", len(digepred_metrics_df))
                digepred_metrics_df.sort_values(by=['digenic score'], ascending=False, inplace=True, ignore_index=True)
            else:
                LOGGER.critical("DiGePred %s missing.  Alarming, as .html file was found", digepred_csv_filename)
                del digepred_csv_filename

            # Now prepare the entries which will go to the html table
            # For now, max 20 or all
            max_rows = min(20, len(digepred_metrics_df))
            for index, row in digepred_metrics_df.head(max_rows).iterrows():
                gene_pair_dict = {}
                for key in ['gene A', 'gene B', 'digenic score']:
                    gene_pair_dict[key] = row[key]
                digepred_gene_pairs.append(gene_pair_dict)

        case_report_template_location = os.path.dirname(os.path.realpath(__file__))
        env = Environment(loader=FileSystemLoader(case_report_template_location))
        # if len(found_graphics_files) == 2:
        #     LOGGER.info("Integrating DiGePred svg graphics")
        #     digenic_graphics_score_only_filename = os.path.join('.', found_graphics_files['digenic_score'])
        #     template = env.get_template("html/DigenicInteractionsReportTemplate.html")
        #     # Probably a goof - but the template file restates the full graphics filenames
        #     html_out = template.render({'case': args.projectORstructures})
        #     digenic_graphics_html_filename = os.path.join(collaboration_dir, 'casewide',
        #         #         #         #         #         #         #   'DigenicGraphics.html')  # The argument is an entire project UDN124356
        #     with open(digenic_graphics_html_filename, "w") as f:
        #         # f.write(html_out)
        #     website_filelist.append(digenic_graphics_html_filename)
        # else:
        #     LOGGER.warning("DiGePred svg files are missing.  Did you run DiGePred")

        # pprint.pformat(mutation_summaries)
        # Grab Souhrids gene interaction information
        # geneInteractions_generic, html_table_generic, html_report_generic, text_report_generic = \
        #     gene_interaction_report(collaboration_dir, args.projectORstructures, False)
        # geneInteractions_familial, html_table_familial, html_report_familial, text_report_familial = \
        #     gene_interaction_report(collaboration_dir, args.projectORstructures, True)
        # if geneInteractions_generic:
        #     for genome_header in genome_headers:
        #         if genome_header['Gene'] in geneInteractions_generic:
        #             genome_header['Gene Interactions Generic'] = ' '.join(
        #                 geneInteractions_generic[genome_header['Gene']])
        # if geneInteractions_familial:
        #     for genome_header in genome_headers:
        #         if genome_header['Gene'] in geneInteractions_familial:
        #             genome_header['Gene Interactions Familial'] = ' '.join(
        #                 geneInteractions_familial[genome_header['Gene']])
        case_report_template = env.get_template("case_report_template.html")
        # print html_table
        final_gathered_info = {'variant_isoform_summaries': variant_isoform_summaries,
                               'genome_headers': genome_headers,
                               # 'firstGeneTable': html_table_generic,
                               # 'firstGeneReport': html_report_generic,
                               # 'secondGeneTable': html_table_familial,
                               # 'secondGeneReport': html_report_familial,
                               'case': args.projectORstructures,
                               "date": time.strftime("%Y-%m-%d"),
                               'disease1_variant_short_description': config_pathprox_dict[
                                   'disease1_variant_short_description'],
                               'disease2_variant_short_description': config_pathprox_dict[
                                   'disease2_variant_short_description'],
                               'digepred_html_filename': digepred_html_filename,
                               'digepred_gene_pairs': digepred_gene_pairs
                               }
        if LOGGER.isEnabledFor(logging.DEBUG):
            pp = pprint.PrettyPrinter(indent=1)
            LOGGER.debug("Dictionary final_gathered_info to render:\n%s" % pp.pformat(final_gathered_info))

        html_out = case_report_template.render(final_gathered_info)
        case_summary_filename = os.path.join(collaboration_dir,
                                             "%s.html" % args.projectORstructures)  # The argument is an entire project UDN124356
        website_filelist.append(case_summary_filename)
        with open(case_summary_filename, "w") as f:
            f.write(html_out)
        lastmsg = "The case summary report is: " + case_summary_filename
        case_summary_json = os.path.join(collaboration_dir,
                                         "%s.json" % args.projectORstructures)  # The argument is an entire project UDN124356
        with open(case_summary_json, 'w') as fp:
            json.dump(final_gathered_info, fp)
    else:
        lastmsg = "No mutation summaries - bummer"

    print("")
    print("")

    if args.verbose:
        LOGGER.info(lastmsg)
    else:
        print(lastmsg)

    if website_filelist:
        try:
            os.remove('index.html')  # First time through will gen an exception.  That's A-OK
        except:
            pass
        os.symlink(case_summary_filename, 'index.html')
        website_filelist.append('index.html')

        website_filelist_filename = os.path.join(collaboration_dir,
                                                 "%s_website_files.list" % args.projectORstructures)  # The argument is an entire project UDN124356
        with open(website_filelist_filename, "w") as f:
            f.write('\n'.join((os.path.relpath(
                website_file,  # os.path.realpath(website_file),
                (os.path.realpath(os.path.join(config_dict['output_rootdir'], config_dict['collaboration']))))
                for website_file in website_filelist)))
        LOGGER.info("A filelist for creating a website is in %s", website_filelist_filename)
        website_zip_filename = "%s.zip" % args.projectORstructures
        pkzip_maker = 'rm -f %s; cd ..; cat %s | zip -r@ %s > %s.stdout; cd -' % (website_zip_filename,
                                                                                  os.path.join(args.projectORstructures,
                                                                                               website_filelist_filename),
                                                                                  os.path.join(args.projectORstructures,
                                                                                               website_zip_filename),
                                                                                  os.path.join(args.projectORstructures,
                                                                                               website_zip_filename))
        LOGGER.info("Executing: %s", pkzip_maker)
        subprocess.call(pkzip_maker, shell=True)

        website_tar_filename = "%s.tar.gz" % args.projectORstructures
        tar_maker = 'rm -f %s; cd ..; tar cvzf %s --files-from %s --mode=\'a+rX,go-w,u+w\' > %s.stdout; cd -' % (
            website_tar_filename,
            os.path.join(args.projectORstructures, website_tar_filename),
            os.path.join(args.projectORstructures, website_filelist_filename),
            os.path.join(args.projectORstructures, website_tar_filename))
        LOGGER.info("Executing: %s", tar_maker)
        subprocess.call(tar_maker, shell=True)

        print("Compressed website files are in %s and %s" % (website_zip_filename, website_tar_filename))

sys.exit(0)
