#!/usr/bin/env python3
#
# Project        : PSB Pipeline
# Filename       : psb_rep.py
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
import pandas as pd
from jinja2 import Environment, FileSystemLoader
from io import StringIO
import json

import subprocess
from lib import PDBMapSQLdb
from lib import PDBMapTranscriptBase
from lib import PDBMapTranscriptUniprot
from lib import PDBMapTranscriptEnsembl

from vustruct import VUstruct
from report_modules.rep_one_variant import report_one_variant_one_isoform
from report_modules.load_one_variant import CalculationResultsLoader



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

required_config_items = ['output_rootdir', 'collaboration', 'ddg_config', 'rate4site_dir', 'cosmis_dir']

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

# The case_root_dir is the master directory for the case, i.e. for one patient
# Example: /dors/capra_lab/projects/psb_collab/UDN/UDN532183
udn_root_directory = os.path.join(config_dict['output_rootdir'], config_dict['collaboration'])

if oneMutationOnly:
    collaboration_absolute_dir = os.path.dirname(os.path.dirname(args.projectORstructures))
else:
    collaboration_absolute_dir = os.path.join(udn_root_directory, args.projectORstructures)
    vustruct = VUstruct("report", args.projectORstructures, __file__)

collaboration_absolute_dir = os.path.realpath(collaboration_absolute_dir)
LOGGER.info("Collaboration_absolute_dir = %s", collaboration_absolute_dir)

case_root_dir = os.path.relpath(os.getcwd(), collaboration_absolute_dir)
LOGGER.info("Relative to cwd=%s, case_root_directory= %s", os.getcwd(), case_root_dir)


def copy_html_css_javascript():
    # The source html directory, and it's large set of supporting files, is located beneath the directory
    # where this file, psb_rep.py, is located.  All mutations share this one large html/css/javascript
    # sourcee directory
    src = os.path.join(os.path.dirname(os.path.realpath(__file__)), "html")
    dest = "%s/html" % case_root_dir

    nglviewer_filename = os.path.join(dest, "nglviewer", "examples", "dist", "ngl.js")
    if os.path.exists(nglviewer_filename):
        LOGGER.warning("%s found.  Skipping copy_html_css_javascript", nglviewer_filename)
        return
    else:
        LOGGER.info("%s not found.  Assuming first time copy", nglviewer_filename)

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

    pairs = {'direct': [], 'pathways': [], 'phenotypes': [], 'other': []}

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







# =============================================================================
# Main logic here.  A period in the argument means the user wants to launch one mutation only,
# Whether one mutation, or the more typical set, we must have html/css/javascript statuc
# resources property installed in the destination tree
# =============================================================================

PDBMapSQLdb.set_access_dictionary(config_dict)

copy_html_css_javascript()

# CREATE the log/ webpage first, which includes details of everything done so far...
case_report_template_location = os.path.dirname(os.path.realpath(__file__))
jinja2_environment = Environment(loader=FileSystemLoader(case_report_template_location))
vustruct_logs_template = jinja2_environment.get_template("vustruct_logs.html")
vustruct_dict_for_jinja2 = vustruct.dict_for_jinja2()
# We need to get rid of log/ prefixes, because we generate index.html
# down in the log directory
for command_line_module in vustruct_dict_for_jinja2.keys():
    if 'logfile' in vustruct_dict_for_jinja2[command_line_module]:
        vustruct_dict_for_jinja2[command_line_module]['logfile'] = \
            os.path.basename(vustruct_dict_for_jinja2[command_line_module]['logfile'])

vustruct_logs_info = {
    'vustruct': vustruct_dict_for_jinja2
}

html_out = vustruct_logs_template.render(vustruct_logs_info)
vustruct_html_filename = os.path.join("log/", "index.html")
website_filelist.append(vustruct_html_filename)
with open(vustruct_html_filename, "w") as f:
    f.write(html_out)

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
    udn_csv_filename = os.path.join(case_root_dir,
                                    "%s_missense.csv" % args.projectORstructures)  # The argument is an entire project UDN124356
    msg = "Retrieving project mutations from %s" % udn_csv_filename
    if not infoLogging:
        print(msg)
    LOGGER.info(msg)
    df_all_mutations = pd.read_csv(udn_csv_filename, sep=',', index_col=None, keep_default_na=False, encoding='utf8',
                                   comment='#', skipinitialspace=True)
    df_all_mutations.fillna('NA', inplace=True)

    if args.slurm:
        slurm_directory = os.path.join(case_root_dir, "slurm")
        slurm_file = os.path.join(slurm_directory, "psb_reps.slurm")
        slurm_stdout_directory = os.path.join(slurm_directory, "stdout")
        if not os.path.exists(slurm_directory):  # make the containing directory if needed
            os.makedirs(slurm_directory)
        if not os.path.exists(slurm_stdout_directory):  # make the stdout directory if needed
            os.makedirs(slurm_stdout_directory)
        jobCount = len(df_all_mutations)
        msg = "Slurm script to run %d reports for %s is %s" % (jobCount, args.projectORstructures, slurm_file)
    else:  # not a slurm run - so just do one report after another- normal case
        msg = "Reporting on all %d case=%s mutations" % (len(df_all_mutations), args.projectORstructures)
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


    # When we have chrom+pos+change columns in the input, then make collapsible rows
    # for each entry, where the transcripts are shown when the row is expanded.
    if {'chrom', 'pos', 'change'}.issubset(df_all_mutations.columns) and not args.slurm:
        for index, row in df_all_mutations.iterrows():
            chrom_pos_change = (row['chrom'], row['pos'], row['change'])
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
            genome_header['variant_isoform_summaries'] = variant_isoform_summaries

            for index, genome_variant_row in genome_variants[chrom_pos_change]:
                msg = "%s %-10s %-10s %-6s" % (
                    "      ", genome_variant_row['refseq'], genome_variant_row['mutation'], genome_variant_row['unp'])
                if not infoLogging:
                    print(msg)
                LOGGER.info(msg)
                variant_directory = "%s_%s_%s" % (
                    genome_variant_row['gene'], genome_variant_row['refseq'], genome_variant_row['mutation'])
                variant_isoform_summary, additional_website_files = report_one_variant_one_isoform(
                    case_root_dir= case_root_dir,
                    variant_directory_segment= variant_directory,
                    parent_report_row=genome_variant_row,
                    local_log_level= "debug" if args.debug else "info" if args.verbose else "warn",
                    vustruct_pipeline_launched= bool(vustruct_dict_for_jinja2['launch']['executable']),
                    config_PathProx_dict= config.items('PathProx')
                )
                if variant_isoform_summary:
                    website_filelist.extend(additional_website_files)
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

                    # Go to SQL to get the Gene ID for this gene.....
                    genome_header['gene_id'] = variant_isoform_summary['gene_id']

                    genome_header['ddG Monomer Min'] = min(
                        [variant_isoform_summary['ddG Monomer Min'] for variant_isoform_summary in
                         variant_isoform_summaries \
                         if 'ddG Monomer Min' in variant_isoform_summary and variant_isoform_summary[
                             'ddG Monomer Min'] is not None],
                        default=None)
                    genome_header['ddG Monomer Max'] = max(
                        [variant_isoform_summary['ddG Monomer Max'] for variant_isoform_summary in
                         variant_isoform_summaries \
                         if 'ddG Monomer Max' in variant_isoform_summary and variant_isoform_summary[
                             'ddG Monomer Max'] is not None],
                        default=None)

                    genome_header['ddG Cartesian Min'] = min(
                        [variant_isoform_summary['ddG Cartesian Min'] for variant_isoform_summary in
                         variant_isoform_summaries \
                         if 'ddG Cartesian Min' in variant_isoform_summary and variant_isoform_summary[
                             'ddG Cartesian Min'] is not None],
                        default=None)
                    genome_header['ddG Cartesian Max'] = max(
                        [variant_isoform_summary['ddG Cartesian Max'] for variant_isoform_summary in
                         variant_isoform_summaries \
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
                        [variant_isoform_summary['disease1_pp Min'] for variant_isoform_summary in
                         variant_isoform_summaries \
                         if 'disease1_pp Min' in variant_isoform_summary and variant_isoform_summary[
                             'disease1_pp Min'] is not None],
                        default=None)
                    genome_header['disease1_pp Max'] = max(
                        [variant_isoform_summary['disease1_pp Max'] for variant_isoform_summary in
                         variant_isoform_summaries \
                         if 'disease1_pp Max' in variant_isoform_summary and variant_isoform_summary[
                             'disease1_pp Max'] is not None],
                        default=None)

                    genome_header['disease2_pp Min'] = min(
                        [variant_isoform_summary['disease2_pp Min'] for variant_isoform_summary in
                         variant_isoform_summaries \
                         if 'disease2_pp Min' in variant_isoform_summary and variant_isoform_summary[
                             'disease2_pp Min'] is not None],
                        default=None)
                    genome_header['disease2_pp Max'] = max(
                        [variant_isoform_summary['disease2_pp Max'] for variant_isoform_summary in
                         variant_isoform_summaries \
                         if 'disease2_pp Max' in variant_isoform_summary and variant_isoform_summary[
                             'disease2_pp Max'] is not None],
                        default=None)

                    genome_header['Error'] = max(
                        [variant_isoform_summary['Error'] for variant_isoform_summary in variant_isoform_summaries \
                         if 'Error' in variant_isoform_summary and variant_isoform_summary['Error'] is not None],
                        default=None)

                    genome_headers.append(genome_header)



                else:
                    msg = "Due to lack of pathprox or ddG outputs, %s %s %s has no html (or pdf) report" % (
                        genome_variant_row['gene'], genome_variant_row['refseq'], genome_variant_row['mutation'])
                    if not infoLogging:
                        print(msg)
                    LOGGER.info(msg)

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
        #      if os.path.exists(os.path.join(case_root_dir, digenic_graphic_filename)):
        #         found_graphics_files[digenic_graphic_type] = digenic_graphic_filename
        #         website_filelist.append(os.path.join('.', digenic_graphic_filename))

        # Was there a generated .html file from DiGePred
        digepred_html_filename = os.path.join('casewide', 'DiGePred',
                                              '%s_all_gene_pairs_summary.html' % (
                                                  args.projectORstructures,))

        if os.path.exists(os.path.join(case_root_dir, digepred_html_filename)):
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
            if os.path.exists(os.path.join(case_root_dir, digepred_csv_filename)):
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
        jinja2_environment = Environment(loader=FileSystemLoader(case_report_template_location))
        # if len(found_graphics_files) == 2:
        #     LOGGER.info("Integrating DiGePred svg graphics")
        #     digenic_graphics_score_only_filename = os.path.join('.', found_graphics_files['digenic_score'])
        #     template = jinja2_environment.get_template("html/DigenicInteractionsReportTemplate.html")
        #     # Probably a goof - but the template file restates the full graphics filenames
        #     html_out = template.render({'case': args.projectORstructures})
        #     digenic_graphics_html_filename = os.path.join(case_root_dir, 'casewide',
        #         #         #         #         #         #         #   'DigenicGraphics.html')  # The argument is an entire project UDN124356
        #     with open(digenic_graphics_html_filename, "w") as f:
        #         # f.write(html_out)
        #     website_filelist.append(digenic_graphics_html_filename)
        # else:
        #     LOGGER.warning("DiGePred svg files are missing.  Did you run DiGePred")

        # pprint.pformat(mutation_summaries)
        # Grab Souhrids gene interaction information
        # geneInteractions_generic, html_table_generic, html_report_generic, text_report_generic = \
        #     gene_interaction_report(case_root_dir, args.projectORstructures, False)
        # geneInteractions_familial, html_table_familial, html_report_familial, text_report_familial = \
        #     gene_interaction_report(case_root_dir, args.projectORstructures, True)
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
        case_report_template = jinja2_environment.get_template("case_report_template.html")
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

        # args.projectORstructures is an entire project UDN124356
        case_summary_filename = os.path.join(case_root_dir,
                                             "%s.html" % args.projectORstructures)
        website_filelist.append(case_summary_filename)
        with open(case_summary_filename, "w") as f:
            f.write(html_out)
        last_message = "The case summary report is: " + case_summary_filename
        case_summary_json = os.path.join(case_root_dir,
                                         "%s.json" % args.projectORstructures)
        with open(case_summary_json, 'w') as fp:
            json.dump(final_gathered_info, fp)
    else:
        last_message = "No mutation summaries - bummer"

    print("")
    print("")

    if args.verbose:
        LOGGER.info(last_message)
    else:
        print(last_message)

    if website_filelist:
        try:
            os.remove('index.html')  # First time through will gen an exception.  That's A-OK
        except OSError:  # A-OK if file is not arleady there
            pass
        try:
            # The symlink fails when we are running in some vm environments - so we just copy in those cases
            os.symlink(case_summary_filename, 'index.html')
        except PermissionError:
            LOGGER.info("Creating index.html as symlink failed in VM.  Attempting simpler cp operations")
            shutil.copy(src=case_summary_filename, dst='index.html')
            pass
        website_filelist.append('index.html')

        website_filelist_filename = os.path.join(case_root_dir,
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
