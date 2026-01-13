#!/usr/bin/env python
#
# Project        : ddg_run_cartesian
# Filename       : run_DIEP.py - based on Alican's run_DIEP and plot_DIEP scripts
# Authors        : Chris Moth, based on ddg_cartesian.py
# Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : chris.moth@vanderbilt.edu
# Date           : 2023-
# =============================================================================#

"""
Launch DIEP for a case genelist, and capture outputs in a manner
compatible with containerization and pipeline reporting.
"""

import os
import time
import stat
import grp
import sys
from io import StringIO
import subprocess
import logging
from logging.handlers import RotatingFileHandler
import configparser
import gzip
import tempfile
import warnings
import datetime
import pandas as pd
import matplotlib
# Force matplotlib to not use any Xwindows backend
matplotlib.use('Agg')
import seaborn

# from time import strftime
from typing import Dict

from psb_shared import psb_config
from psb_shared.psb_progress import PsbStatusManager

# from lib import PDBMapTranscriptUniprot
# from lib import PDBMapAlphaFold

# Initiate logging
ch = logging.StreamHandler()
LOGGER = logging.getLogger()
LOGGER.addHandler(ch)

log_format_string = '%(asctime)s %(levelname)-4s [%(filename)16s:%(lineno)d] %(message)s'
date_format_string = '%H:%M:%S'
log_formatter = logging.Formatter(log_format_string, date_format_string)
ch.setFormatter(log_formatter)

LOGGER.setLevel(logging.DEBUG)
ch.setLevel(logging.INFO)

cmdline_parser = psb_config.create_default_argument_parser(__doc__, os.path.dirname(os.path.dirname(__file__)),
                                                           ".",
                                                           add_help=True)

cmdline_parser.add_argument("-o", "--outdir", type=str,
                            help="Directory to use for output and results")
cmdline_parser.add_argument("--uniquekey", type=str, required=False,
                            help="A gene/refseq/mutation/structure/chain/flavor unique identifer for this pathprox job")
cmdline_parser.add_argument("--genesfile", type=str, required=True,
                            help="A comma delimited file of gene names to process in DIEPi.  The first row header must include gene as a column")
cmdline_parser.add_argument("--collate_only", required=False, action='store_true', default=False,
                            help="Skip the calling of DIEP, and jump directly to processing output files")

args = cmdline_parser.parse_args()

# cache_dir = "/tmp/cache"
log_filename = os.path.join(args.outdir, args.uniquekey + ".log")

LOGGER.info("Log file is %s" % log_filename)
need_roll = os.path.isfile(log_filename)
LOGGER.info("Command: %s" % ' '.join(sys.argv))

if args.debug:
    LOGGER.setLevel(logging.DEBUG)
    ch.setLevel(logging.DEBUG)

required_config_items = [
    "singularity_images_dir"
    ]

config, config_dict = psb_config.read_config_files(args, required_config_items)

config_dict_reduced = {x: config_dict[x] for x in required_config_items}

config_dict = config_dict_reduced

def _open(self):
    """
    Open the current base file with the (original) mode and encoding.
    Return the resulting stream.
    """
    open_func = self._builtin_open
    return open_func(self.baseFilename, self.mode,
                     encoding=self.encoding, errors=self.errors)

logger_fh = RotatingFileHandler(log_filename, backupCount=5)

logger_fh.setFormatter(log_formatter)
logger_fh.setLevel(logging.DEBUG if args.debug else logging.INFO)
LOGGER.addHandler(logger_fh)

if need_roll:
    logger_fh.doRollover()

psb_status_manager = PsbStatusManager(args.outdir)

genes_df = pd.read_csv(args.genesfile,sep=',')
if not 'gene' in genes_df.columns:
    LOGGER.error("%s does not include a gene column.  Aborting", args.genesfile)
    sys.exit(0)

gene_list_with_possible_duplicates = genes_df['gene'].tolist()
gene_set = set()
gene_list = []
for gene in gene_list_with_possible_duplicates:
    if gene not in gene_set:
        gene_set.add(gene)
        gene_list.append(gene)



singularity_command_list = [
    'apptainer',
    'exec',
    os.path.join(config_dict['singularity_images_dir'],'DIEP.simg'])
    '/bin/bash']

LOGGER.info("Job status directory: %s" % psb_status_manager.status_dir)
DIEP_outfile = os.path.join(args.outdir,"DIEP_scores.csv")
if args.collate_only:
    LOGGER.info("Skipping containerized calls to DIEP, assuming already run")
else: # We're running musite deep as usual
    LOGGER.info("Launching container with %s " % ' '.join(singularity_command_list))
    container_stdin_list = []
    for i in range(0,len(gene_list)):
        for j in range(0,len(gene_list)):
            if i != j:
                gene_i = gene_list[i]
                gene_j = gene_list[j]
                if gene_i != gene_j:
                    container_stdin_list.append("/bin/java -jar /transfer.jar --read /Coding_predict_fixed/Coding_predict_fixed.txt -gn %s %s | grep -A1 GeneA | grep -v GeneA" % (
                        gene_i, gene_j))

    container_stdin_list.append("exit\n")

    container_stdin_str = "\n".join(container_stdin_list)

    LOGGER.info("stdin=%s" % container_stdin_str)

    completed_process = subprocess.run(
        singularity_command_list,
        input=container_stdin_str,
        shell=False,
        text=True,
        # env=environment_override, 
        capture_output=True)

    LOGGER.info("stdout = %s" % completed_process.stdout);
    LOGGER.info("stderr = %s" % completed_process.stderr);
    LOGGER.info("Writing output to %s" % DIEP_outfile)
    with open(DIEP_outfile,'w') as DIEP_f:
        DIEP_f.write(completed_process.stdout)
    LOGGER.info("Successful end of run_DIEP.py")
    psb_status_manager.mark_complete()

# Whether we collate or not, we need to create the .png file


DIEP_df = pd.read_csv(DIEP_outfile, sep='\t', names=['X', 'Y', 'Z'])
LOGGER.info("DIEP graphic will be created from\n%s" % DIEP_df)

data_pivoted = DIEP_df.pivot_table(index='X', columns='Y', values='Z')

ax = seaborn.heatmap(data_pivoted, cmap='viridis', vmin=0, vmax=1.0)
matplotlib.pyplot.title('DIEP')
matplotlib.pyplot.tight_layout()
DIEP_graphics_filename = os.path.join(args.outdir,"DIEP.png")
LOGGER.info("Saving graphic to %s", DIEP_graphics_filename)
matplotlib.pyplot.savefig(DIEP_graphics_filename, dpi=600)


sys.exit(0)
