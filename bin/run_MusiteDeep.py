#!/usr/bin/env python
#
# Project        : ddg_run_cartesian
# Filename       : ddg_run_cartesian.py
# Authors        : Chris Moth, based on ddg_cartesian.py
# Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : chris.moth@vanderbilt.edu
# Date           : 2023-
# =============================================================================#

"""
Launch MusiteDeep for a specific uniprot ID, and capture outputs in a manner
compatible with containerization and pipeline reporting.
"""

import os
import time
import stat
import grp
import sys
import subprocess
import logging
from logging.handlers import RotatingFileHandler
import configparser
import gzip
import tempfile
import warnings
import datetime
import pandas as pd
from tabulate import tabulate

# from time import strftime
from typing import Dict
from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser
from Bio import BiopythonWarning
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from psb_shared import psb_config
from psb_shared.psb_progress import PsbStatusManager

from lib import PDBMapTranscriptUniprot


#=============================================================================#
## Function Definitions ##
try:
    capra_group = grp.getgrnam('capra_lab').gr_gid
except KeyError:
    capra_group = os.getegid()

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

cmdline_parser.add_argument(
                "--unp", type=str,
                help="Uniprot ID for fasta sequence load")

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
    "pdb_dir",
    "modbase2020_dir",
    "modbase2020_summary",
    "alphafold_dir",
    "modbase2020_dir",
    "swiss_dir",
    "swiss_summary"]

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

LOGGER.info("Job status directory: %s" % psb_status_manager.status_dir)

transcript_uniprot = PDBMapTranscriptUniprot(args.unp)
# transcript_uniprot.load_aa_seq_from_sql()
fasta_filename = args.unp + '.fasta'
fasta_fullpath = os.path.join(args.outdir,fasta_filename)
LOGGER.info("Writing " + args.unp + " in fasta format to " + fasta_fullpath)
with open(fasta_fullpath,'w') as fasta_f:
    fasta_f.write(transcript_uniprot.fasta_aa_seq)

singularity_command_list = [
    'singularity',
    'exec',
    '/dors/capra_lab/users/mothcw/VUStruct/musitedeep.simg',
    '/bin/bash']

LOGGER.info("Launching container with %s " % ' '.join(singularity_command_list))

musite_model_prefix_list = [
    "Hydroxylysine",
    "Hydroxyproline",
    "Methylarginine",
    "Methyllysine",
    "N6-acetyllysine",
    "N-linked_glycosylation",
    "O-linked_glycosylation",
    "Phosphoserine_Phosphothreonine",
    "Phosphotyrosine",
    "Pyrrolidone_carboxylic_acid",
    "S-palmitoyl_cysteine",
    "SUMOylation",
    "Ubiquitination"]


container_stdin_list = [
   'PATH=/opt/miniconda/bin:$PATH',
   'source /opt/miniconda/etc/profile.d/conda.sh',
   'conda activate MusiteDeep',
   'cd /MusiteDeep_web/MusiteDeep']


for musite_model_prefix in musite_model_prefix_list:
    container_stdin_list.append(" ".join([
        'python3',
        'predict_multi_batch.py',
         '-input', fasta_fullpath,
         '-output', os.path.join(args.outdir, args.unp + '_' + musite_model_prefix ),
         '-model-prefix', os.path.join("models", musite_model_prefix)
        ]))

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

    


"""
        fail_message = None
        if completed_process.returncode == 0:
            LOGGER.info("%s completed successfully (exit 0)", binary_program_basename)
        else:
            fail_message = "%s failed with exit %d" % (binary_program_basename, completed_process.returncode)
            LOGGER.critical(fail_message)
"""


# Resume by figuring out that we should keep everything or whatever
# this case of a SHEEP protein...

LOGGER.info("Successful end of MusiteDeep.py")
psb_status_manager.mark_complete()
sys.exit(0)
