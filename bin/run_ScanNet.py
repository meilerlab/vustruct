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
Launch ScanNet for a specific uniprot ID, and capture outputs in a manner
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
from Bio.PDB import PDBIO
from Bio import BiopythonWarning
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from psb_shared import psb_config
from psb_shared.psb_progress import PsbStatusManager

from lib import PDBMapTranscriptUniprot
from lib import PDBMapAlphaFold

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
                help="Uniprot ID for alphafold model selection")

cmdline_parser.add_argument(
                "--transcript_mutation", type=str,
                help="Transcript offset of variant of unknown significance")

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
alphafold = PDBMapAlphaFold(config_dict)

# For now, we don't have a trans_mut_pos - but we might integrate this later...
trans_mut_pos = int(args.transcript_mutation[1:-1])
if trans_mut_pos:
    model_seq_start, model_seq_end, alphafold_modelid = alphafold.best_covering_model(
        args.unp,
        len(transcript_uniprot.aa_seq),
        trans_mut_pos)
else:
    alphafold_modelid = alphafold.first_window_modelid(args.unp)
    model_seq_start = 1

alphafold_cif_filename = alphafold.get_coord_filename(alphafold_modelid)

if not os.path.exists(alphafold_cif_filename):
    LOGGER.critical("Alphafold model %s does not exist", alphafold_cif_filename)
    exit(1)

# This is so annoying - but we need to save the .cif as a .pdb to feed into scannet
cif_pos = alphafold_cif_filename.find(".cif")
# LOGGER.info("cifpos = %s" % cif_pos)
alphafold_local_pdb = alphafold_cif_filename[0:cif_pos] + ".pdb"
alphafold_local_pdb = os.path.join(args.outdir,os.path.basename(alphafold_local_pdb))

def _open_for_extension(file_name: str):
    _, _coordinates_filename_extension = os.path.splitext(file_name)
    if _coordinates_filename_extension == '.gz':
        _open_function, _open_mode = (gzip.open, 'rt')
    elif _coordinates_filename_extension == '.xz':
        _open_function, _open_mode = (lzma.open, 'rt')
    else:
        _open_function, _open_mode = (open, 'r')

    return _open_function(file_name, _open_mode)

def load_structure_and_mmcif_dict(alphafold_cif_filename, alphafold_modelid):
    # HACKED TO NOT consider PDB - just cif for alphafold
    # With self.source_coordinates_filename initialized beforehand, call this routine to load
    # a structure file in pdb or cif format.  On return, self.structure will be
    # populated with a Biopython structure.  If the source file is cif format
    # we will also init self.mmcif_dict

    mmcif_dict = {}
    structure = None

    with _open_for_extension(alphafold_cif_filename) as fin:
        warnings.filterwarnings('ignore', category=PDBConstructionWarning)
        # source_coordinates_filename_pieces = self.source_coordinates_filename.split('.')
        _mmcif_parser = MMCIFParser(QUIET=True)
        structure = _mmcif_parser.get_structure(alphafold_modelid, fin)
        mmcif_dict = _mmcif_parser._mmcif_dict
        warnings.resetwarnings()

    return [structure, mmcif_dict]


LOGGER.info("struct filename is %s %s" % (alphafold_cif_filename, alphafold_local_pdb) )

alphafold_structure, alphafold_mmcif_dict = load_structure_and_mmcif_dict(alphafold_cif_filename,alphafold_modelid)

# Get a count of residues for the display.  Assumes (always) only one model and chain in the alphafold
for model in alphafold_structure:
    for chain in model:
        break

LOGGER.info("%s loaded with %d residues" % (alphafold_cif_filename, len(chain)))

io=PDBIO()
io.set_structure(alphafold_structure)
io.save(alphafold_local_pdb)
LOGGER.info("structure saved as PDB in %s" % alphafold_local_pdb)


singularity_command_list = [
    'singularity',
    'exec',
    '/dors/capra_lab/users/mothcw/VUStruct/ScanNet.simg',
    '/bin/bash']
container_stdin_list = [
   'PATH=/usr/local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin',
   'python3 /ScanNet/predict_bindingsites.py %s --noMSA --predictions_folder %s --mode interface' % (
              alphafold_local_pdb,
              os.path.join(args.outdir,"ScanNet"))
    ]

container_stdin_list.append("exit\n")

container_stdin_str = "\n".join(container_stdin_list)

LOGGER.info("Launching container with %s " % ' '.join(singularity_command_list))
LOGGER.info("stdin=%s" % container_stdin_str)

completed_process = subprocess.run(
    singularity_command_list,
    input=container_stdin_str,
    shell=False,
    text=True, 
    # env=environment_override, 
    capture_output=True)
LOGGER.info("returncode = %d" % completed_process.returncode);
LOGGER.info("stdout = %s" % completed_process.stdout);
LOGGER.info("stderr = %s" % completed_process.stderr);
LOGGER.info("Successful end of ScanNet.py")
if completed_process.returncode == 0:
    psb_status_manager.write_info("Success")
    psb_status_manager.mark_complete()
else:
    psb_status_manager.write_info(completed_process.stderr)
    psb_status_manager.mark_failed()
exit(0)





# transcript_uniprot.load_aa_seq_from_sql()
fasta_filename = args.unp + '.fasta'
fasta_fullpath = os.path.join(args.outdir,fasta_filename)
LOGGER.info("Writing " + args.unp + " in fasta format to " + fasta_fullpath)
with open(fasta_fullpath,'w') as fasta_f:
    fasta_f.write(transcript_uniprot.fasta_aa_seq)



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




for musite_model_prefix in musite_model_prefix_list:
    container_stdin_list.append(" ".join([
        'python3',
        'predict_multi_batch.py',
         '-input', fasta_fullpath,
         '-output', os.path.join(args.outdir, args.unp + '_' + musite_model_prefix ),
         '-model-prefix', os.path.join("models", musite_model_prefix)
        ]))



    


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

LOGGER.info("Successful end of ScanNet.py")
psb_status_manager.mark_complete()
sys.exit(0)
