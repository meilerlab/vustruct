#!/usr/bin/env python
#
# Project        : VUStruct
# Filename       : run_PTMGP2.py
# Authors        : Chris Moth
# Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : chris.moth@vanderbilt.edu
# Date           : 2023-
# =============================================================================#

"""
Launch PTMGP2 for a specific uniprot ID, and capture outputs in a manner
compatible with containerization and pipeline reporting.

The prediction algoritm is strict sequence based.
After the predictions are made, we report all those residues with heavy atoms 
with 8 Angstroms of the mutation pos, using alphafold models
"""

import os
import time
import stat
import grp
import sys
import glob
import re
import json
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
from tabulate import tabulate

# from time import strftime
from typing import Dict
from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser
from Bio.PDB import NeighborSearch
from Bio.PDB import Selection
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
                help="Uniprot ID for fasta sequence load")

cmdline_parser.add_argument(
                "--transcript_variant", type=str,
                help="Transcript offset of variant of unknown significance")

cmdline_parser.add_argument("--collate_only", required=False, action='store_true', default=False,
                            help="Skip the calling of PTMGP2, and jump directly to processing output files")

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

trans_mut_pos = int(args.transcript_variant[1:-1])

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
    '--nv',
    '--bind',
    '/home/resv146/',
    '--bind',
    '/data/p_csb_meiler/resv146/PTMGP2/',
    '/home/resv146/vustruct/external_apps/PTMGP2/PTMGP2.simg', # The singularity container
    '/home/resv146/vustruct/external_apps/PTMGP2/GPT2-Inference.py', # The code that will be run "inside" the container
    '--fasta',
    fasta_fullpath,
    '--model_root',
    '/data/p_csb_meiler/resv146/PTMGP2/', # Location of the data directories
    '--output_root',
    args.outdir
    ]


if args.collate_only:
    LOGGER.info("Skipping containerized calls to PTMGP2, assuming already run")
else: # We're running musite deep as usual
    LOGGER.info("Launching container with %s " % ' '.join(singularity_command_list))

    completed_process = subprocess.run(
        singularity_command_list,
        shell=False,
        text=True, 
        # env=environment_override, 
        capture_output=True)

    LOGGER.info("stdout = %s" % completed_process.stdout);
    LOGGER.info("stderr = %s" % completed_process.stderr);
    assert completed_process.returncode ==0, "Exiting with returncode %d" % completed_process.returncode
       

# Now none of this works unless we have an alphafold model - so get that now
alphafold = PDBMapAlphaFold(config_dict)
model_seq_start, model_seq_end, alphafold_modelid = alphafold.best_covering_model(
    args.unp,
    len(transcript_uniprot.aa_seq),
    trans_mut_pos)

alphafold_struct_filename = alphafold.get_coord_filename(alphafold_modelid)
LOGGER.info("Using Alphafold model: %s  seq_start=%d", alphafold_struct_filename, model_seq_start)

if alphafold_struct_filename.endswith(".gz"):
    alphafold_structure_fin = gzip.open(alphafold_struct_filename, 'rt')
elif alphafold_struct_filename.endswith(".cif"):
    alphafold_structure_fin = open(alphafold_struct_filename, 'rt')
else:
    sys.exit("Alphafold model has unrecognized file extension %s" % alphafold_struct_filename)

alphafold_structure_buffer = alphafold_structure_fin.read()
alphafold_structure_fin.close()

alphafold_struct_fin = StringIO(alphafold_structure_buffer)
alphafold_MMCIFParser = MMCIFParser()
alphafold_structure = alphafold_MMCIFParser.get_structure(alphafold_modelid, alphafold_struct_fin)
alphafold_mmCIF_dict = alphafold_MMCIFParser._mmcif_dict
alphafold_struct_fin.close()
alphafold_struct_fin = None

# We need to make sure we can read all the created JSON files from 
# the singularity container run
# the filenames will have forms like:
# output/sp_O-linked Glycosylation (S,T).json
# output/sp_Ubiquitination (K).json
# and the contents are NEGATIVE POSITIVE strings for all the positions contianing the 
# amino acid of interest
ptm_filename_regex_string = os.path.join(args.outdir,"sp_*.json")
json_ptm_prediction_files = glob.glob(ptm_filename_regex_string)

# We extract the name/type of the PTM from this filename using a regex:

ptm_regex_compiled = re.compile(ptm_filename_regex_string.replace("*",r"(.*) \(.*\)") + '$')

# We build out a predictions dictionry of form ptm_predictions_dict[res_no][ptm_type]: 'POSITIVE'
# This paralallels the same dict structure in run_MuSite - except there is not a numerical
# probability for PTMGP2 - just POSITIVE vs NEGATIVE strings
ptm_predictions_dict = {}
for json_ptm_prediction_filename in json_ptm_prediction_files:
    ptm_regex_match = ptm_regex_compiled.match(json_ptm_prediction_filename)
    if ptm_regex_match is None:
        LOGGER.error("Bad json filename: %s", json_ptm_prediction_filename)
        continue
    ptm_type = ptm_regex_match.group(1)
    # print("SUCCESS: %s" % ptm_type)
    json_ptm_dict = {}
    with open(json_ptm_prediction_filename,'r') as json_f:
        try:
            json_ptm_dict = json.load(json_f)
        except json.JSONDecodeError as e:
            json_ptm_dict = {}
            LOGGER.error(f"Invalid json from file {json_ptm_prediction_filename} Error:\n{e}")

    if 'Results' in json_ptm_dict:
        for residue_prediction in json_ptm_dict['Results']:
            # Each residue_prediction is weirdly a residue number key
            # and then either POSITIVE or NEGATIVE
            for res_no_string, ptm_prediction in residue_prediction.items():
                if ptm_prediction == "POSITIVE": # We don't have a numerical value with this PTMGP2
                    res_no = int(res_no_string)
                    if not res_no in ptm_predictions_dict:
                        ptm_predictions_dict[res_no] = {}
                    ptm_predictions_dict[res_no][ptm_type] = 'POSITIVE'
    
LOGGER.debug("ptm_predctions_dict %s" , ptm_predictions_dict)
af_center_residue_id = (' ',1+trans_mut_pos-model_seq_start,' ')
af_center_residue=alphafold_structure[0]['A'][af_center_residue_id]

af_center_atoms=Selection.unfold_entities(af_center_residue, 'A')
af_center_heavy_atoms =[atom for atom in af_center_atoms if atom.name[0] != 'H']

# Create a list of PTM-predicted resdidues that are found in the alphafold
# structure.  Then, winnow to those within 8A of the VUS position...
af_residues  = Selection.unfold_entities(alphafold_structure[0], 'R')

af_ptm_predicted_residues = [res for res in af_residues if (res.id[1]+model_seq_start-1) in ptm_predictions_dict]
af_ptm_predicted_atoms = []
for res in af_ptm_predicted_residues:
    af_ptm_predicted_atoms += [atom for atom in res if atom.name[0] != 'H']

ns=NeighborSearch(af_ptm_predicted_atoms)

# We are building up a set of all alphafold residues here
# That are within 8 angstroms of the central one
af_residues_near_af_center_residue = {res for af_center_heavy_atom in af_center_heavy_atoms
                                for res in ns.search(af_center_heavy_atom.coord, 8.0, 'R')}

LOGGER.info("PTM-predicted residues within 8 angstroms of %d: %s" , trans_mut_pos, af_residues_near_af_center_residue)

af_sorted_near_ptm_residues = sorted(af_res.id for af_res in af_residues_near_af_center_residue)

ptm_neighborhood_df = pd.DataFrame(columns=['unp','Native Residue','PTM Residue','PTM Type', 'PTM Probability'])

for af_near_res_id in af_sorted_near_ptm_residues:
    ptm_dict_for_near_residue = ptm_predictions_dict[af_near_res_id[1]+model_seq_start-1]
    for ptm_type in ptm_dict_for_near_residue:
        ptm_neighborhood_row = {
            'unp': [args.unp],
            'Native Residue': [trans_mut_pos],
            'PTM Residue': [model_seq_start+af_near_res_id[1]-1],
            'PTM Type': [ptm_type],
            'PTM Probability': [ptm_dict_for_near_residue[ptm_type]] # A list of 'POSITIVE' strings
            }
        ptm_neighborhood_df = pd.concat([ptm_neighborhood_df, pd.DataFrame(ptm_neighborhood_row)], ignore_index=True)
    LOGGER.info("Nearby res AF%s: %s" % (af_near_res_id, ptm_predictions_dict[af_near_res_id[1]+model_seq_start-1]))

LOGGER.info("Compare to ALL AF predictions")
for all_res in af_ptm_predicted_residues:
    LOGGER.info("transcript %s = AF %s: %s" % (all_res.id, model_seq_start-1 + all_res.id[1], ptm_predictions_dict[model_seq_start-1 + all_res.id[1]]))


ptm_neighborhood_filename = os.path.join(args.outdir,"PTM_neighborhood.csv")
LOGGER.info("Writing %d PTM predictions within 8A to %s", len(ptm_neighborhood_df), ptm_neighborhood_filename)
ptm_neighborhood_df.to_csv(ptm_neighborhood_filename,sep='\t')

    
# Resume by figuring out that we should keep everything or whatever
# this case of a SHEEP protein...

LOGGER.info("Successful end of run_PTMGP2.py")
psb_status_manager.mark_complete()
sys.exit(0)
