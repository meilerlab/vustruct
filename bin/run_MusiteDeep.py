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

The prediction algoritm is strict sequence based.
After the predictions are made, we report all those residues with heavy atoms 
with 8 Angstroms of the mutation pos, using alphafold models
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
                            help="Skip the calling of MusiteDeep, and jump directly to processing output files")

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
    '--bind',
    '/dors/capra_lab',
    '/dors/capra_lab/users/mothcw/VUStruct/musitedeep.simg',
    '/bin/bash']


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

if args.collate_only:
    LOGGER.info("Skipping containerized calls to MusiteDeep, assuming already run")
else: # We're running musite deep as usual
    LOGGER.info("Launching container with %s " % ' '.join(singularity_command_list))
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

############################
# With all the processing run, now gather the outputs
# Step one is to gather the outputs into a dictionary indexed on residue number
# I iterate over all the _possible_ output files.
# The format of the several output files from the MusiteDeep is very weird: 
# And I parse out the "hits" by emulating Alican's check_best.sh code
"""Position	Residue	PTMscores	Cutoff=0.5
>sp|P51688|SPHM_HUMAN N-sulphoglucosamine sulphohydrolase OS=Homo sapiens OX=9606 GN=SGSH PE=1 SV=1
103	K	Ubiquitination:0.489	None
123	K	Ubiquitination:0.292	None
124	K	Ubiquitination:0.595	Ubiquitination:0.595
156	K	Ubiquitination:0.201	None
161	K	Ubiquitination:0.207	None
196	K	Ubiquitination:0.381	None
303	K	Ubiquitination:0.371	None
339	K	Ubiquitination:0.327	None
393	K	Ubiquitination:0.187	None
425	K	Ubiquitination:0.298	None
470	K	Ubiquitination:0.317	None
490	K	Ubiquitination:0.629	Ubiquitination:0.629
"""

ptm_predictions_dict: dict[int, dict[str, float]] = {}
for musite_model_prefix in musite_model_prefix_list:
    predictions_filename = os.path.join(args.outdir, args.unp + '_' + musite_model_prefix + "_results.txt")
    with open(predictions_filename, 'r') as predictions_f:
        header1 = predictions_f.readline()
        # Do a sanity check on the first line
        if header1.find("Position") == 0:
            header2 = predictions_f.readline()
            if len(header2) > 0 and header2[0] == ">":
                for line in predictions_f:
                    cols = line.strip().split('\t')
                    if len(cols) == 4:
                        # A None in the last column means that we did not meet the 0.5 threshold
                        # given to the prediction program command line
                        if "None" not in cols[3]:
                            res_no = int(cols[0])
                            ptm_and_probability = cols[3].split(':')
                            if len(ptm_and_probability) == 2:
                                ptm_type = ptm_and_probability[0]
                                probability = float(ptm_and_probability[1])
                                # print ("%s %s %s" % (res_no, ptm, probability))
                                if not res_no in ptm_predictions_dict:
                                    ptm_predictions_dict[res_no] = {}
                                ptm_predictions_dict[res_no][ptm_type] = probability

print("%s" % ptm_predictions_dict)

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




LOGGER.info("ptm_predctions_dict %s" , ptm_predictions_dict)
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
            'PTM Probability': [ptm_dict_for_near_residue[ptm_type]]
            }
        ptm_neighborhood_df = pd.concat([ptm_neighborhood_df, pd.DataFrame(ptm_neighborhood_row)], ignore_index=True)
    LOGGER.info("Nearby res AF%s: %s" % (af_near_res_id, ptm_predictions_dict[af_near_res_id[1]+model_seq_start-1]))

print("Compare to ALL AF predictions")
for all_res in af_ptm_predicted_residues:
    LOGGER.info("transcript %s = AF %s: %s" % (all_res.id, model_seq_start-1 + all_res.id[1], ptm_predictions_dict[model_seq_start-1 + all_res.id[1]]))


ptm_neighborhood_filename = os.path.join(args.outdir,"PTM_neighborhood.csv")
LOGGER.info("Writing %d PTM predictions within 8A to %s", len(ptm_neighborhood_df), ptm_neighborhood_filename)
ptm_neighborhood_df.to_csv(ptm_neighborhood_filename,sep='\t')

    


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
