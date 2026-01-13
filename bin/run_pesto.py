#!/usr/bin/env python3
#
# Project        : VUSTRUCT: PeSTo integration for Protein-Prot/Ligand/DNA/etc
#                interaction :           prediction
# Filename       : run_pesto.py
#                :
# Author         : Chris Moth
# Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : chris.moth@vanderbilt.edu
# Date           : 2025 July
# Description    : Renumbers an input PDB or model to match transcript number
#                : Then calls pesto to calculate PPIs of different forms
#                : spatial position relative to known pathogenic and neutral
#                : variation in protein structure.
# =============================================================================#
## Package Dependenecies ##
import sys
import os
import shutil
import logging
import gzip
import lzma
import csv
import platform
import grp
import stat
import json
import re
import string
import subprocess
from time import strftime
from itertools import combinations
from collections import OrderedDict
import pprint
from subprocess import Popen, PIPE

from enum import Enum

class ResidueInteractionType(Enum):
    PROTEIN = 0
    DNA_RNA = 1
    ION     = 2
    LIGAND  = 3
    LIPID   = 4

from typing import Dict



import configparser

# Numerical
import pandas as pd
import numpy as np
import math

TOL = 1e-5  # zero-tolerancection with a given name return the same logger instance. This means that logger instances never need to be passed between d

# Stats
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.stats import norm, percentileofscore, mannwhitneyu
from scipy.spatial import KDTree

import Bio

from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.SeqUtils import seq3
from Bio.SeqUtils import seq1
# Import our PDBMap library
from lib import PDBMapGlobals
from lib import PDBMapAlignment
from lib import PDBMapProtein
from lib import PDBMapSQLdb
from lib import PDBMapTranscriptBase
from lib import PDBMapTranscriptUniprot
from lib import PDBMapTranscriptFasta
from lib import PDBMapTranscriptEnsembl
from lib import PDBMapGenomeVariants
from lib import PDBMapGnomad
from lib import PDB36Parser
from lib import PDBMapComplex


if __name__ == "__main__":
    # Logging setup can vary a lot depending on whether pathprox is used as import module or mainline route
    ch = logging.StreamHandler()
    LOGGER = logging.getLogger()
    LOGGER.addHandler(ch)

    log_format_string = '%(asctime)s %(levelname)-4s [%(filename)16s:%(lineno)d] %(message)s'
    date_format_string = '%H:%M:%S'
    log_formatter = logging.Formatter(log_format_string, date_format_string)
    ch.setFormatter(log_formatter)

    LOGGER.setLevel(logging.DEBUG)
    ch.setLevel(logging.INFO)

    from psb_shared import psb_config, psb_progress, psb_capra_group
    from psb_shared.psb_progress import PsbStatusManager

    # Setup the Command Line Argument Parser
    cmdline_parser = psb_config.create_default_argument_parser(__doc__, os.path.dirname(os.path.dirname(__file__)),
                                                               add_help=True)

    cmdline_parser.add_argument("-o", "--outdir", type=str, required=True,
                                help="Directory to use for output and results")
    cmdline_parser.add_argument("--uniquekey", type=str, required=False,
                                help="A gene/refseq/mutation/structure/chain/flavor unique identifer for this PeSTO job")
    cmdline_parser.add_argument("--transcript_variant", type=str, required=True,
                                help="Variant of interest in form AnnnB, for reporting in final json file")
    # Input parameters
    group = cmdline_parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--pdb', type=str,
                       help="4 character PDB ID with optional .chain suffix")
    group.add_argument('--biounit', type=str,
                       help="4 character PDB ID with optional .chain suffix")
    group.add_argument('--modbase', type=str,
                       help="Modbase 20 model ID with optional .chain suffix")
    group.add_argument('--swiss', type=str,
                       help="Swissmodel ID with optional .chain suffix")
    group.add_argument('--alphafold', type=str, metavar='FILE',
                       help="ID of an alphafold model with optional .chain suffix")
    group.add_argument('--usermodel', type=str, metavar='FILE',
                       help="Filename of a user model.  Requires explicit transcript specifications")
    # cmdline_parser.add_argument("entity",type=str,
    #                   help="Gene ID, UniProt AC, PDB ID, or PDB filename")
    args, args_remaining = cmdline_parser.parse_known_args()

    if args.debug:
        LOGGER.setLevel(logging.DEBUG)


    psb_capra_group.makedirs_capra_lab(args.outdir, 'Main outdir creation')

    # cache_dir = "/tmp/cache"
    log_filename = os.path.join(args.outdir, args.uniquekey + ".log")

    LOGGER.info("Log file is %s" % log_filename)
    needRoll = os.path.isfile(log_filename)

    # Write current version of this script to output directory
    # pw_name = pwd.getpwuid( os.getuid() ).pw_name # example jsheehaj or mothcw
    from logging.handlers import RotatingFileHandler
    from logging import handlers

    fh = RotatingFileHandler(log_filename, maxBytes=(1048576 * 5), backupCount=7)
    # formatter = logging.Formatter('%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s')
    fh.setFormatter(log_formatter)
    fh.setLevel(logging.INFO)
    LOGGER.addHandler(fh)

    if needRoll:
        fh.doRollover()
else:
    LOGGER = logging.getLogger(__name__)

import inspect  # Update status with stack trace info if we crash

# Warnings
import warnings

if __name__ == "__main__":
    curdir = os.getcwd()  # Note the current directory

    # cmdline_parser.add_argument("refseq",type=str,
    #                  help="NM_.... refseq transcript ID")

    # Filter parameters
    cmdline_parser.add_argument("--chain", type=str,
                                help="Limit the analysis to one particular chain")
    cmdline_parser.add_argument("--unp", type=str,
                                help="Explicit declaration of the reference uniprot id")

    LOGGER.info("Command: %s" % ' '.join(sys.argv))

    # cmdline_parser.add_argument("collaboration",type=str,help="Collaboration ID (ex. UDN)")
    cmdline_parser.add_argument("mutation", nargs='?', type=str, default='', help="HGVS mutation string (ex S540A)")

    args, args_remaining = cmdline_parser.parse_known_args()  # args_remaining)
    # if not os.path.exists(cache_dir):
    #  os.makedirs(cache_dir)


    LOGGER.info("Command: %s" % ' '.join(sys.argv))

    psb_status_manager = PsbStatusManager(args.outdir)
    LOGGER.info("Job status directory: %s" % psb_status_manager.status_dir)


    def sys_exit_failure(info):
        psb_status_manager.sys_exit_failure(info)
        sys.exit(info)


    psb_status_manager.clear_status_dir()

    psb_status_manager.write_info('Begun')

    required_config_items = [# "dbhost", "dbname", "dbuser", "dbpass",
                             "collaboration",
                             "pdb_dir",
                             "sec2prim",
                             "collaboration",
                             "idmapping",
                             "interpro_dir",
                             "modbase2020_dir",
                             "modbase2020_summary",
                             "alphafold_dir",
                             "gnomad_dir",
                             "gnomad_filename_template",
                             "output_rootdir",
                             "sprot",
                             "swiss_dir",
                             "swiss_summary"]

    config, config_dict = psb_config.read_config_files(args, required_config_items)
    PDBMapSQLdb.set_access_dictionary(config_dict)

    config_dict_reduced = {x: config_dict[x] for x in required_config_items}

    config_dict = config_dict_reduced

    PDBMapGlobals.config = config_dict
    PDBMapGlobals.exit_function = sys_exit_failure

    config_dict_shroud_password = {x: config_dict[x] for x in required_config_items}
    dbpass = config_dict.get('dbpass', '?')
    config_dict_shroud_password['dbpass'] = '*' * len(dbpass)

    variant_resno = int(args.transcript_variant[1:-1])
    # The biopython variant residue ID in the apply_model.py written output
    # Which will be the renumbered PDB of course
    # We need to keep thinking through this!
    variant_resid = (' ', variant_resno, ' ')

    pesto_predictions_dict = {args.transcript_variant:{}}

    # pathprox_config_dict = dict(config.items('PathProx'))

    LOGGER.info("PATH=%s", os.getenv("PATH"))
    LOGGER.info("PYTHONPATH=%s", os.getenv("PYTHONPATH"))
    LOGGER.info("sys.path=%s", sys.path)

    # LOGGER.info("Command Line Arguments")
    # pprint.pprint(vars(args))
    LOGGER.info("Command Line Arguments:\n%s" % pprint.pformat(vars(args)))
    LOGGER.info("Configuration File parameters:\n%s" % pprint.pformat(config_dict_shroud_password))
    # LOGGER.info("PathProx parameters:\n%s"%pprint.pformat(pathprox_config_dict))

    # Step 1 is to parse out the assignments of chains to ids
    PDBMapProtein.load_idmapping(PDBMapGlobals.config['idmapping'])

    # PROCESS USER-SPECIFIED transcript<->chain mappings
    # if the user has any transcript->chain assignments, parse them from command line

    # transcript_parser = argparse.ArgumentParser()
    alpha_fold_local_metrics = None

    # Parse the --chainBunp=P12435, --chainXenst=ENST00000012432 type command line arguments:
    if args_remaining:
        LOGGER.info("Additional comamnd line arguments for chain assignments: %s" % args_remaining)

    user_chain_to_transcript: dict[str, str] = PDBMapComplex.parse_user_chain_to_transcript(args_remaining)

    if args.pdb:  # Then this is going to be a single chain calculation - so do try a biounit which may lack chain of interest
        complex = PDBMapComplex('pdb', args.pdb, args.chain, user_chain_to_transcript)

    elif args.biounit:
        complex = PDBMapComplex('biounit', args.biounit, args.chain, user_chain_to_transcript)

    elif args.usermodel:
        complex = PDBMapComplex('usermodel',
                                args.usermodel,
                                args.chain,
                                user_chain_to_transcript,
                                "PeSTo")
    elif args.swiss:  # This is a lenghty-ish swiss model ID, NOT the file location
        complex = PDBMapComplex('swiss', args.swiss, args.chain, user_chain_to_transcript=user_chain_to_transcript)

    elif args.modbase:  # This is a ENSP.... modbase file.  Currently we look only in modbase20
        complex = PDBMapComplex('modbase',args.modbase, user_chain_to_transcript=user_chain_to_transcript)

    elif args.alphafold:
        complex = PDBMapComplex('alphafold', args.alphafold, user_chain_to_transcript=user_chain_to_transcript)

    complex.assign_chain_letters_where_blank()

    # residue_coms = {}

    # for chain in model:
    #    for residue in chain:
    #        # Compute the center of mass for all residues
    #        natom = float(len([atom for atom in residue]))
    #        com = sum([atom.coord for atom in residue]) / natom
    #        residue_coms[(chain.id, residue.id)] = com
    renumbered_structure = complex.structure_renumber_per_alignments_and_lose_disorder()

    if args.verbose:
        LOGGER.info("Active options:")
        for arg in vars(args):
            try:
                LOGGER.info("  %15s:  %s" % (arg, getattr(args, arg).name))
            except:
                LOGGER.info("  %15s:  %s" % (arg, getattr(args, arg)))
        LOGGER.info("")

    AAseq = None

    psb_status_manager.write_info('Configured')

    # =============================================================================#
    ## Begin Analysis ##

    # Init PDBMap by load the summary information dictionaries that it will need to use:
    LOGGER.info("Initializing PDBMAP by loading idmapping sec2prom, sprot, modbase, and swiss model meta dictionaries")
    PDBMapProtein.load_sec2prim(config_dict['sec2prim'])
    PDBMapProtein.load_sprot(config_dict['sprot'])
    psb_status_manager.write_info('Initialized')

    #
    # Write the complex out using the renumbered chains in  output directory
    #
    LOGGER.info("Writing renumbered PDB for PeSTo in %s", args.outdir)
    os.chdir(args.outdir)
    complex.write_renumbered(file_format='pdb')
    os.chdir(curdir)

    # Sanity check the written ...renum.pdb before jumping in
    from glob import glob
    input_renum_pdb_filename = glob(os.path.join(args.outdir, "*renum.pdb"), recursive=False)

    # glob returns list - we need one filename
    assert len(input_renum_pdb_filename) == 1
    input_renum_pdb_filename = input_renum_pdb_filename[0]

    LOGGER.info("Verifying correct transcript residue %s AA is found in\npesto input structure:%s", args.transcript_variant, input_renum_pdb_filename)
    pdb_parser=PDBParser()
    pesto_input_pdb = pdb_parser.get_structure("renum_check", input_renum_pdb_filename)

    # We sadly have to fiddle with traversing this way because pesto outputs 1-based residue numbers.:w
    def find_nth_residue() -> int:
        variant_at_nth_residue = 0
        for _model in pesto_input_pdb:
            for _chain in _model:
                assert variant_resid in _chain, "The variant %s is not in chain %s of the structure %s" % (
                    args.transcript_variant, _chain.id, input_renum_pdb_filename)
                for residue in _chain:
                    if residue.id == variant_resid:
                        return variant_at_nth_residue
                    else:
                        variant_at_nth_residue += 1

    variant_at_nth_residue = find_nth_residue()
    LOGGER.info("Transcript residue %s is %d residues into the first chain", args.transcript_variant, variant_at_nth_residue)
    singularity_command_list = [
        'singularity',
        'exec',
        # '--cleanenv',
        '--env',
        'PYTHONPATH=/opt/PeSTo',
        '--bind',
        '/home/resv146',
        '/home/resv146/containers/PeSTo.simg',
        '/home/resv146/vustruct/external_apps/PeSTo/apply_model.py',
        '--pytorch_device=cpu',
        '--data_dir=%s' % args.outdir
        ]

    LOGGER.info("Launching container with %s " % ' '.join(singularity_command_list))

    environment_without_PYTHONNOUSERSITE = os.environ.copy()
    if 'PYTHONNOUSERSITE' in environment_without_PYTHONNOUSERSITE:
        LOGGER.warn("Removing PYTHONNOUSERSITE from subprocess environment")
        del environment_without_PYTHONNOUSERSITE['PYTHONNOUSERSITE']
    

    completed_process = subprocess.run(
        singularity_command_list,
        shell=False,
        text=True,
        env=environment_without_PYTHONNOUSERSITE,
        capture_output=True)
    LOGGER.info("returncode = %d" % completed_process.returncode);
    LOGGER.info("stdout = %s" % completed_process.stdout);
    LOGGER.info("stderr = %s" % completed_process.stderr);
    completion_message = "PeSTo apply_model.py failed with returncode = %d" % completed_process.returncode
    if completed_process.returncode != 0:
        LOGGER.error(completion_message)
        psb_status_manager.write_info(completed_process.stderr)
        psb_status_manager.mark_failed()
        psb_status_manager.sys_exit_failure(completion_message)
    
    LOGGER.info("PeSTo apply_model.py ended with returncode = 0")
    psb_status_manager.write_info("Success")
    psb_status_manager.mark_complete()

    from glob import glob
    pdb_filepaths = glob(os.path.join(args.outdir, "*i[0-4].pdb"), recursive=False)
    assert len(pdb_filepaths) == 5, "Somehow the needed *i[0-4].pdb files are not output"
    LOGGER.info("PeSTo apply_model.py outputs pdbs are: %s.....", pdb_filepaths)

    pdb_filebase = pdb_filepaths[0][:-5] # Get rid of last 5 chars [0-4].pdb

    for residue_interaction in ResidueInteractionType:
        pdb_filename = pdb_filebase + str(residue_interaction.value) + ".pdb"
        LOGGER.info("%-.15s: %s", residue_interaction.name, pdb_filename)
        assert os.path.isfile(pdb_filename)
        pdb_parser=PDBParser()
        pesto_predictions_pdb = pdb_parser.get_structure(residue_interaction.name, pdb_filename)

        first_chain = next(iter(pesto_predictions_pdb[0]))

        residue_counter = 0
        variant_res = None
        for nth_residue in first_chain:
            if residue_counter >= variant_at_nth_residue:
                variant_res = nth_residue
                break
            residue_counter += 1 
            

        LOGGER.info("Original transcript_variant %s is residue %s in the pesto output", args.transcript_variant, nth_residue.id[1])
        # Get the three-letter code
        _three_letter_code = variant_res.resname
        # Convert to one-letter code
        _one_letter_code = seq1(_three_letter_code)
        assert _one_letter_code == args.transcript_variant[0], "transcript_variant expects %s but %s found in structure" % (args.transcript_variant[0], _one_letter_code)
        for atom in variant_res:
            bfactor = atom.get_bfactor()
            pesto_predictions_dict[args.transcript_variant][residue_interaction.name] = float(bfactor)
            break # No need to read every atom - as they'll all have the same bfactor

    pesto_predictions_filename = os.path.join(args.outdir,"pesto_predictions.json")
    with open(pesto_predictions_filename, 'w') as pesto_predictions_f:
        pesto_predictions_f.write(json.dumps(pesto_predictions_dict))

    LOGGER.info("PeSTo predictions for %s written to file %s", args.transcript_variant, pesto_predictions_filename)
        

