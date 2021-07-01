#!/usr/bin/env python
#
# Project        : ddg_run
# Filename       : ddg_run.py
# Authors        : Chris Moth, based on udn_pipeline2.py
# Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : chris.moth@vanderbilt.edu
# Date           : 2020-May
# =============================================================================#

"""ddg_run
   - organizes command line execution of ddg_monomer, ddg_complex
   - Strips (preps) mmCIF or PDB input structures for Rosetta format, leaving
     behind a translation dictionary
   - maintains a repository (repo) of past ddg results

   References: Prediction of protein mutational free energy: benchmark and
               sampling improvements increase classification accuracy
               Brandon Frenz, Steven Lewis, Indigo King, Hahnbeom Park, Frank DiMaio, Yifan Song

   """

import os
import time
import stat
import grp
import sys
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
from configparser import SafeConfigParser
from typing import Dict
from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser
from Bio import BiopythonWarning
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from psb_shared import psb_config
from psb_shared.ddg_load_structure import ddg_load_structure,LoadStructureError
from psb_shared.ddg_repo import DDG_repo
from psb_shared import ddg_clean
from psb_shared.ddg_monomer import DDG_monomer
from psb_shared.psb_progress import PsbStatusManager

from lib import PDBMapSwiss
from lib import PDBMapModbase2020
# from lib import PDBMapModbase2016
# from lib import PDBMapModbase2013

RESOLUTION_QUALITY_MAX=2.5  # Only structures with resolution < 2.5 are routed to the more rigid rosetta "high quality" algorithm

#=============================================================================#
## Function Definitions ##
capra_group = grp.getgrnam('capra_lab').gr_gid



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

cmdline_parser.add_argument(
    "--ddg_config",
    help="ddG Configuration File specifying binary programs, parameters, and rep location",
    required=False, metavar="FILE")
# Input parameters
group = cmdline_parser.add_mutually_exclusive_group(required=True)
group.add_argument('--pdb', type=str,
                   help="4 character PDB ID with optional .chain suffix")
group.add_argument('--biounit', type=str,
                   help="4 character PDB ID with optional .chain suffix")
group.add_argument('--modbase', type=str,
                   help="Modbase 13 or Modbase 16 model ID with optional .chain suffix")
group.add_argument('--swiss', type=str,
                   help="Swissmodel ID with optional .chain suffix")
group.add_argument('--usermodel', type=str, metavar='FILE',
                   help="Filename of a user model.  Requires explicit transcript specifications")
# cmdline_parser.add_argument("entity",type=str,
#                   help="Gene ID, UniProt AC, PDB ID, or PDB filename")

cmdline_parser.add_argument("--chain", type=str,
                            help="Limit the ddG processing to one particular chain")

cmdline_parser.add_argument(
            "-o", "--UDNoutdir", type=str,
            help="Optional directory to echo output and results back to UDN pipeline")
cmdline_parser.add_argument(
            "--UDNuniquekey", type=str, required=False,
            help="Optional gene/refseq/mutation/structure/chain/flavor unique identifer for this pipeline ddG run")
cmdline_parser.add_argument(
            "--variant", type=str,
            help="Amino Acid Variant in modified HGVS format.  PDB insertion codes can be added to residue #")
cmdline_parser.add_argument(
            "--species_filter", type=str, required=False, default=None,
            help="Optionally retain only residues with DBREF of, example 'HUMAN'")

cmdline_parser.add_argument("--label", type=str, default='',
                            help="Optional analysis label (overrides entity inference)")
# cmdline_parser.add_argument("--no-timestamp", "-nt", action="store_true", default=False,
#                            help="Disables output directory timestamping")
# cmdline_parser.add_argument("--overwrite", action="store_true", default=False,
#                             help="Overwrite previous results. Otherwise, exit with error")

LOGGER.info("Command: %s" % ' '.join(sys.argv))

args = cmdline_parser.parse_args()

if args.debug:
    LOGGER.setLevel(logging.DEBUG)
    ch.setLevel(logging.DEBUG)

if not args.UDNoutdir:
    if args.UDNuniquekey:
        args.UDNoutdir = args.uniquekey
    else:
        LOGGER.info("ddg_run.py proceeding independently of UDN")

required_config_items = [
    "pdb_dir",
    "modbase2020_dir",
    "modbase2020_summary",
    # "modbase2013_dir",
    # "modbase2013_summary",
    # "modbase2016_dir",
    # "modbase2016_summary",
    "swiss_dir",
    "swiss_summary"]

if not args.ddg_config:  # Usually ddg_config will come from a config file - but only require it if not on command line
    required_config_items.append("ddg_config")

config, config_dict = psb_config.read_config_files(args, required_config_items)

config_dict_reduced = {x: config_dict[x] for x in required_config_items}

config_dict = config_dict_reduced

if not args.ddg_config:  # Usually ddg_config will come from a config file - but only require it if not on command line
    args.ddg_config = config_dict['ddg_config']

ddg_repo = DDG_repo(args.ddg_config,
                    calculation_flavor='ddG_monomer')
# ddg_repo_dir must be extended below
if args.pdb:
    ddg_structure_dir = ddg_repo.set_pdb(args.pdb.lower(),args.chain)
elif args.swiss:
    ddg_structure_dir = ddg_repo.set_swiss(args.swiss,args.chain)
elif args.modbase:
    ddg_structure_dir = ddg_repo.set_modbase(args.modbase,args.chain)
elif args.usermodel:
    ddg_structure_dir = ddg_repo.set_usermodel(args.usermodel,args.chain)
else:
    message="One of --pdb or --swiss or --modbase or --usermodel required on command line"
    LOGGER.critical(message)
    sys.exit(message)

ddg_repo.set_variant(args.variant)

ddg_repo.make_variant_directory_heirarchy()

# Add a rotating file handler log in the calculation directory
need_roll = os.path.isfile(ddg_repo.log_filename)

logger_fh = RotatingFileHandler(ddg_repo.log_filename, backupCount=5)

logger_fh.setFormatter(log_formatter)
logger_fh.setLevel(logging.DEBUG if args.debug else logging.INFO)
LOGGER.addHandler(logger_fh)

if need_roll:
    logger_fh.doRollover()

LOGGER.info("Job status directory: %s" % ddg_repo.psb_status_manager.status_dir)

def test_completed_earlier():
    complete_timestamp = ddg_repo.psb_status_manager.complete_file_present
    if complete_timestamp:
        local_mtime = time.ctime(complete_timestamp)
        errstr = 'This ddg was previously marked complete at {0}'.format(local_mtime)
        LOGGER.info(errstr)
        print(errstr, file=sys.stderr)
        sys.exit(0)

test_completed_earlier()



ddg_repo.psb_status_manager.clear_status_dir()
ddg_repo.psb_status_manager.write_info('Begun')

# If the repository has the cleaned structure and the xref files, then just use them.

residue_to_clean_xref = ddg_repo.residue_to_clean_xref

structure_info_dict = {}
if residue_to_clean_xref:
    structure_info_dict = ddg_repo.structure_config
else:
    try:
        structure_id,structure,structure_info_dict,mmcif_dict = ddg_load_structure(args,config_dict)
    except LoadStructureError as exception:
        ddg_repo.psb_status_manager.sys_exit_failure(exception.message)


    ddg_cleaner = ddg_clean.DDG_Cleaner(structure,mmcif_dict,True,species_filter=args.species_filter)
    cleaned_structure, residue_to_clean_xref = ddg_cleaner.clean_structure_for_ddg([args.chain])


    from Bio.PDB import PDBIO
    pdbio = PDBIO()
    cleaned_structure_temp_filename = None
    with tempfile.NamedTemporaryFile(delete=False,dir=ddg_structure_dir,mode="w") as cleaned_structure_temp:
        pdbio.set_structure(cleaned_structure)
        cleaned_structure_temp_filename = cleaned_structure_temp.name
        pdbio.save(cleaned_structure_temp_filename , write_end=True, preserve_atom_numbering=False)

    # Save the cleaned structure in the repository
    ddg_repo.mv_cleaned_structure_in(cleaned_structure_temp_filename)

    # Add our cross reference to the repository
    ddg_repo.residue_to_clean_xref = residue_to_clean_xref
    ddg_repo.structure_config = structure_info_dict

# Resume by figuring out that we should keep everything or whatever
# this case of a SHEEP protein...
cleaned_structure_filename = ddg_repo.cleaned_structure_filename
ddg_monomer = DDG_monomer(ddg_repo, mutations=args.variant)
ddg_outcome,ddg_results_df = ddg_monomer.run(
    high_resolution=('resolution' in structure_info_dict) and (float(structure_info_dict['resolution']) < 2.5))

if not ddg_outcome: # Really we should never take this branch
    LOGGER.info("ddg_monomer -FAILED- with message: %s" % ddg_results_df)
    sys.exit(ddg_results_df)

LOGGER.info("Successful end of ddg_run:\n%s"%tabulate(ddg_results_df,headers='keys',tablefmt='psql'))

sys.exit(0)

import pdb; pdb.set_trace()



"""
# if args.biounit:
#  is_biounit,
# chain_to_transcript)
if args.usermodel:
    if not args.label:
        exitmsg = "--label is required on the command line when loading a --usermodel"
        LOGGER.critical(exitmsg)
        ddg_repo.psb_status_manager.sys_exit_failure(exitmsg)
    structure = load_structure(args.label, args.usermodel)
elif args.swiss:  # This is a lenghty-ish swiss model ID, NOT the file location
    PDBMapSwiss.load_swiss_INDEX_JSON(config_dict['swiss_dir'], config_dict['swiss_summary']);
    structure = load_structure(args.swiss, PDBMapSwiss.get_coord_file(args.swiss))

    assigned_chain_id = next(iter(chain_to_transcript))

    # Swiss _could_ be a homo-oliomer - and we want to capture that!
    # by pointing the additional chain IDs at the same alignment
    for chain in structure[0]:
        if chain.id not in chain_to_transcript:
            first_residue = next(iter(structure[0][chain.id]))
            # Swiss model seems to put some HETATMs in some of the chains.  DO NOT align to those!
            if first_residue.id[0] == ' ':
                chain_to_transcript[chain.id] = chain_to_transcript[assigned_chain_id]


    # Finally - if the chain left has no id, set the id to A for sanity
for chain in list(structure.get_chains()):
    if chain.id in ('', ' '):
        LOGGER.info('Renaming blank/missing chain ID to A')
        chain.id = 'A'



assert structure, statusdir_info(
    "A structure file must be specified via --pdb, --biounit, --swiss, --modbase, or --usermodel")

print("AWESOME - %s"%structure)

# config.items() returns a list of (name, value) pairs - convert to dict
config_dict = dict(config.items("Genome_PDB_Mapper"))





print('Test me')
sys.exit(1)
"""
