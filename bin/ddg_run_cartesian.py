#!/usr/bin/env python
#
# Project        : ddg_run_cartesian
# Filename       : ddg_run_cartesian.py
# Authors        : Chris Moth, based on ddg_cartesian.py
# Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : chris.moth@vanderbilt.edu
# Date           : 2020-May
# =============================================================================#

"""ddg_relax
   - Strips (preps) mmCIF or PDB input structures for Rosetta format,
   - relax structure in preparation for ddg_cartesian.
   - calls ddg_Cartesian on the lowest energy relax output structure.
   - maintains a repository (repo) of past ddg results

   References:
   J. Chem. Theory Comput. 2016, 12, 6201−6212
   Park et al.
     Simultaneous Optimization of Biomolecular Energy Functions on
      Features from Small Molecules and Macromolecules

   Frontiers in Bioengineering and Biotechnology   October 2020 | Volume 8 | Article 558247
   Frenz et al.
     Prediction of Protein Mutational Free Energy: Benchmark and Sampling
     Improvements Increase Classification Accuracy

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
from typing import Dict
from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser
from Bio import BiopythonWarning
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from psb_shared import psb_config
from psb_shared.ddg_load_structure import ddg_load_structure,LoadStructureError
from psb_shared.ddg_repo import DDG_repo,DDG_repo_RotatingFileHandler
from psb_shared import ddg_clean
from psb_shared.ddg_commandline_parser import ddg_commandline_parser
from psb_shared.ddg_cartesian import DDG_cartesian
from psb_shared.psb_progress import PsbStatusManager

from lib import PDBMapSwiss
from lib import PDBMapModbase2020
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

ddg_commandline_parser(cmdline_parser)

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
    "alphafold_dir",
    "modbase2020_dir",
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
                    calculation_flavor='ddG_cartesian')
# ddg_repo_dir must be extended below
if args.pdb:
    ddg_structure_dir = ddg_repo.set_pdb(args.pdb.lower(),args.chain)
elif args.swiss:
    ddg_structure_dir = ddg_repo.set_swiss(args.swiss,args.chain)
elif args.modbase:
    ddg_structure_dir = ddg_repo.set_modbase(args.modbase,args.chain)
elif args.alphafold:
    ddg_structure_dir = ddg_repo.set_alphafold(args.alphafold,args.chain)
elif args.usermodel:
    ddg_structure_dir = ddg_repo.set_usermodel(args.usermodel,args.chain)
else:
    message="One of --pdb, --swiss, --modbase, --alphafold, or --usermodel required on command line"
    LOGGER.critical(message)
    sys.exit(message)

ddg_repo.set_variant(args.variant)

ddg_repo.make_variant_directory_heirarchy()


# Add a rotating file handler log in the calculation directory
need_roll = os.path.isfile(ddg_repo.log_filename)




def _open(self):
    """
    Open the current base file with the (original) mode and encoding.
    Return the resulting stream.
    """
    open_func = self._builtin_open
    return open_func(self.baseFilename, self.mode,
                     encoding=self.encoding, errors=self.errors)

logger_fh = DDG_repo_RotatingFileHandler(ddg_repo, backupCount=5)

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
# Unlike DDG Monomer, residues are not renumbered for DDG_Cartesian. However this file
# is still used to answer whether a structure has been cleaned, for purposes
residue_to_clean_xref = ddg_repo.residue_to_clean_xref

structure_info_dict = {}
if residue_to_clean_xref:
    structure_info_dict = ddg_repo.structure_config
else:
    try:
        structure_id, structure,structure_info_dict, mmcif_dict = ddg_load_structure(args,config_dict)
    except LoadStructureError as exception:
        ddg_repo.psb_status_manager.sys_exit_failure(exception.message)

    ddg_cleaner = ddg_clean.DDG_Cleaner(structure,
                                        mmcif_dict,
                                        True,
                                        species_filter=args.species_filter,
                                        renumber_residues=True)
    cleaned_structure, residue_to_clean_xref = ddg_cleaner.clean_structure_for_ddg([args.chain])


    from Bio.PDB import PDBIO
    pdbio = PDBIO()
    cleaned_structure_temp_filename = None
    with tempfile.NamedTemporaryFile(delete=False,dir=ddg_structure_dir,mode="w") as cleaned_structure_temp:
        cleaned_structure_temp_filename = cleaned_structure_temp.name
        ddg_repo.os_file_chmod(cleaned_structure_temp_filename)
        ddg_repo.set_group(cleaned_structure_temp_filename)
        pdbio.set_structure(cleaned_structure)
        # with os.fdopen(ddg_repo.os_open(cleaned_structure_temp_filename,'w'),'w') as cleaned_structure_f:
        pdbio.save(cleaned_structure_temp, write_end=True, preserve_atom_numbering=False)

    # Save the cleaned structure in the repository
    ddg_repo.mv_cleaned_structure_in(cleaned_structure_temp_filename)

    # Add our cross reference to the repository
    ddg_repo.residue_to_clean_xref = residue_to_clean_xref
    ddg_repo.structure_config = structure_info_dict

# Resume by figuring out that we should keep everything or whatever
# this case of a SHEEP protein...
cleaned_structure_filename = ddg_repo.cleaned_structure_filename
LOGGER.info("Cleaned structure is %s", ddg_repo.cleaned_structure_filename)
ddg_cartesian = DDG_cartesian(ddg_repo, mutations=args.variant)
ddg_outcome,ddg_results_df = ddg_cartesian.run()

if not ddg_outcome:  # Really we should never take this branch
    LOGGER.info("ddg_cartesian -FAILED- with message: %s" % ddg_results_df)
    sys.exit(ddg_results_df)

LOGGER.info("Successful end of ddg_run_cartesian:\n%s"%tabulate(ddg_results_df,headers='keys',tablefmt='psql'))
sys.exit(0)
