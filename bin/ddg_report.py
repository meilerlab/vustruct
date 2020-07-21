#!/usr/bin/env python
#
# Project        : ddg_report
# Filename       : ddg_report.py
# Authors        : Chris Moth, based on udn_pipeline2.py
# Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : chris.moth@vanderbilt.edu
# Date           : 2020-July
# =============================================================================#

"""ddg_report

   Retrieve ddg calculation results created in the repositry by ddg_run

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
from psb_shared.ddg_repo import DDG_repo
from psb_shared import ddg_clean
from psb_shared.ddg_monomer import DDG_monomer
from psb_shared.psb_progress import PsbStatusManager

from lib import PDBMapSwiss
from lib import PDBMapModbase2016
from lib import PDBMapModbase2013

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


if not args.ddg_config:  # Usually ddg_config will come from a config file - but only require it if not on command line
    required_config_items = [
        # "pdb_dir",
        # "modbase2013_dir",
        # "modbase2013_summary",
        # "modbase2016_dir",
        # "modbase2016_summary",
        # "swiss_dir",
        # "swiss_summary"
        ]

    if not args.ddg_config:  # Usually ddg_config will come from a config file - but only require it if not on command line
        required_config_items.append("ddg_config")

    config, config_dict = psb_config.read_config_files(args, required_config_items)

    config_dict_reduced = {x: config_dict[x] for x in required_config_items}
  
    config_dict = config_dict_reduced

    args.ddg_config = config_dict['ddg_config']

ddg_repo = DDG_repo(args.ddg_config,
                    calculation_flavor='ddg_monomer')
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

ddg_all_results = None
variant_list = None
if args.variant[-1] == '*':
    variant_list = [args.variant[0:-1] + amino_acid for amino_acid in "ACDEFGHIKLMNPQRSTVWY" if amino_acid != args.variant[0] ]
else:
    variant_list = [args.variant]

for variant in variant_list:
    ddg_repo.set_variant(variant)
    ddg_monomer = DDG_monomer(ddg_repo, variant)

    # Get None or a dataframe with results
    ddg_result = ddg_monomer.retrieve_result()

    if ddg_all_results is None:
        ddg_all_results = ddg_result
    else:
       ddg_all_results = pd.concat([ddg_all_results,ddg_result],ignore_index=True)

    print("DDG Retrieved for %s chain=%s %s:\n%s"%(
        ddg_repo.structure_id,
        ddg_repo.chain_id,
        ddg_repo.variant,
        tabulate(ddg_all_results,headers='keys',tablefmt='psql') if ddg_all_results is not None else "Not Found in Repository"))

