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
from psb_shared.ddg_repo import DDG_repo
from psb_shared import ddg_clean
from psb_shared.ddg_monomer import DDG_monomer
from psb_shared.psb_progress import PsbStatusManager

from lib import PDBMapSwiss
from lib import PDBMapModbase2016
from lib import PDBMapModbase2013

RESOLUTION_QUALITY_MAX=2.5  # Only structures with resolution < 2.5 are routed to the more rigid rosetta "high quality" algorithm

#=============================================================================#
## Function Definitions ##
capra_group = grp.getgrnam('capra_lab').gr_gid

def load_structure(id,coord_filename):
    with gzip.open(coord_filename,'rt') if coord_filename.split('.')[-1]=="gz" else open(coord_filename,'r') as fin:
        warnings.filterwarnings('ignore',category=PDBConstructionWarning)
        structure = PDBParser().get_structure(id,fin)
        warnings.resetwarnings()
    return structure

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


# def __info_update(info):
#     new_info_filename = os.path.join(statusdir, "info.new")
#     with open(new_info_filename, 'w') as f:
#         f.write(info + '\n')
#     final_info_filename = os.path.join(statusdir, "info")
#     os.rename(new_info_filename, final_info_filename)
#     LOGGER.info("%s now contains: %s" % (final_info_filename, info))
# 
# 
# def sys_exit_failure(info):
#     __info_update(info)
#     new_progress_filename = os.path.join(statusdir, "progress.new")
#     with open(new_progress_filename, 'w') as f:
#         f.write("%s: %s\n" % (__file__, inspect.currentframe().f_back.f_lineno))
#     os.rename(new_progress_filename, '%s/progress' % statusdir)
#     # Mark this job as failed
#     fail_filename = os.path.join(statusdir, "FAILED")
#     open(fail_filename, 'w').close()
#     LOGGER.critical("Creating FAILURE file %s" % fail_filename)
#     sys.exit(info)
# 
# 
# def statusdir_info(info):
#     __info_update(info)
#     new_progress_filename = os.path.join(statusdir, "progress.new")
#     with open(new_progress_filename, 'w') as f:
#         f.write("%s: %s\n" % (__file__, inspect.currentframe().f_back.f_lineno))
#     os.rename(new_progress_filename, '%s/progress' % statusdir)
#     return info




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
    "modbase2013_dir",
    "modbase2013_summary",
    "modbase2016_dir",
    "modbase2016_summary",
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

ddg_repo.set_variant(args.variant)

ddg_repo.make_calculation_directory_heirarchy()

# Add a rotating file handler log in the calculation directory
need_roll = os.path.isfile(ddg_repo.log_filename)

logger_fh = RotatingFileHandler(ddg_repo.log_filename, backupCount=5)

logger_fh.setFormatter(log_formatter)
logger_fh.setLevel(logging.DEBUG if args.debug else logging.INFO)
LOGGER.addHandler(logger_fh)

if need_roll:
    logger_fh.doRollover()

LOGGER.info("Job status directory: %s" % ddg_repo.psb_status_manager)

def test_completed_earlier():
    complete_timestamp = ddg_repo.psb_status_manager.complete_file_present
    if complete_timestamp:
        local_mtime = time.ctime(complete_timestamp)
        errstr = 'This ddg was previously marked complete at {0}'.format(local_mtime)
        LOGGER.info(errstr)
        print(errstr, file=sys.stderr)
        sys.exit(0)

test_completed_earlier()

def mine_structure_info_from_mmCIF(mmCIF_dict):
    method = mmCIF_dict['_exptl.method'][0]
    pdb_id = method,mmcif_dict['_entry.id']
    deposition_date = None

    raw_deposition_date = mmCIF_dict.get('_pdbx_database_status.recvd_initial_deposition_date',None)
    if raw_deposition_date: # Split into YYYY-MM-DD
        raw_YYYY_MM_DD = raw_deposition_date[0].split('-')
        if len(raw_YYYY_MM_DD) == 3:
            deposition_date = datetime.datetime(
                int(raw_YYYY_MM_DD[0]),int(raw_YYYY_MM_DD[1]),int(raw_YYYY_MM_DD[2]))

        if not deposition_date:
            LOGGER.warning("Deposition date %s apparently not YYYY-MM-DD",raw_deposition_date)

    resolution = None # Not applicable is correct for NMR structures
    if method == 'X-RAY DIFFRACTION':
        resolution_keys = ['_refine.ls_d_res_high','_reflns.d_resolution_high','_refine_hist.d_res_high']
        ethod = 'X-RAY'
    elif method == 'ELECTRON MICROSCOPY':
        method = 'EM'
        resolution_keys = ['_em_3d_reconstruction.resolution']
    elif method.find('NMR') > -1: # resolution is not applicable for nmr
        resolution_keys = []
    else:
        resolution_keys = ['_refine.ls_d_res_high','_reflns.d_resolution_high','_refine_hist.d_res_high']
        LOGGER.critical('%s Experimental method for %s is unknown'%(method,pdb_id))

    for resolution_key in resolution_keys:
        if resolution_key in mmCIF_dict:
            # Congratulations.  You found the old REMARK 2 RESOLUTION entry.... or so we hope.
            resolution = float(mmCIF_dict[resolution_key][0])
            break
    if not resolution and method.find('X-RAY') > 0:
        LOGGER.warning("pdb %s is method=%s with no resolution entry"%(pdb_id,method))

    structure_info_dict =  {'method': method, 'deposition_date': deposition_date}
    if resolution:
        structure_info_dict['resolution'] = resolution
    return structure_info_dict


ddg_repo.psb_status_manager.clear_status_dir()
ddg_repo.psb_status_manager.write_info('Begun')

# If the repository has the cleaned structure and the xref files, then just use them.

residue_to_clean_xref = ddg_repo.residue_to_clean_xref

structure_info_dict = {}
if residue_to_clean_xref:
    structure_info_dict = ddg_repo.structure_config
if not structure_info_dict:
    # We need to load the structure, clean it, and save the xref
    if args.pdb:
        # Load mmCIF from pdb repository
        structure_id = args.pdb.lower()
        original_structure_filename = os.path.join(
            config_dict['pdb_dir'],
            'structures',
            'divided',
            'mmCIF',
            structure_id[1:3],
            "%s.cif.gz" % structure_id)
        original_structure_type = 'mmCIF'
    elif args.swiss:
        PDBMapSwiss.load_swiss_INDEX_JSON(config_dict['swiss_dir'],config_dict['swiss_summary']);
        swiss_info = PDBMapSwiss.get_info(args.swiss)
        assert swiss_info,"%s was not found in the swiss INDEX JSON file"%args.swiss
        original_structure_type = 'pdb'
        original_structure_filename = PDBMapSwiss.get_coord_file(args.swiss)
        structure_id = swiss_info['modelid']
        structure_info_dict['pdb_template'] = swiss_info['template']
        # The swiss template info has source PDB chain in the last position
        pdb_template_chain = structure_info_dict['pdb_template'][-1]
        assert pdb_template_chain == args.chain,\
            "With swiss models, chain on command line %s must match modelled chain %s"%(
                args.chain,pdb_template_chain)
        remark3_metrics = PDBMapSwiss.load_REMARK3_metrics(structure_id)

        if not DDG_monomer.evaluate_swiss(structure_id,remark3_metrics):
            ddg_repo.psb_status_manager.sys_exit_failure("Terminating as %s of insufficient quality for ddg monomer"%args.swiss)

        structure_info_dict['method'] = remark3_metrics['mthd']
        structure_info_dict['template_identity'] = float(remark3_metrics['sid'])
    elif args.modbase:  # This is a ENSP.... modbase file.  Could be either Modbase13 or 16 - so look in both places
        original_structure_type = 'pdb'
        modbase16 = PDBMapModbase2016(config_dict)
        original_structure_filename = modbase16.get_coord_file(args.modbase)
        if os.path.exists(original_structure_filename):
            structure = load_structure(args.modbase, original_structure_filename)
        else:
            modbase13 = PDBMapModbase2013(config_dict)
            original_structure_filename = modbase13.get_coord_file(args.modbase)
            if os.path.exists(original_structure_filename):
                structure = load_structure(args.modbase, original_structure_filename)
            else:
                ddg_repo.psb_status_manager.sys_exit_Failure(
                   "Modbase model id %s not found in modbase2013/2016 directories"%\
                        args.modbase)

        if not DDG_monomer.evaluate_modbase(original_structure_filename):
            ddg_repo.psb_status_manager.sys_exit_failure("Terminating as %s of insufficient quality for ddg monomer"%args.modbase)


    if not os.path.exists(original_structure_filename):
        exit_str = "%s: No structure file %s"%(ddg_repo.structure_id, original_structure_filename)
        ddg_repo.psb_status_manager.sys_exit_failure(exit_str)

    # Open a file handle to our structure
    # LOGGER.debug("Loading orignal structure %s", original_structure_filename)
    LOGGER.info("Loading structure id=%s from %s"%(structure_id,original_structure_filename))
    structure_fin = None
    if original_structure_filename.endswith(".gz"):
        structure_fin = gzip.open(original_structure_filename, 'rt')
    else:
        structure_fin = open(original_structure_filename, 'rt')

    mmcif_dict = None
    if original_structure_type == 'mmCIF':
        mmCIF_parser = MMCIFParser() # QUIET=True)
        structure = mmCIF_parser.get_structure(structure_id, structure_fin)
        structure_fin.close()

        # DDG_Cleaner uses the mmcif information to remove non-species tags in PDB files
        mmcif_dict = mmCIF_parser._mmcif_dict  # << This is a little dangerous but the get_structure code above looks quite clean

        # Gather resolution information, etc, from the 
        structure_info_dict = mine_structure_info_from_mmCIF(mmcif_dict)
    else:
        pdb_parser = PDBParser()
        structure = pdb_parser.get_structure(structure_id,structure_fin)
        # If we ever return to deposited PDBs for ddg calculations (we now use mmcif), then this has to
        # be worked out.  For now, the only .pdb files we are processing are models that will not have
        # spurious residues that need to be cleaned out
        mmcif_dict = None


    ddg_cleaner = ddg_clean.DDG_Cleaner(structure,mmcif_dict,True,species_filter=args.species_filter)
    cleaned_structure, residue_to_clean_xref = ddg_cleaner.clean_structure_for_ddg(False, [args.chain])
        
    from Bio.PDB import PDBIO
    pdbio = PDBIO()
    cleaned_structure_temp_filename = None
    with tempfile.NamedTemporaryFile(delete=False,dir=ddg_structure_dir,mode="w") as cleaned_structure_temp:
        pdbio.set_structure(cleaned_structure)
        cleaned_structure_temp_filename = cleaned_structure_temp.name
        pdbio.save(cleaned_structure_temp_filename , write_end=True, preserve_atom_numbering=True)

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
