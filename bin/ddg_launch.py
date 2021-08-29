#!/usr/bin/env python
#
# Project        : ddg_repo
# Filename       : ddg_launch.py
# Authors        : Chris Moth
# Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : chris.moth@vanderbilt.edu
# Date           : 2020-August
# =============================================================================#

"""ddg_launch
   - Create and launch slurm files to call ddg_run.py on all (or many) variants in a PDB structure
   - Strips (preps) mmCIF or PDB input structures for Rosetta format, leaving
     behind a translation dictionary to enumaerate all residues on which calculations can be formed,
     and to map original residues numbering (and insertion codes) to Rosetta simplified 1..N residues)

   References: Prediction of protein mutational free energy: benchmark and
               sampling improvements increase classification accuracy
               Brandon Frenz, Steven Lewis, Indigo King, Hahnbeom Park, Frank DiMaio, Yifan Song

   """

import os
import grp
import sys
import logging
from logging.handlers import RotatingFileHandler
import pprint
import tempfile
import datetime

# from time import strftime
from typing import Dict,Tuple,List

from Bio.SeqUtils import seq1

from Bio.PDB import MMCIFParser
from Bio.PDB import Residue
from Bio.PDB import PDBParser
from Bio import BiopythonWarning
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from configparser import SafeConfigParser
from psb_shared.ddg_load_structure import ddg_load_structure,LoadStructureError

from psb_shared import psb_config
from psb_shared.ddg_repo import DDG_repo
from psb_shared import ddg_clean
from psb_shared.ddg_monomer import DDG_monomer
# from # psb_shared.ddg_monomer import DDG_cartesian
from psb_shared.psb_progress import PsbStatusManager

from slurm import slurm_submit

from lib import PDBMapSwiss
from lib import PDBMapModbase2016
from lib import PDBMapModbase2013


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
ch.setLevel(logging.DEBUG)

ddg_flavor='ddg_monomer'

cmdline_parser = psb_config.create_default_argument_parser(__doc__,
                                                           os.path.dirname(os.path.dirname(__file__)),
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
        "--pdbres", type=str,default='all',
        help="Either 'all' or a comma separate list of NNNA specific residues where A is an optional pdb insertion code")

cmdline_parser.add_argument(
            "--species_filter", type=str, required=False, default=None,
            help="Optionally retain only residues with DBREF of, example 'HUMAN'")

cmdline_parser.add_argument(
            "-s", "--residues_per_slurm", type=int, required=False, default=1,
            help="Increase to reduce the number of slurm files submitted, and incase the array jobs per slurm file")

cmdline_parser.add_argument("-n", "--nosbatch",
                            help="Create slurm files, but do NOT exec sbatch to launch them.", default=False,
                            action="store_true")


LOGGER.info("Command: %s" % ' '.join(sys.argv))

args = cmdline_parser.parse_args()


ddg_launch_log_filename = "ddg_launch.log"

# Add a rotating file handler log in the calculation directory
need_roll = os.path.isfile(ddg_launch_log_filename)

logger_fh = RotatingFileHandler(ddg_launch_log_filename, backupCount=5)

logger_fh.setFormatter(log_formatter)
logger_fh.setLevel(logging.DEBUG if args.debug else logging.INFO)
LOGGER.addHandler(logger_fh)

if need_roll:
    logger_fh.doRollover()

if args.debug:
    LOGGER.setLevel(logging.DEBUG)
    ch.setLevel(logging.DEBUG)







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


# We make it anyway below residue_to_clean_xref = ddg_repo.residue_to_clean_xref

slurm_directory = ddg_repo.slurm_directory_makedirs()


slurm_stdout_directory = os.path.join(slurm_directory,"out")
os.makedirs(slurm_stdout_directory,exist_ok=True)

structure_info_dict = {}

try:
    structure_id, structure, structure_info_dict, mmcif_dict = ddg_load_structure(args, config_dict)
except LoadStructureError as exception:
    sys.exit(exception.message)

# if residue_to_clean_xref:
#    structure_info_dict = ddg_repo.structure_config
#    LOGGER.info("ddg xref dictionary read from prior run")
# else:
ddg_cleaner = ddg_clean.DDG_Cleaner(structure, mmcif_dict, True, species_filter=args.species_filter)
cleaned_structure, residue_to_clean_xref = ddg_cleaner.clean_structure_for_ddg([args.chain])

# Not sure we need to save the cleaned file at this point.
# I suppose it's good to know if we're going to have a problem before we launch all the slurm scripts
from Bio.PDB import PDBIO

pdbio = PDBIO()
cleaned_structure_temp_filename = None
with tempfile.NamedTemporaryFile(delete=False, dir=ddg_structure_dir, mode="w") as cleaned_structure_temp:
    pdbio.set_structure(cleaned_structure)
    cleaned_structure_temp_filename = cleaned_structure_temp.name
    pdbio.save(cleaned_structure_temp_filename, write_end=True, preserve_atom_numbering=True)

# Save the cleaned structure in the repository
ddg_repo.mv_cleaned_structure_in(cleaned_structure_temp_filename)

# Add our cross reference to the repository
ddg_repo.residue_to_clean_xref = residue_to_clean_xref
ddg_repo.structure_config = structure_info_dict

# Resume by figuring out that we should keep everything or whatever
# this case of a SHEEP protein...
cleaned_structure_filename = ddg_repo.cleaned_structure_filename

try:
  slurmParametersAll = dict(config.items("SlurmParametersAll"))
except Exception as ex:
  slurmParametersAll = {}

import copy
slurmParametersDDGMonomer = copy.deepcopy(slurmParametersAll)

try:
  slurmParametersDDGMonomer.update(dict(config.items("SlurmParametersDDGMonomer")))
except Exception as ex:
  pass

slurm_required_settings = ['mem', 'account', 'ntasks', 'time']
for slurmSettingsDesc, slurmSettings in [('slurmParametersDDGMonomer', slurmParametersDDGMonomer)]:
    for req in slurm_required_settings:
        if not req in slurmSettings:
            LOGGER.error(
                "Can't launch jobs because you have not provided the required\nslurm setting [%s] for section %s (or slurmParametersAll)\n in either file %s or %s\n" %
                (req, slurmSettingsDesc, args.config, args.userconfig))
            sys.exit(1)

# print "Command Line Arguments"
LOGGER.info("Command Line Arguments:\n%s" % pprint.pformat(vars(args)))
LOGGER.info("slurmParametersDDGMonomer:\n%s" % pprint.pformat(slurmParametersDDGMonomer))

def slurm_file_create(residue_id_list: List[Tuple[str,int,str]], resname_list: List[str]) -> int:
    """

    :residues_in_slurm_file is a List of pairs, each containing:
        1) Residue id from original (not cleaned) structure that is to have ddg run
        2) Three letter standard amino acid, as standardized by the cleaner
    :return:         launched slurm jobid
    """

    first_residue_id = residue_id_list[0]
    pdb_residue_range_for_jobname = str(first_residue_id[1]) + str(first_residue_id[2]).strip()

    # Make the slurm job filename, and jobname, include the end residue
    # if more than one residue per slurm file
    if len(residue_id_list) > 1:
        last_residue_id = residue_id_list[-1]
        pdb_residue_range_for_jobname += '_' + str(last_residue_id[1]) + str(last_residue_id[2]).strip()


    # Write out the slurm file meta/docs
    slurm_filename = os.path.join(ddg_repo.slurm_dir,"%s_%s.slurm"%(structure_id,pdb_residue_range_for_jobname))
    with open(slurm_filename, 'w') as slurmf:
        slurmf.write("""\
#!/bin/bash
#
# Project        : ddG Repository
# Filename       : %s
# Generated By   : %s
# Organization   : Vanderbilt Genetics Institute,
#                : Program in Personalized Structural Biology,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University
# Generated on   : %s
# Description    : Generate a slurm file for each position in a structure file
#                  which executes ddg calculations on all possible variants
#===============================================================================
# Slurm Parameters
""" % (slurm_filename, __file__, str(datetime.datetime.now())))

        slurmf.write("\n")

        # It is important to make shallow copies of the default slurm confugration dictionaries, 
        # else mutations #7 can pick up dictionary entries set by mutation #4
        # slurmDict contains all the entries that will fill the header of our generated .slurm files
        slurmDict = dict(slurmParametersDDGMonomer)

        slurmDict['output'] = os.path.join(
            slurm_stdout_directory,
            "%s_%s_%%A_%%a.out" % (structure_id,pdb_residue_range_for_jobname))

        slurmDict['job-name'] = "%s_%s_%s" % (ddg_flavor,structure_id, pdb_residue_range_for_jobname)

        jobCount = len(residue_id_list) * 20 # Launch all 19 + 1 do nothing ddG calculations for each of N residues in this .slurm file
        slurmDict['array'] = "0-%d" % (jobCount - 1)


        # We can set a max on the number of simultaneously running jobs in the array.  For now, don't do this
        # if jobCount > 10:
        #   slurmDict['array'] += "%%10"

        for slurm_field in ['job-name', 'mail-user', 'mail-type', 'ntasks', 'time', 'mem', 'account', 'output',
                            'array']:
            if slurm_field in slurmDict:
                slurmf.write("#SBATCH --%s=%s\n" % (slurm_field, slurmDict[slurm_field]))

        # The job outputs start with a dump of runtime job->node assignment info 
        slurmf.write("""
echo "SLURM_ARRAY_TASK_ID="$SLURM_ARRAY_TASK_ID
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
# echo "SLURMTMPDIR="$SLURMTMPDIR
echo "SLURM_SUBMIT_DIR = "$SLURM_SUBMIT_DIR
""")

        all_amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
        slurmf.write("""
# For each n Job Array # 0 to 19, start ddg_run.py on the nth amino acid variant   
amino_acids=(%s)
"""%' '.join(all_amino_acids))

        bash_formatted_residue_ids = [
                seq1(resname) + str(resid[1])+resid[2].strip() for (resid,resname) in zip(residue_id_list,resname_list)]

        slurmf.write("""
# AminoAcid, PDB residue number, optional insertion code of each residues on which ddG is being run
pdb_positions=(\n %s)
"""%' '.join([bash_formatted_residue_id + '\n' * (1 if (i+1)%10 == 0 else 0) for i,bash_formatted_residue_id in enumerate(bash_formatted_residue_ids)]))


        slurmf.write("""
residue_subscript=$((SLURM_ARRAY_TASK_ID/20))
amino_acid_subscript=$((SLURM_ARRAY_TASK_ID%20))

pdb_position=${pdb_positions[$residue_subscript]:1}
native_aa=${pdb_positions[$residue_subscript]:0:1}
variant_aa=${amino_acids[$amino_acid_subscript]}

variant=$native_aa$pdb_position$variant_aa

echo Running ddG for variant $variant
""")

        slurmf.write("""
source %s
cmd="ddg_run.py --%s=%s --chain=%s --variant=$variant --species=human"
echo Running $cmd
`$cmd`

if [ $? != 0 ]; then
echo Failure to run:
echo  $cmd 
exit 1
fi
""" % ( os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "psb_prep.bash"),
        ddg_repo.structure_source,  # Will fill in '--pdb' or '--swiss' etc
        structure_id,
        args.chain,
        ))


    LOGGER.info("Slurm file created: %s", slurm_filename)

    jobid = -1
    if args.nosbatch:
        LOGGER.warning("sbatch will not be called to launch %s because args.nosbatch=%s", slurm_filename,
                       args.nosbatch)
    else:
        jobid = slurm_submit(['sbatch', slurm_filename])
        try:
            jobid = int(jobid)
        except ValueError:
            jobid = -1
            
        LOGGER.info("sbatch %s returned jobid=%d  WARNING CLEARING STATUS DIRECTORIES",slurm_filename,jobid)

        for bash_formatted_residue_id in bash_formatted_residue_ids:
            native_aa_letter = bash_formatted_residue_id[0]
            for variant_aa in all_amino_acids:
                if variant_aa != native_aa_letter:
                    ddg_repo.set_variant("%s%s"%(bash_formatted_residue_id,variant_aa))
                    ddg_repo.psb_status_manager.clear_status_dir()
                    ddg_repo.psb_status_manager.write_info('Launched Job %d'%jobid)

    return jobid

residue_list = []
if args.pdbres.lower() == 'all':
    residue_list = [residue.id for residue in structure[0][args.chain]]
else:
    raw_residue_list = args.pdbres.split(',')
    for raw_residue in raw_residue_list:
        if raw_residue[-1].isalpha():
            residue_id = (' ',int(raw_residue[0:-1]),raw_residue[-1])
        else:
            residue_id = (' ',int(raw_residue),' ')
        residue_list.append(residue_id)
    # If the residue is in the cleaned protein submitted
    # to rosetta - then we need to process all variants in it

residue_ids_in_slurm_file = []
resnames_in_slurm_file = []
for residue_id in residue_list:
    if residue_id in residue_to_clean_xref:
        cleaned_residue_id = residue_to_clean_xref[residue_id]
        residue_ids_in_slurm_file.append(residue_id)
        resnames_in_slurm_file.append(cleaned_structure[0][args.chain][cleaned_residue_id].get_resname())
        if len(residue_ids_in_slurm_file) == args.residues_per_slurm:
            slurm_file_create(residue_ids_in_slurm_file,resnames_in_slurm_file)
            residue_ids_in_slurm_file = [] # Clear out the list for the next slurm file
            resnames_in_slurm_file = [] # Clear out the list for the next slurm file
    else:
        if residue_id[0] != 'W': # No need to report on skipping waters
            LOGGER.info("DDG will not be run on residue %s"%str(residue_id))

assert len(residue_ids_in_slurm_file) == len(resnames_in_slurm_file)
if len(residue_ids_in_slurm_file) > 0:
    slurm_file_create(residue_ids_in_slurm_file,resnames_in_slurm_file)

