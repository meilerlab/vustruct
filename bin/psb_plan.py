#!/usr/bin/env python3
#
# Project        : PDBMap
# Filename       : psb_plan.py
# Authors        : R. Michael Sivley
#project Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2017-01-22
# Description    : Identifies available structures for a specific mutation
#                : and generates a table of work for the computer cluster
#=============================================================================#

"""\
Given a case on the command line, and either:

- a UDN project ID (ex. UDN123456) which will cause reading of all mutants in:
  ..../UDN/UDN124356/UDN123456_missense.csv

or

- single gene name or UniProt identifier, refseq transcript ID, and mutation

this script will analyze structural coverage and genetic datasets to 
generate a detailed execution plan for analyzing various structural 
and spatial properties of the available structural and genetic datasets.

The three output files for each gene include:
  structure_report.csv
  dropped_structures.csv
  workplan.csv
"""

print("%s: Pipeline execution plan generator.  -h for detailed help"%__file__)

import os
import sys
import grp
import stat
import pwd
import gzip
import lzma
import logging
import datetime

from typing import Dict

from logging.handlers import RotatingFileHandler
from logging import handlers

sh = logging.StreamHandler()
LOGGER = logging.getLogger()
global_fh = None # Setup global fileHandler log when we get more info

# Now that we've added streamHandler, basicConfig will not add another handler (important!)
log_format_string = '%(asctime)s %(levelname)-4s [%(filename)16s:%(lineno)d] %(message)s'
date_format_string = '%H:%M:%S'
log_formatter = logging.Formatter(log_format_string,date_format_string)

def gene_specific_set_log_formatter(log_handler: logging.Handler, entry: int, gene:str) -> None:
    log_format_string = '%3d %-7s'%(entry,gene) + ' %(asctime)s %(levelname)-4s [%(filename)16s:%(lineno)d] %(message)s'
    date_format_string = '%H:%M:%S'
    
    log_handler.setFormatter(logging.Formatter(log_format_string,date_format_string))

LOGGER.setLevel(logging.DEBUG)
sh.setLevel(logging.INFO)
sh.setFormatter(log_formatter)
LOGGER.addHandler(sh)

import pandas as pd
pd.options.mode.chained_assignment = 'raise'
pd.set_option("display.max_columns",100)
pd.set_option("display.width",1000)
import numpy as np
import time
from get_structures import get_pdbs,get_modbase_swiss
import argparse,configparser
from glob import glob
from lib import PDBMapSIFTSdb
from lib import PDBMapProtein
from lib import PDBMapSwiss
# from lib import PDBMapModbase2016
# from lib import PDBMapModbase2013
from lib import PDBMapModbase2020
from lib.PDBMapTranscriptUniprot import PDBMapTranscriptUniprot
from lib.PDBMapTranscriptEnsembl import PDBMapTranscriptEnsembl
from lib.PDBMapAlignment import PDBMapAlignment
from lib.amino_acids import longer_names
import pprint
import shutil
import warnings
# from jinja2 import Environment, FileSystemLoader

from psb_shared import psb_config

cmdline_parser = psb_config.create_default_argument_parser(__doc__,os.path.dirname(os.path.dirname(__file__)))

cmdline_parser.add_argument("-m","--maybe",type=int,default=10,metavar='maybe_threshold',
help='Sequence distance "outside of structure" threshold, below which structures which lack coverage are marked as coverage="Maybe" instead of "No"')
# cmdline_parser.add_argument("collaboration",type=str,help="Collaboration ID (ex. UDN)")
cmdline_parser.add_argument("project",type=str,help="Project ID (ex. UDN124356)",default=os.path.basename(os.getcwd()),nargs='?')
cmdline_parser.add_argument("gene",nargs='?',type=str,help="optional Gene ID (ex TTN)")
cmdline_parser.add_argument("refseq",nargs='?',type=str,help="Omit OR NM_... refseq transcript identifier")
cmdline_parser.add_argument("mutation",nargs='?',type=str,help="Omit OR HGVS mutation string (ex S540A)")
args,remaining_argv = cmdline_parser.parse_known_args()

infoLogging = False

if args.debug:
  infoLogging = True
  sh.setLevel(logging.DEBUG)
elif args.verbose:
  infoLogging = True
  sh.setLevel(logging.INFO)
# If neither option, then logging.WARNING (set at top of this code) will prevail for stdout

# If we can reload from RAM, it saves a ton of time
from copy import deepcopy
MMCIF_DICT_CACHE = {}
MMCIF_STRUCTURE_CACHE = {}

required_config_items = ["dbhost","dbname","dbuser","dbpass",
  "collaboration",
  "idmapping",
  "interpro_dir",
  "pdb_dir",
  "sec2prim",
  # "modbase2013_dir",
  # "modbase2013_summary",
  # "modbase2016_dir",
  # "modbase2016_summary",
  "modbase2020_dir",
  "modbase2020_summary",
  "output_rootdir",
  "swiss_dir",
  "swiss_summary"]

# parser.add_argument("udn_excel",type=str,help="Raw input UDN patient report (format: xls/xlsx)")
# parser.add_argument("udn_csv",type=str,help="Parsed output pipeline filename (e.g. filename.csv)")

config,config_dict = psb_config.read_config_files(args,required_config_items)

keep_pdbs = []
if 'KeepPDBs' in config:
    keep_pdbs = dict(config.items('KeepPDBs'))['keeppdbs'].upper().split(',')

from psb_shared import psb_perms
psb_permissions = psb_perms.PsbPermissions(config_dict)

# The collaboration_dir is the master directory for the case, i.e. for one patient
# Example: /dors/capra_lab/projects/psb_collab/UDN/UDN532183
udn_root_directory = os.path.join(config_dict['output_rootdir'],config_dict['collaboration'])
collaboration_dir = os.path.join(udn_root_directory,args.project)
psb_permissions.makedirs(collaboration_dir)
# collaboration_log_dir = os.path.join(collaboration_dir,"log")

# Whether we made it above or not, we want our main directory for this run to be group capra_lab and sticky!
psb_permissions.set_dir_group_and_sticky_bit(collaboration_dir)

# psb_permissions.makedirs(collaboration_log_dir)

LOGGER.info("Loading UniProt ID mapping...")

"""
# This detect_entity code seems entirely unused in 2019
def detect_entity(io,entity):
  LOGGER.info("\nAttempting to identify entity: %s..."%entity)
  etype = io.detect_entity_type(entity)
  if etype == "unp":
    unp  = entity
    hgnc = PDBMapProtein.unp2hgnc(unp)
    LOGGER.info("UniProt: %s => HGNC: %s..."%(unp,hgnc))
  elif etype and etype!="structure" and etype!="model": # HGNC gene ID
    hgnc = entity
    unp  = etype
    LOGGER.info("HGNC: %s => UniProt: %s..."%(hgnc,unp))
  else:
    LOGGER.error("TERMINATING: Unrecognized ID %s in detect_entity(). Gene name or UniProt/Swiss-Prot AC required."%entity)
    sys.exit(1)
  return unp,hgnc
"""


config_dict_reduced = {x:config_dict[x] for x in required_config_items}

config_dict = config_dict_reduced

config_dict_shroud_password = {x:config_dict[x] for x in required_config_items}
dbpass = config_dict.get('dbpass','?')
config_dict_shroud_password['dbpass'] = '*' * len(dbpass)

oneMutationOnly = args.gene or args.refseq or args.mutation

# If this is a "global" run of psb_plan, then make a global log file in the collaboration directory
if not oneMutationOnly: 
  # pw_name = pwd.getpwuid( os.getuid() ).pw_name # example jsheehaj or mothcw
  log_filename = os.path.join(collaboration_dir,"%s_psb_plan.log"%args.project)

  sys.stderr.write("Collaboration-wide  psb_plan log file is %s\n"%log_filename)
  needRoll = os.path.isfile(log_filename)

  global_fh = RotatingFileHandler(log_filename, backupCount=7)
  formatter = logging.Formatter('%(asctime)s %(levelname)-4s [%(filename)20s:%(lineno)d] %(message)s',datefmt="%H:%M:%S")
  global_fh.setFormatter(formatter)
  global_fh.setLevel(logging.INFO)
  LOGGER.addHandler(global_fh)

  if needRoll:
    global_fh.doRollover()

# print "Command Line Arguments"
# pprint.pprint(vars(args))
LOGGER.info("Command Line Arguments:\n%s"%pprint.pformat(vars(args)))
LOGGER.info("Configuration File parameters:\n%s"%pprint.pformat(config_dict_shroud_password))

config_pathprox_dict = dict(config.items("PathProx"))  
LOGGER.info("Pathprox config entries:\n%s"%pprint.pformat(config_pathprox_dict))

def pathprox_config_to_argument(disease1_or_2_or_neutral,default):
  variants_filename_key = "%s_variant_filename"%disease1_or_2_or_neutral
  if variants_filename_key in config_pathprox_dict:
    config_str = config_pathprox_dict[variants_filename_key]
    variants_pathprox_args = ''
    if 'neutral' in variants_filename_key:
       variants_pathprox_args = "--neutral %s --neutral_label %s"%(config_str,config_pathprox_dict["%s_variant_sql_label"%disease1_or_2_or_neutral])
    else:
       variants_pathprox_args = "--pathogenic %s --pathogenic_label %s"%(config_str,config_pathprox_dict["%s_variant_sql_label"%disease1_or_2_or_neutral])
    return variants_pathprox_args
  else:
    variants_sql_label_key = "%s_variant_sql_label"%disease1_or_2_or_neutral
    if variants_sql_label_key not in config_pathprox_dict:
      LOGGER.critical("A configuration file must define a variant set for %s or %s in the [PathProx] section"%(variants_sql_label_key,variants_filename_key))
    config_str = config_pathprox_dict[variants_sql_label_key]
    if config_str == "clinvar38":
      return "--add_pathogenic38"
    if config_str == "clinvar":
      return "--add_pathogenic"
    if config_str == "tcga":
      return "--add_tcga"
    if config_str == "cosmic":
      return "--add_cosmic"
    if config_str == "cosmic38":
      return "--add_cosmic38"
    if config_str == "1kg3":
      return "--add_1kg3"
    if config_str == "gnomad":
      return "--add_gnomad"
    if config_str == "gnomad38":
      return "--add_gnomad38"
    if config_str == "exac":
      return "--add_exac"
    if config_str == "ERROR":
      return "BAD VALUE IN PATHPROX CONFIG FILE SECTION OR CALLER"
    LOGGER.critical("pathprox variant config of [%s] was not found. Applying default of [%s]",config_str,default)
    return pathprox_config_to_argument(default,"ERROR")

pathprox_arguments = {
  'disease1': pathprox_config_to_argument('disease1',"clinvar"),
  'disease2': pathprox_config_to_argument('disease2',"COSMIC"),
  'neutral' : pathprox_config_to_argument('neutral',"exac")
}

LOGGER.info("Pathprox Disease 1 command line argument: %s"%pathprox_arguments['disease1'])
LOGGER.info("Pathprox Disease 2 command line argument: %s"%pathprox_arguments['disease2'])
LOGGER.info("Pathprox Neutral   command line argument: %s"%pathprox_arguments['neutral'])


LOGGER.info("Results for patient case %s will be rooted in %s"%(args.project,collaboration_dir))

swiss_filename = os.path.join(config_dict['swiss_dir'],config_dict['swiss_summary'])
if not infoLogging:
  print("Loading swiss model JSON metadata from %s"%swiss_filename)
LOGGER.info("Loading swiss model JSON metadata from %s"%swiss_filename)

# import cProfile
# cProfile.run("PDBMapSwiss.load_swiss_INDEX_JSON(config_dict['swiss_dir'],config_dict['swiss_summary'])")
# sys.exit(0)
PDBMapSwiss.load_swiss_INDEX_JSON(config_dict['swiss_dir'],config_dict['swiss_summary'])
# In order to set the directories, we need the gene name worted out, interestingly enough
if not infoLogging:
  print("Loading idmapping file from %s"%config_dict['idmapping'])

LOGGER.info("Loading idmapping")
PDBMapProtein.load_idmapping(config_dict['idmapping'])
LOGGER.info("Loading done")

# print "Establishing connection with the PDBMap database..."
# if (not io):
#  io = PDBMapIO(config_dict['dbhost'],config_dict['dbuser'],config_dict['dbpass'],config_dict['dbname']) # ,slabel=args.slabel)

def repr_subset(ci_df_excerpt):
    """Discard models and pdbs that duplicate sequence coverage
       and prioritize pdbs over models.  Note, the way this is called, this is
       not really happening."""
    # Sort order: Mutation Coverage, PDBs, Length (rounded), Resolution/Identity
    tdf = ci_df_excerpt.copy()
    # Convert these sequence start/end fields to ints
    tdf['rndlen'] = 1+tdf['trans_last_aligned'] - tdf['trans_first_aligned']
    tdf['mutcov'] = tdf['distance_from_aligned']==0
    tdf['ispdb']  = tdf['label'] == 'pdb'
    tdf = tdf.sort_values(by=['mutcov','ispdb','rndlen','resolution','template_identity'],ascending=[False,False,False,True,False])
    # Create the empty dataframe of minimally overlapping (quasi) structures
    # with the same column names and column types as the input
    nr_cov = pd.DataFrame(columns=tdf.columns).astype(tdf.dtypes.to_dict())
    cov = set([])
    template_set = set([])
    for _,row in tdf.iterrows():
      s_cov = list(range(int(row['trans_first_aligned']),1+int(row['trans_last_aligned'])))
      # Test if <10% of residues in *this model* overlap current coverage
      if row['pdb_template'] is None:
        templatePDB = row['structure_id'] 
      else:
        templatePDB = row['pdb_template']
      if len(templatePDB) > 4:
        templatePDB = templatePDB[:4]
      # if ( (templatePDB and (not templatePDB in template_set)) or
      #     (not cov or len(set(s_cov).intersection(cov))/float(len(s_cov))<0.1)):
      if True:
        nr_cov = nr_cov.append(row)
        # Add the sequence range to the cumulative coverage
        cov |= set(s_cov)
        template_set |= set([templatePDB])
    # Return the minimally overlapping subset of structures
    return nr_cov

def get_pdb_pos(*args):
  label,sid,chain,transcript,trans_pos,d2b = args[0].values
  if d2b:
    # PDB does not cover this residue. No associated PDB position.
    return np.nan
  global io

  Pdb_sql = """
  SELECT trans_seqid,chain_seqid
  FROM Alignment
  WHERE label = 'pdb' AND
        structid=%(sid)s AND chain=%(chain)s AND 
        transcript=%(trans)s AND trans_seqid=%(trans_pos)s
  """

  ModbaseSwiss_sql = """
  SELECT trans_seqid,actr.chain_res_num as chain_seqid
  FROM AlignChainTranscript AS act
  LEFT OUTER JOIN AlignChainTranscriptResidue actr ON actr.al_id = act.al_id
  WHERE label in ('modbase','swiss') AND
        act.structid=%(sid)s AND act.chain=%(chain)s AND 
        act.transcript=%(trans)s AND actr.trans_seqid=%(trans_pos)s
  """

  df = pd.DataFrame() # Init to empty dataframe
  if len(sid) <=4 or label == 'pdb':
    df = pd.read_sql(Pdb_sql,io._con,params={'sid':sid,'chain':chain,
                            'trans':transcript,'trans_pos':trans_pos})

  if df.empty:
    df = pd.read_sql(ModbaseSwiss_sql,io._con,params={'sid':sid,'chain':chain,
                            'trans':transcript,'trans_pos':trans_pos})
    if df.empty:
      LOGGER.info("WARNING: Residue %s is missing in %s %s.%s."%(trans_pos,label,sid,chain))
      return np.nan

  return df["chain_seqid"].values[0]

def makejob(flavor,command,params,options: str,cwd=os.getcwd()):
    job = {}
    job['flavor'] = flavor
    job['config'] = args.config
    job['userconfig'] = args.userconfig
    job['command'] = command
    job['options'] = options
    job['gene'] = params['gene']
    job['refseq'] = params['refseq']
    job['method'] = params.get('method','')
    job['pdbmut'] = params.get('pdb_mutation','')
    job['mutation'] = params['mutation']
    job['project'] = args.project
    job['unp'] = params['unp']
    
    job['pdbid'] = params.get('structure_id',params.get('pdbid',''))

    # If our structure name ends in .pdb, that's fine - but trim that out for the formation
    # of our output directory and unique key
    dot_pdb = job['pdbid'].find(".pdb")
    if dot_pdb > 2:
        job['pdbid'] = job['pdbid'][0:dot_pdb]
    job['chain'] = params.get('chain_id',params.get('chain',''))
    job['mers'] = params['mers']


    if "SequenceAnnotation" in flavor and job['pdbid'] == 'N/A':
        job['uniquekey'] = "%s_%s_%s_%s"%(job['gene'],job['refseq'],job['mutation'],job['flavor'])
        job['outdir'] = params['mutation_dir']
    elif "DigenicAnalysis" in flavor and job['pdbid'] == 'N/A':
        job['uniquekey'] = job['flavor']
        job['outdir'] = params['mutation_dir']
    else: # Most output directories have pdbid and chain suffixes
        if not job['chain'].strip() or job['chain'] == "''" or job['chain'] == "' '":
            pdb_chain_segment = job['pdbid']
        else:
            pdb_chain_segment = "%s_%s"%(job['pdbid'],job['chain'])
        if "ddG" in flavor:
            job['outdir'] = "ddG repository"
        else:
            job['outdir'] = os.path.join(params['mutation_dir'] , pdb_chain_segment )
        job['uniquekey'] = "%s_%s_%s_%s_%s"%(job['gene'],job['refseq'],job['mutation'],pdb_chain_segment,job['flavor'])
    job['cwd'] = cwd
    # too much detail LOGGER.debug(pprint.pformat(job))
    return pd.Series(job)

# Jobs for PathProx analyses
def makejobs_pathprox_df(params: Dict, ci_df: pd.DataFrame, multimer: bool) -> pd.DataFrame:
    pathprox_jobs_to_run_df = pd.DataFrame()
    _params = params
    for index,ci_row in ci_df.iterrows():  # For each chain_info record...
        params = _params.copy() # Dont affect the dictionary passed in  Shallow copy fine
        params['method'] = ci_row['method']
        params['structure_id'] = ci_row['structure_id']
        params['chain'] = "%s"%ci_row['chain_id']  # Quote it in case the IDentifier is a space, as in many models
        params['mers'] = ci_row['mers']
        params['pdb_mutation'] = ci_row['pdb_mutation']

        # sqlcache was a brief idea - no longer needed
        # params['sqlcache'] = os.path.join(params['mutation_dir'],"sqlcache")
        params['pathprox_permutations'] = 99999

        params_transcript = ''
        # Maybe we resurrect ENST transcript later.
        # However, pathprox should be more in uniprot space
        # Though, it needs complete transcript info to load variants from SQL
        # if 'transcript' in params and params['transcript']:
        #    params_transcript = " --isoform '%s'"%params['transcript']

        # ClinVa params['disease1argument'] = pathprox_disease1argument

        #  Make jobs for Disease 1 (maybe ClinVar) PathProx and Disease 2 (maybe COSMIC) PathProx
        disease_key_list = ['disease1', 'disease2' ]
        # sometimes we set disease 1 and 2 to same - that's OK - but don't run them both
        if config_pathprox_dict['disease1_variant_sql_label'] == config_pathprox_dict['disease2_variant_sql_label']:
            disease_key_list = ['disease1']
            LOGGER.warning("Only one Pathprox disease type sql label has been specified %s"%disease_key_list[0])
        for disease_1or2 in disease_key_list:
            disease_variant_sql_label = config_pathprox_dict['%s_variant_sql_label'%disease_1or2]
            params['diseaseArgument'] = pathprox_arguments[disease_1or2]
            params['neutralArgument'] = pathprox_arguments['neutral']
            pathprox_flavor="PP_%s_%s"%("complex" if multimer else "monomer",disease_variant_sql_label)
            structure_designation = ''
            if ci_row['label'] == 'usermodel':
                structure_designation = "--usermodel=%s --label=%s_%s "%(ci_row['struct_filename'],ci_row['structure_id'],params['chain'])
            elif ci_row['label'] == 'swiss':
                structure_designation = "--swiss=%s "%ci_row['structure_id'] # <- Pathprox will make label from this and --chain
            elif ci_row['label'] == 'modbase':
                structure_designation = "--modbase=%s "%ci_row['structure_id'] # <- Pathprox will make label from this and --chain
            else:
                structure_designation = "--%s=%s "%("biounit" if multimer else "pdb",params['structure_id'].lower())
            j = makejob(pathprox_flavor,"pathprox3.py",params,
                ("-c %(config)s -u %(userconfig)s --variants='%(chain)s:%(transcript_mutation)s' " +
                 "--chain%(chain)sunp='%(unp)s' " +
                 structure_designation +
                 ("--chain='%s' "%ci_row['chain_id'] if not multimer else "") + 
                 "%(diseaseArgument)s %(neutralArgument)s --radius=D " + 
                 ## deprecated -> "--sqlcache=%(sqlcache)s " + 
                 (("--permutations=%s "%params['pathprox_permutations']) if 'pathprox_permutations' in params else '') + 
                 "--overwrite")%params + params_transcript)
            pathprox_jobs_to_run_df = pathprox_jobs_to_run_df.append(j,ignore_index=True)

    return pathprox_jobs_to_run_df

def makejob_udn_sequence(params):
    # Secondary structure prediction and sequence annotation
    return makejob("SequenceAnnotation","udn_pipeline2.py",params,   #command
          ("--config %(config)s --userconfig %(userconfig)s " +
           "--project %(collab)s --patient %(project)s --gene %(gene)s --transcript %(transcript_mutation)s")%params)

def makejob_ddg_monomer(params):
    return makejob("ddG_monomer","ddg_run.py",params,   #command
          ("--config %(config)s --userconfig %(userconfig)s " +
           "--%(label)s %(structure_id)s " + 
           "--chain %(chain)s --variant %(pdb_mutation)s")%params) # options

# Save ourselves a lot of trouble with scoping of df_dropped by declaring it global with apology
# Fix this in python 3 - where we have additional local scoping options
df_dropped = None

def plan_one_mutation(index:int, gene: str,refseq: str,mutation: str,user_model:str =None,unp: str=None):
    # Use the global database connection
    global io
    global args
    global oneMutationOnly
    global df_dropped

    LOGGER.info("Planning mutation for %s %s %s %s",args.project,gene,refseq,mutation)
    if user_model:
        LOGGER.info("Additional User model %s requested",user_model)

    # 2017-10-03 Chris Moth modified to key off NT_ refseq id
    unps = PDBMapProtein.refseqNT2unp(refseq)
    if unps:
        newunp = unps[0]
        # If the refseq - looked-up-unp is materially different from the supplied unp, make some noise
        if unp and ( 
            (newunp.split('-')[0] != unp.split('-')[0])  # Clearly A12345 is NOT Same as J54321
            or # We know both match to first 6 positions....  Do their Uniparc IDs match
            PDBMapProtein.unp2uniparc(unp) != PDBMapProtein.unp2uniparc(newunp)):
                LOGGER.warning("unp %s determined from refseq %s will override provided %s",newunp,refseq,unp)
        unp = newunp
        gene = PDBMapProtein.unp2hgnc(unp)
        del unps
    else: # Sort this out the old way
        if unp:
            LOGGER.warning("Setting refseq to NA and using unp=%s",unp) 
            gene = PDBMapProtein.unp2hgnc(unp)
            refseq = 'NA'
        else:
            msg = "Neither refseq nor unp enable transcript identification"
            LOGGER.critical(msg)
            sys.exit(msg)
            LOGGER.warning("Refseq %s does not map to a uniprot identifier.  Attempting map of gene %s\n"%(refseq,gene))
            # io = PDBMapIO(config_dict['dbhost'],config_dict['dbuser'],config_dict['dbpass'],config_dict['dbname']) # ,slabel=args.slabel)
            # unp,gene = detect_entity(io,entity)

    if not gene:
        LOGGER.critical("A gene name was not matched to the refseq or uniprot ID.  Gene will be set to %s",gene)
    
    mutation_dir = os.path.join(collaboration_dir,"%s_%s_%s"%(gene,refseq,mutation))
    psb_permissions.makedirs(mutation_dir)
    mutation_log_dir = mutation_dir
    psb_permissions.makedirs(mutation_log_dir)

    # pw_name = pwd.getpwuid( os.getuid() ).pw_name # example jsheehaj or mothcw
    log_filename = os.path.join(mutation_log_dir,"psb_plan.log")


    if oneMutationOnly:
        sys.stderr.write("psb_plan log file is %s\n"%log_filename)
    else:
        LOGGER.info("Additionally logging to file %s"%log_filename)
    needRoll = os.path.isfile(log_filename)

    local_fh = RotatingFileHandler(log_filename, backupCount=7)
    formatter = logging.Formatter('%(asctime)s %(levelname)-4s [%(filename)20s:%(lineno)d] %(message)s',datefmt="%H:%M:%S")
    local_fh.setFormatter(formatter)
    local_fh.setLevel(logging.INFO)
    LOGGER.addHandler(local_fh)

    if needRoll:
        local_fh.doRollover()

    gene_specific_set_log_formatter(global_fh,index,gene)
    gene_specific_set_log_formatter(sh,index,gene)


    # LOGGER.info("MySQL: Marking psb_collab.MutationSummary incomplete for %s %s %s %s %s"%
    #         (args.project,gene,refseq,mutation,unp))
    
    # Ensure that this mutation is tracked in the database
    # Even though typically already added in Mtuation Summary
    """c = io._con.cursor()
    # No-op if this mutation is already recorded in the database
    sql  = "INSERT IGNORE INTO psb_collab.MutationSummary "
    sql += "(collab,project,gene,refseq,mutation,uniprot) "
    sql += "VALUES (%s,%s,%s,%s,%s,%s) "
    c.execute(sql,(config_dict['collaboration'],args.project,gene,refseq,mutation,unp))
    c.close()"""

    """ 
    # Update the database to include the UniProt AC for this gene
    # Update to incomplete status, in case this is a rerun of a failed report
    c   = io._con.cursor()
    sql  = "UPDATE psb_collab.MutationSummary SET uniprot=%s, complete=0, comments='' WHERE "
    sql += "collab=%s AND project=%s AND gene=%s AND refseq=%s AND mutation=%s"
    c.execute(sql,(unp,config_dict['collaboration'],args.project,gene,refseq,mutation))
    c.close()
    """
    
    ################################################
    # Load relevant PDBs for consideration         # 
    ################################################
    pdbs_chains = PDBMapSIFTSdb.pdbs_chains_for_all_unp_isoforms(unp)

    LOGGER.info("%d Deposited PDB Chains are for unp=%s"%(len(pdbs_chains),unp.split('-')[0]))

    transcript = PDBMapTranscriptUniprot(unp)
    unp_is_canonical = PDBMapProtein.isCanonicalByUniparc(unp)

    from Bio.PDB import MMCIFParser
    from Bio.PDB import PDBParser
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from Bio.SeqUtils import seq1
    mmCIF_parser = MMCIFParser(QUIET=True)
    pdb_parser = PDBParser(QUIET=True)
    structures_dict = {}
    mmCIFDicts_dict = {}
    biounits_dict = {}

    # Load the pdbs and biounits from filesystem
    # Biosystems give us the best 'complex' multi-chain picture
    # Asymetric units are great for single chain processing, and for
    # Backstory on method, resolution, etc
    for pdb_upper_id,chain_id in pdbs_chains:
        pdb_id = pdb_upper_id.lower()
        biounit = None
        structure = None

        # Don't reload asymetric unit or biounit if we already loaded before for another chain
        if pdb_id in structures_dict:
            structure = structures_dict[pdb_id]
            # We would have also loaded a biounit when we grabbed that pdbid
            if pdb_id in biounits_dict:
                biounit = biounits_dict[pdb_id]
        else:
            # We've a new pdb_id to sort out.
            # to my great annoyance, Biounits are typically only available in pdb format - so we try that first
            biounit_filename = os.path.join(config_dict['pdb_dir'],'biounit',"PDB",'divided',pdb_id[1:3],"%s.pdb1.gz"%pdb_id)
            if os.path.exists(biounit_filename):
                LOGGER.debug("Loading %s",biounit_filename)
                if biounit_filename.endswith(".gz"):
                    biounit_fin = gzip.open(biounit_filename,'rt')
                else:
                    biounit_fin = open(biounit_filename,'rt')
                biounit =  pdb_parser.get_structure(pdb_id,biounit_fin)
                biounit_fin.close()
                biounits_dict[pdb_id] = biounit
            else: # _Maybe_ the biounit is in an mmCIF file (big CRYO-em for example)
                biounit_filename = os.path.join(config_dict['pdb_dir'],'biounit',"mmCIF",'divided',pdb_id[1:3],"%s-assembly1.cif.gz"%pdb_id)
                if os.path.exists(biounit_filename):
                    LOGGER.debug("Loading %s",biounit_filename)
                    if biounit_filename.endswith(".gz"):
                        biounit_fin = gzip.open(biounit_filename,'rt')
                    else:
                        biounit_fin = open(biounit_filename,'rt')
                    biounit =  mmCIF_parser.get_structure(pdb_id,biounit_fin)
                    biounit_fin.close()
                    biounits_dict[pdb_id] = biounit
                else:
                    LOGGER.debug("No biounit PDB or mmCIF file found for pdbid=%s",pdb_id)

            # All assymetric units are available in mmCIF format - so grab that
            if pdb_id in MMCIF_DICT_CACHE:
                mmCIFDicts_dict[pdb_id] = deepcopy(MMCIF_DICT_CACHE[pdb_id])
                structures_dict[pdb_id] = MMCIF_STRUCTURE_CACHE[pdb_id].copy()
            else:
                structure_filename = os.path.join(config_dict['pdb_dir'],'structures','divided','mmCIF',pdb_id[1:3],"%s.cif.gz"%pdb_id)
                if os.path.exists(structure_filename):
                    LOGGER.debug("Loading %s",structure_filename)
                    if structure_filename.endswith(".gz"):
                        structure_fin = gzip.open(structure_filename,'rt')
                    else:
                        structure_fin = open(structure_filename,'rt')
                    structure =  mmCIF_parser.get_structure(pdb_id,structure_fin)
                    structure_fin.close()
                    mmCIF_dict = mmCIF_parser._mmcif_dict # << This is a little dangerous but the get_structure code above looks quite clean
                    mmCIFDicts_dict[pdb_id] = mmCIF_dict
                    MMCIF_DICT_CACHE[pdb_id] = deepcopy(mmCIF_dict)
                    MMCIF_STRUCTURE_CACHE[pdb_id] = structure.copy()
                    structures_dict[pdb_id] = structure
                else:
                    LOGGER.critical("%s: No structure file %s",pdb_id,structure_filename)

    # Life is much easier if we have a class that knows 'everything' relevant about ranking pdb structures
    class ChainInfo(object):
        def __init__(self,_label: str = 'N/A'):
            self.label = _label        # Typically pdb, swiss, modbase, usermodel
            self.structure_id = None   # unique ID (filename base) of the structure within the domain of pdb/swiss/modbase
            self.structure_url = None  # url for more information from rcsb, swissmodel, etc
            self.struct_filename = None# Required for user models as argument to pathprox.  Nice to gather as we can for others
            self.mers = 'monomer'      # Could be dimer/trimer/etc - Processing of an entire complex is a top level idea, with the other chains peer calculations
            self.pdb_template = None   # For models, the deposited PDB structure upon which the model was created
            self.template_identity = None # For models, the Swiss/Modbase provided seq identity number between transcript and pdb template
            self.template_chain = None # The chain ID of the template primarily used to build the model
            self.chain_id = None       # The .cif chain letter A/B/C.... possibly #s and lower letters for large cryoEM
            self.biounit = None        # Boolean to record whether a biounit file is available
            self.method = None         # X-RAY...  CRYO-EM, etc
            self.resolution = None     # Resolution in Angtroms, method-specific.  Also resolution of model templates, if available
            self.deposition_date = None# YYYY-MM-DD format
            self.biounit_chains = None # count of chains in the biounit
            self.nresidues = None      # Number of residues in the biounit or structure
            self.perc_aligned = None   # percent of transcript aligned to the PDB     PDBMapAlignment.perc_aligned 
            self.perc_identity = None  # percent of identical aligned residues        PDBMapAlignment.perc_identity
            self.aln_score = None      # Blosum based score summed for all aligned residues  PDBMapAlignment.aln_score
            self.mut_pdb_res = None    # The pdb residue to which the transcript mutation point is aligned
            self.mut_pdb_icode = None  # The pdb insertion code of the residue to which the transccript mut point aligned
            self.continuous_to_left = -1  # continuous residues to "left" in pdb.  -1 if our residue is not in the pdb at all
            self.continuous_to_right = -1 # continuous residues to "right" in pdb  -1 if our residue of interest not in pdb
            self.distance_from_aligned = None # minimum of the 2 continuous_to_* values
            self.analyzable = None     # 'Yes', 'No', or 'Maybe' - supports legacy understanding
            self.trans_mut_pos = None  # The 1..len(transcript) transcript position of the mutation site of interest
            self.trans_first_aligned = None # The 1..len(transcript) transcript position of the first aligned amino acid
            self.trans_last_aligned = None  # The 1..len(transcript) transcript position of the last aligned amino acid

        def __str__(self):
            return str(vars(self))

        def __repr__(self):
            return "ChainInfo:" +  str(self)

        def set_alignment_profile(self,alignment: PDBMapAlignment,structure): 
            self.perc_aligned = alignment.perc_aligned
            self.perc_identity = alignment.perc_identity
            self.aln_score = alignment.aln_score

            sorted_align_seqs = sorted(alignment.seq_to_resid.keys())
            self.trans_first_aligned = sorted_align_seqs[0]
            self.trans_last_aligned = sorted_align_seqs[-1]

            chain_aa_letter = None
            self.continuous_to_left = -1
            self.continuous_to_right = -1
            if self.trans_mut_pos: # Typically we have a transcript position nutation number
                if self.trans_first_aligned <= self.trans_mut_pos <= self.trans_last_aligned:
                    self.distance_from_aligned = 0 # I.e., our mutation of interest inside the aligned transcript aa_seq range
                elif self.trans_mut_pos < self.trans_first_aligned: # Our mutation of interest to 'left' of first aligned transcript sequence
                    self.distance_from_aligned = self.trans_first_aligned - self.trans_mut_pos 
                else: # Our mutation of interest to 'right' of last aligned transcript aa_seq
                    self.distance_from_aligned = self.trans_mut_pos - self.trans_last_aligned
                if self.distance_from_aligned == 0:
                    self.analyzable = 'Yes'
                elif self.distance_from_aligned <= args.maybe:
                    self.analyzable = 'Maybe'
                else:
                    self.analyzable = 'No'

                resid = alignment.seq_to_resid.get(self.trans_mut_pos,None)
                if resid:
                    self.continuous_to_left = 0
                    self.continuous_to_right = 0
                    self.mut_pdb_res = resid[1]
                    self.mut_pdb_icode = resid[2]
                    chain_aa_letter = seq1(structure[0][self.chain_id][resid].get_resname().lower())
                    while True:
                        test_left = self.trans_mut_pos - self.continuous_to_left - 1
                        if (test_left <= 0) or (test_left not in alignment.seq_to_resid):
                            break
                        self.continuous_to_left += 1
                    while True:
                        test_right = self.trans_mut_pos + self.continuous_to_right + 1
                        if (test_right > transcript.len) or (test_right not in alignment.seq_to_resid):
                            break
                        self.continuous_to_right += 1

        def mers_string(self,chain_count: int) -> str:
            return { # dimer/trimer... or 7-mer/9-mer via .get default
                1: 'mono',
                2: 'di',
                3: 'tri',
                4: 'tetra',
                5: 'penta',
                6: 'hexa',
                8: 'octa',
               10: 'deca'}.get(chain_count,"%d-"%chain_count) + 'mer'

    # Now that we have all the candidate pdbs and biounits in memory.....
    # Create the dataframe that we will use to pick the 'best' pdb(s) for pipeline analysis
    df_columns = vars(ChainInfo()).keys()
    ci_pdbs_df = pd.DataFrame(columns=df_columns).astype({'biounit':bool})
    ci_modbase_swiss_df = ci_pdbs_df.copy(deep = True)
    ci_usermodel_df = ci_pdbs_df.copy(deep = True)

    ################################################
    # Load relevant swiss models for consideration # 
    ################################################
    swiss_modelids = PDBMapSwiss.unp2modelids(unp,IsoformOnly = True)
    for swiss_modelid in swiss_modelids:
        info = PDBMapSwiss.get_info(swiss_modelid)
        ci = ChainInfo('swiss')
        ci.struct_filename= PDBMapSwiss.get_coord_file(swiss_modelid)
        
        if not os.path.exists(ci.struct_filename):
            LOGGER.critical("Swiss model %s in meta data not found at %s",swiss_modelid,ci.struct_filename)
            continue
        ci.structure_id = info['modelid']
        ci.structure_url = "https://swissmodel.expasy.org/repository/%s/report"%info['coordinate_id']
        ci.biounit = False
        ci.pdb_template = info['template']
        # Swiss templates always formatted as pdb_id.1.chain_id
        ci.template_chain = ci.pdb_template[-1]

        if mutation:
            ci.trans_mut_pos = int(mutation[1:-1])
        else:
            ci.trans_mut_pos = None

        remark3_metrics = PDBMapSwiss.load_REMARK3_metrics(swiss_modelid)
        ci.method = remark3_metrics['mthd']

        LOGGER.info("Loading Swiss model from %s"%ci.struct_filename)
        swiss_structure = PDBParser().get_structure(ci.structure_id,ci.struct_filename)
        chain_count = 0
        for chain in swiss_structure[0]:
            chain_count += 1

        first_swiss_chain_id = next(swiss_structure[0].get_chains()).get_id()

        # We have the task of figuring out the best chain ID to run ddG on - typically
        # we run best on the template chain  in the structure name... but only bother 
        # to try to figure this out if we have a multimeric swiss model

        if chain_count == 1: 
            ci.chain_id = first_swiss_chain_id # Easy when swiss model is a monomer
        elif 'smtle' in remark3_metrics: # Grab the chain ID from structure name
            ci.chain_id = remark3_metrics['smtle'][-1]
        elif 'chain' in remark3_metrics:
            ci.chain_id = remark3_metrics['chain']
            LOGGER.warning("No chain ID was found in the REMARK 3 header of %s",ci.struct_filename)
        else: # Strange - just grab first 
            ci.chain_id = first_swiss_chain_id
        assert ci.chain_id not in [' ',''],"UUGH - this swiss model has blank chain ID"

        ci.nresidues = len(list(swiss_structure[0][ci.chain_id].get_residues()))
        ci.biounit_chains = 1 # This _could_ hopefully change - but not for now
        ci.mers = ci.mers_string(ci.biounit_chains)
        alignment = PDBMapAlignment(swiss_structure,ci.chain_id,None,'trivial')
        success,errormsg = alignment.align_trivial(transcript,swiss_structure,ci.chain_id)
        assert success,errormsg
        ci.method = 'swiss'
        ci.set_alignment_profile(alignment,swiss_structure)

        ci_modbase_swiss_df = ci_modbase_swiss_df.append(vars(ci),ignore_index=True)
        if chain_count > 1:
            LOGGER.info("Adding swiss model as %d-mer"%chain_count) 
            ci.biounit_chains = chain_count
            ci.mers = ci.mers_string(ci.biounit_chains)
            ci_modbase_swiss_df = ci_modbase_swiss_df.append(vars(ci),ignore_index=True)

    ################################################
    # Load relevant modbase models for consideration # 
    ################################################
    def chain_info_from_modbase_summary(modbase,modbase_summary) -> ChainInfo:
        """Because we may want to support not only modbase current models, but also older ones
        It is desirable to centralize processing of the models for all 3 class types"""
   
        ci = ChainInfo('modbase')
        modbase_modelid = modbase_summary['database_id'] # This is the ENSP.... model filename
        ci.struct_filename= modbase.get_coord_file(modbase_modelid)

        if not os.path.exists(ci.struct_filename):
            LOGGER.critical("Modbase model %s in meta data not found at %s",modbase_modelid,ci.struct_filename)
            return None
        ci.structure_id = modbase_modelid
        # URL per Alex at Sali lab - thanks!
        # This method (column model_id) not available for modbase 2013
        if "model_id" in modbase_summary:
            ci.structure_url = "https://salilab.org/modbase/searchbyid?modelID=%s&displaymode=moddetail"%modbase_summary['model_id']
        ci.template_identity = float(modbase_summary['sequence identity'])

        ci.biounit = False
        ci.pdb_template = modbase_summary['pdb_code']

        if mutation:
            ci.trans_mut_pos = int(mutation[1:-1])
        else:
            ci.trans_mut_pos = None

        ci.method = None

        LOGGER.info("Opening Modbase Model File %s"%ci.struct_filename)
        if ci.struct_filename.endswith(".gz"):
            modbase_structure_fin = gzip.open(ci.struct_filename,'rt')
        elif ci.struct_filename.endswith(".xz"):
            modbase_structure_fin = lzma.open(ci.struct_filename,'rt')
        elif ci.struct_filename.endswith(".pdb"):
            modbase_structure_fin = open(ci.struct_filename,'rt')
        else:
            sys.exit("Modbase model has unrecognized file extension %s"%ci.struct_filename)
        modbase_structure = PDBParser().get_structure(ci.structure_id,modbase_structure_fin)
        modbase_structure_fin.close()
        modbase_structure_fin = None
        # Modbase models do not enclode chain IDs to match template - so just grab what we have
        # I.e. get_chains() is an iterator, and we want to first element that comes back
        # without setting up a for loop - so next() is perfect
        modbase_chainid = next(modbase_structure[0].get_chains()).get_id()
        if modbase_chainid in [' ','']:
            modbase_structure[0][modbase_chainid].id= 'A'
            ci.chain_id = 'A'
        else:
            ci.chain_id = modbase_chainid
        
        # Modbase super-annoying because the chain_id is often missing
        ci.nresidues = len(list(modbase_structure[0][ci.chain_id].get_residues()))
        ci.biounit_chains = 1 # This _could_ hopefully change - but not for now
        ci.mers = ci.mers_string(ci.biounit_chains)
        alignment = PDBMapAlignment(modbase_structure,ci.chain_id,None,'trivial')
        success,errormsg = alignment.align_trivial(transcript,modbase_structure,ci.chain_id)
        if not success:
            LOGGER.critical("Modbase Model %s fails to align to transcript - skipping.  Error:\n%s",
                ci.structure_id,errormsg)
            return None
        ci.method = 'modbase'
        ci.set_alignment_profile(alignment,modbase_structure)
        return ci

    modbase_2020 = PDBMapModbase2020(config_dict)
    modbase_summary_rows = modbase_2020.transcript2summary_rows(transcript)
    for modbase_summary in modbase_summary_rows:
        ci = chain_info_from_modbase_summary(modbase_2020,modbase_summary)
        if ci is not None:
            ci_modbase_swiss_df = ci_modbase_swiss_df.append(vars(ci),ignore_index=True)
        
    ################################################################
    # Load relevant Modbase2016  AND Modbase2013 for consideration # 
    ################################################################

    # ensp_proteins = PDBMapProtein.unp2ensp(unp)
    # enst_transcripts = PDBMapProtein.unp2enst(unp)
    # if ensp_proteins:
    #    for enst_transcript in [PDBMapTranscriptEnsembl(enst_id) for enst_id in enst_transcripts]:
    #        if enst_transcript.aa_seq == transcript.aa_seq:
    #            break
    #        LOGGER.critical("Ensembl transcript for %s does not match Uniprot transcript sequence for %s:\nENST:%s\nUniprot:%s",
    #             enst_transcript.id,transcript.id,enst_transcript.aa_seq,transcript.aa_seq)
    #        LOGGER.critical("No Ensembl models will be considered for %s",enst_transcript.id)
    #    else:
    #        enst_transcript = None # We looped through and no transcript had aa_seq like
    #    # for modbase_class in [PDBMapModbase16, PDBMapModbase13]
    #    for modbase in [PDBMapModbase2016(config_dict), PDBMapModbase2013(config_dict)] if enst_transcript else []:
    #        for ensp_protein in ensp_proteins:
    #            modbase_modelids = modbase.ensp2modelids(ensp_protein)
    #            for modbase_modelid in modbase_modelids:
    #                info = modbase.get_info(modbase_modelid)
    #                ci = chain_info_from_modbase_info(info)
    #                if ci is not None:
    #                    ci_modbase_swiss_df = ci_modbase_swiss_df.append(vars(ci),ignore_index=True)
        
    LOGGER.info("%d modbase swiss identified."%len(ci_modbase_swiss_df))

    ###########################################################################################################################

    full_complexes = set() # A set of pdbids for which we have noted a full complex is to be considered
    for pdb_upper_id,_chain_id in pdbs_chains:
        LOGGER.info("Processing the next tuple: %s %s",pdb_upper_id,_chain_id)
        ci = ChainInfo('pdb')
        # Extract the mutation position
        if mutation:
            ci.trans_mut_pos = int(mutation[1:-1])
        else:
            ci.trans_mut_pos = None
        ci.structure_id = pdb_upper_id.lower()
        ci.chain_id = _chain_id
        if ci.structure_id not in structures_dict:
            LOGGER.warning("Not considering SIFTS suggestion %s.%s as missing from filesystem"%(ci.structure_id,ci.chain_id))
            continue
        ci.structure_url = "https://www.rcsb.org/structure/%s"%ci.structure_id.upper()
        ci.biounit = (ci.structure_id in biounits_dict)
        ci.biounit_chains = 1
        structure = structures_dict[ci.structure_id]

        # UUGH - trying to duplicate the old REMARK 2 RESOLUTION entry in the old pdb format - still in process
        # At first Ezra.Peisach@rcsb.org wrote that the replacement is _reflns.d_resolution_high.
        # But, when I hit a wall of that entry being missing in some structures, he said sorry
        # and that _refine.ls_d_res_high is the right replacement
        mmCIF_dict = mmCIFDicts_dict[ci.structure_id]
        ci.method = mmCIF_dict['_exptl.method'][0]
        raw_deposition_date = mmCIF_dict.get('_pdbx_database_status.recvd_initial_deposition_date',None)
        if raw_deposition_date: # Split into YYYY-MM-DD
            raw_YYYY_MM_DD = raw_deposition_date[0].split('-')
            if len(raw_YYYY_MM_DD) == 3:
                ci.deposition_date = datetime.datetime(
                    int(raw_YYYY_MM_DD[0]),int(raw_YYYY_MM_DD[1]),int(raw_YYYY_MM_DD[2]))

            if not ci.deposition_date:
                LOGGER.warning("Deposition date %s apparently not YYYY-MM-DD",raw_deposition_date)

        resolution_keys = []
        ci.resolution = None # Not applicable is correct for NMR structures
        if ci.method == 'X-RAY DIFFRACTION':
            resolution_keys = ['_refine.ls_d_res_high','_reflns.d_resolution_high','_refine_hist.d_res_high']
            ci.method = 'X-RAY'
        elif ci.method == 'ELECTRON CRYSTALLOGRAPHY':
            ci.method = 'ELECTRON CRYSTALLOGRAPHY'
            resolution_keys = ['_em_3d_reconstruction.resolution']
        elif ci.method == 'ELECTRON MICROSCOPY':
            ci.method = 'EM'
            resolution_keys = ['_em_3d_reconstruction.resolution']
        elif ci.method.find('NMR') > -1: # resolution is not applicable for nmr
            resolution_keys = []
        elif ci.method.find('SOLUTION SCATTERING') > -1: # resolution is not applicable for SOLUTION SCATTERING
            LOGGER.warning('SOLUTION SCATTERING may be of low resolution: %s',ci.structure_id)
            resolution_keys = []
        else:
            resolution_keys = ['_refine.ls_d_res_high','_reflns.d_resolution_high','_refine_hist.d_res_high']
            LOGGER.critical('%s Experimental method for %s is unknown'%(ci.method,ci.structure_id))

        for resolution_key in resolution_keys:
            if resolution_key in mmCIF_dict:
                # Congratulations.  You found the old REMARK 2 RESOLUTION entry.... or so we hope.
                # HOWEVER, sometimes this key is a '.' or '?' and we should skip that
                resolution_value = mmCIF_dict[resolution_key]
                if resolution_value != '.' and resolution_value != '?':
                    ci.resolution = float(mmCIF_dict[resolution_key][0])
                break
        if not ci.resolution and ci.method.find('X-RAY') > 0:
            LOGGER.warning("pdb %s is method=%s with no resolution entry"%(ci.structure_id,method))

        if ci.biounit: # In case of biounit we have good chance of a multi-chain complex to process
            chain_ATOM_residues = {}
            for chain in biounits_dict[ci.structure_id][0]:
               for residue in chain:
                   if 0 == len(residue.id[0].strip()): # Then it is an ^ATOM entry (not HETERO or WAT)
                       count = chain_ATOM_residues.get(chain.id,0)
                       chain_ATOM_residues[chain.id] = count+1
            ci.biounit_chains = len(chain_ATOM_residues)
            ci.nresidues = sum(chain_ATOM_residues.values())
        else:
            ci.nresidues = len(list(structure[0][ci.chain_id].get_residues()))
            
        alignment = PDBMapAlignment(structure,ci.chain_id,None,'trivial')
        assert ci.chain_id in structure[0],"Chain ID %s not found in file %s"%(ci.chain_id,structure_filename)
        success = False
        success,errormsg = alignment.align_trivial(transcript,structure,ci.chain_id)
        if not success:
            LOGGER.debug("Trivial alignment failed - which is common: %s",errormsg)
            if unp_is_canonical:
                success,errormsg = alignment.align_sifts_canonical(transcript,structure,ci.chain_id)
            if not success:
                pdb_seq_resid_xref = PDBMapAlignment.create_pdb_seq_resid_xref(mmCIF_dict)
                success,errormsg = alignment.align_sifts_isoform_specific(transcript,structure,pdb_seq_resid_xref,chain_id=ci.chain_id,is_canonical=unp_is_canonical)
    
            if not success:
                LOGGER.warning("Unable to align %s.%s with sifts: %s",ci.structure_id,ci.chain_id,errormsg)
                LOGGER.warning("Attempting to align chain %s with biopython call",ci.chain_id)
                success,errormsg = alignment.align_biopython(transcript,structure,ci.chain_id)
                if not success:
                    LOGGER.critical("Also Unable to align %s.%s with biopython: %s",ci.structure_id,ci.chain_id,errormsg)

        if not success:
            LOGGER.critical("Skipping %s:%s because of alignment failure"%(ci.structure_id,ci.chain_id))
            continue

        ci.set_alignment_profile(alignment,structure)


        ci_pdbs_df = ci_pdbs_df.append(vars(ci),ignore_index=True)

        # If this is a multi-chain structure, add a row to process as complex
        if ci.biounit_chains > 1 and ci.structure_id not in full_complexes:
            full_complexes.add(ci.structure_id)
            ci.mers = ci.mers_string(ci.biounit_chains)
            full_complex = True
            ci_pdbs_df = ci_pdbs_df.append(vars(ci),ignore_index=True)
           
    # ci_pdbs_df.to_csv('/tmp/ci_pdbs_df.tsv',sep='\t',index=True)
    #    # alignments_dict[(pdb_id,ci.chain_id)] = alignment
    
    LOGGER.info("Selecting 'best' chains from %d unique pdb chains based on length, coverage, alignments"%len(ci_pdbs_df))

    # LOGGER.info("Calling get_structures() to load PDB/ModBase/Swiss structures of %s %s unp=%s from MySQL..."%(gene,refseq,unp))
    # LOGGER.info("Identifying PDB/ModBase/Swiss structures of %s %s unp=%s..."%(gene,refseq,unp))
    # ci_pdbs_df = get_pdbs(io,unp,ci.trans_mut_pos,maybe_threshold = args.maybe)
    #LOGGER.info("%d pdbs identified."%len(ci_pdbs_df))

   


    # Funny bug - but if ci_pdbs_df is empty, the pd.concat turns all the second dF_modbase_swiss columns from ints to floats
    if len(ci_pdbs_df) == 0:
        ci_df = ci_modbase_swiss_df.copy()
    elif len(ci_modbase_swiss_df) == 0:
        ci_df = ci_pdbs_df.copy()
    else:
        ci_df = pd.concat([ci_pdbs_df,ci_modbase_swiss_df],ignore_index = True)

    # df = get_structures(io,unp,ci.trans_mut_pos,maybe_threshold = args.maybe)
    # LOGGER.info("%d structures identified."%len(df))
    
    # Use shorter header description for tighter htm output
    # 2020 Feb - Let the report module LATER worry about renaming columns.  Otherwise, I lose my mind
    # df.rename(columns={'Analysis Possible?': 'Analyzable?','Sequence Start': 'Seq Start','Sequence End': 'Seq End','Sequence Identity': 'template_identity'}, inplace=True)
    
    # Add in Swiss-model quality metrics from the REMARK 3 entries of the model file
    # This probably should be moved off to tables from load-time
    # LOGGER.warning("DF columns = %s"%str(df.columns))
    for i in range(len(ci_df)):
        ci_row = ci_df.iloc[i]
        if ci_row['label'] == 'swiss':
            remark3_metrics = PDBMapSwiss.load_REMARK3_metrics(ci_row['structure_id'])
            ci_df.iat[i,ci_df.columns.get_loc('template_identity')] = float(remark3_metrics['sid'])
            
    if not ci_df.empty: #ChrisMoth added this line and indented as the code below that makes no sense for empty df
        # df["Gene"] = gene
        # df = df[["gene","label","structure_id","Chain","method","Resolution (PDB)","template_identity",
        #                "PDB Template","Residues","Seq Start","Seq End",
        #                "Transcript Pos","PDB Pos","Distance to Boundary","Analyzable?"]]
    
        df_dropped = pd.DataFrame(columns=ci_df.columns.values.tolist() + ["drop_reason"])
    
        # We now take pains to rid our candidate structures of models which are clearly redundant/degenerate
        # The algorithm retains all experimental structures
        def drop_df_row(k,reason,keeping = -1):
            global df_dropped
            if keeping >= 0: # Add the other structure ID being kept, to the reason text, if asked 
                reason += ' vs %s.%s'%(ci_df.iloc[keeping]['structure_id'],ci_df.iloc[keeping]['chain_id'])
            if ci_df.iloc[k]['structure_id'].upper() in keep_pdbs:
                LOGGER.warning( "Retaining %s.%s because it is the keep list, even though:\n%s",
                    ci_df.iloc[k]['structure_id'],ci_df.iloc[k]['chain_id'],reason)
                return False
            LOGGER.info( "Dropping ci_df duplicate row %d because %s %s"%(k,ci_df.iloc[k]['structure_id'],reason))
 
            # Record the removed model for the drop report
            dropped_model = ci_df.iloc[k].copy()
    
            dropped_model['drop_reason'] = reason 
            df_dropped = df_dropped.append(dropped_model)
    
            # Remove model from set to process
            ci_df.drop(ci_df.index[k],inplace=True) 
            return True

        # Jonathan requested that all "No" structures be dropped - so let's do that first 
        i = 0
        while i < ci_df.shape[0]:
            row_i = ci_df.iloc[i]
            if row_i['analyzable'] == 'No':
                if not drop_df_row(i,"Not Analyzable"):
                    i += 1
            else:
                i += 1

        i = 0
        j = 0
        while i < ci_df.shape[0]-1:
            row_i = ci_df.iloc[i]
    
            row_i_template = None
            if row_i['pdb_template']: # This is an experimenal structure
                row_i_template = row_i['pdb_template'].split('.')[0].upper() 
    
            j = i+1
            row_i_aligned_len = 1+int(row_i['trans_last_aligned'] or 0) - int(row_i['trans_first_aligned'] or row_i['trans_last_aligned'] or 0)
            while j < ci_df.shape[0]:
                row_j = ci_df.iloc[j]
                row_j_aligned_len = 1+int(row_j['trans_last_aligned'] or 0) - int(row_j['trans_first_aligned'] or row_j['trans_last_aligned'] or 0)

                # Do not discard a polymer for a monomer or vice versa.  Keep both
                if row_i['mers'] != row_j['mers']:
                    pass
                # If both are experimental structures, then unique care is warranted before delete
                elif row_i['label'] == 'pdb' and row_j['label'] == 'pdb':
                    # If we have a "yes",then drop the maybes and nos
                    # If we have a "maybe",then drop the nos
                    if row_i['analyzable'] != row_j['analyzable']:
                        if row_i['analyzable'] == 'No':
                            if drop_df_row(i,"is not as analyzable",j):
                                break
                        elif row_j['analyzable'] == 'No':
                            if drop_df_row(j,"is not as analyzable",i):
                                break
                        elif row_i['analyzable'] == 'Maybe':
                            if drop_df_row(i,"is not as analyzable",j):
                                break
                        elif row_j['analyzable'] == 'Maybe':
                            if drop_df_row(j,"is not as analyzable",i):
                                break
                    else: # Both experimental structures _are_ (or are not) analyzable - same level...
                        if row_i['analyzable'] == 'Yes': # Generally this is good - but some peptides appear 100x in pdb and we do not need all that
                            if row_i['method'] != row_j['method']: # Don't throw away an xray for a cryo-em or nmr
                                pass
                            elif (row_i['trans_first_aligned'] == row_j['trans_first_aligned']) and (row_i['trans_last_aligned'] == row_j['trans_last_aligned']):
                                # Both penptides cover same region.. and this is the 
                                # same method experiment.  So, Delete lower resolution
                                # When comparing 2 NMR structures, these wll both be None
                                row_i_resolution = row_i['resolution']
                                row_j_resolution = row_j['resolution']
                                if row_j_resolution == row_i_resolution: # Two structures with same resolution - let one go
                                    if drop_df_row(i,"is same resolution",j):
                                        break
                                elif float(row_i_resolution) > float(row_j_resolution): # Then keep jth
                                    if drop_df_row(i,"is lower resolution",j):
                                        break
                                else:
                                    # row_j_resolution > row_i_resolution: # Then keep ith
                                    if drop_df_row(j,"is lower resolution",i):
                                        break
                            elif ((row_i_aligned_len > row_j_aligned_len) 
                                         and  (row_i['trans_first_aligned'] <= row_j['trans_first_aligned']) 
                                         and (row_i['trans_last_aligned'] >= row_j['trans_first_aligned'])):
                                if drop_df_row(j,"is shorter",i):
                                    break
                            elif ((row_j_aligned_len > row_i_aligned_len) 
                                         and (row_j['trans_first_aligned'] <= row_i['trans_first_aligned']) 
                                         and (row_j['trans_last_aligned'] >= row_i['trans_first_aligned'])):
                                if drop_df_row(i,"is shorter",j):
                                    break
                            # Otherwise, the structures might overlap - but these are experimental structures - do not throw away! 
    
                        elif row_i['analyzable'] == 'No' or row_i['analyzable'] == 'Maybe': # Neither is very analyzable
                            if int(row_i['distance_from_aligned']) > int(row_j['distance_from_aligned']):
                                if drop_df_row(i,"is farther from mutation",j):
                                    break
                            elif int(row_i['distance_from_aligned']) < int(row_j['distance_from_aligned']):
                                if drop_df_row(j,"is farther from mutation",i):
                                    break
                            if row_i_aligned_len > row_j_aligned_len:
                                if drop_df_row(j,"is shorter",i):
                                    break
                            elif row_j_aligned_len > row_i_aligned_len:
                                if drop_df_row(i,"is shorter",j):
                                    break
                            else: # Both 'Np' structures have same length.   Eliminate by resolution difference
                                # CAREFUL.  In this section, a resolution _could_ be None
                                row_i_resolution = row_i['resolution']
                                row_j_resolution = row_j['resolution']
                                if row_i_resolution == row_j_resolution: # Same or both NMR (None)
                                    row_i_structid = row_i['structure_id']
                                    row_j_structid = row_j['structure_id']
                                    if (row_i_structid > row_j_structid):
                                        if drop_df_row(i,"is same length and resolution, but greater structid ID",j):
                                            break
                                    else:
                                        if drop_df_row(j,"is same length and resolution, but greater than or same structid ID",i):
                                            break
                                elif row_i_resolution and (float(row_i_resolution) > float(row_j_resolution)):
                                    if drop_df_row(i,"is higher resolution",j):
                                        break
                                else:
                                    # row_j_resolution > row_i_resolution:
                                    if drop_df_row(j,"is higher resolution",i):
                                        break
                        else: # Both structures are "yes" to analyzable - and much greater caution before dropping one
                            pass
                # Don't remove anything if comparing an experimental structure to a model
                elif row_i['label'] == 'pdb' or row_j['label'] == 'pdb':
                    pass
                else: # We are going to compare two models (not comparing experimental structures)
                    row_j_template = row_j['pdb_template'].split('.')[0].upper()
                    if (row_i_template == row_j_template): # If both models use same template, we need to drop one of them
                        # If both are ENSP models, then we need to just get rid of the 2013 one, and keep the 2016 one containing a period
                        if row_i['label'] == 'modbase' and row_j['label'] == 'modbase':
                            if (('.' in row_i['structure_id'] and '.' in row_j['structure_id']) or # if Both structures modbase16
                                    ('.' not in row_i['structure_id'] and '.' not in row_j['structure_id'])): # or if Both structures modbase13
                                if row_i_aligned_len > row_j_aligned_len:
                                    if drop_df_row(j,"is shorter modbase",i):
                                        break
                                elif row_j_aligned_len > row_i_aligned_len:
                                    if drop_df_row(i,"is shorter modbase",j):
                                        break
                                else: # The two modbase models are on the same structure - but must be two different chains
                                    # The Modbase structure_id has the chain after the final . in the name
                                    row_i_chain = row_i['structure_id'].split('.')[-1].upper()
                                    row_j_chain = row_j['structure_id'].split('.')[-1].upper()
                                    if (row_i_chain > row_j_chain):
                                        if drop_df_row(i,"is same length, but greater chain ID",j):
                                            break
                                    else:
                                        if drop_df_row(j,"is same length, but greater than or same chain ID",i):
                                            break
                            elif '.' in row_i['structure_id']: # Then row_i is Modbase2016 - so drop row_j
                                if drop_df_row(j,"is modbase13",i):
                                    break
                            elif '.' in row_j['structure_id']: # Then row_j is Modbase2016 - so drop row_i
                                if drop_df_row(i,"is modbase13",j):
                                    break
                        # If one of the same-template models is swiss, other modebase
                        # Kepp swiss
                        elif row_i['label'] == 'swiss' and row_j['label'] == 'modbase':
                            if drop_df_row(j,"is modbase (keeping swiss)",i):
                                break
                        elif row_i['label'] == 'modbase' and row_j['label'] == 'swiss':
                            if drop_df_row(i,"is modbase (keeping swiss)",j):
                                break
                        # If they are both swiss, then keep the longer chain... or lower chain ID if same length chains
                        elif row_i['label'] == 'swiss' and row_j['label'] == 'swiss':
                            if row_i_aligned_len > row_j_aligned_len:
                                if drop_df_row(j,"is shorter swissmodel",i):
                                    break
                            elif row_j_aligned_len > row_i_aligned_len:
                                if drop_df_row(i,"is shorter swissmodel",j):
                                    break
                            else: # The two swiss models are on the same structure - but must be two different chains
                                # The Swiss structure_id has the chain after the final . in the name
                                row_i_chain = row_i['structure_id'].split('.')[-1].upper()
                                row_j_chain = row_j['structure_id'].split('.')[-1].upper()
                                if (row_i_chain > row_j_chain):
                                    if drop_df_row(i,"is same length, but greater chain ID",j):
                                        break
                                else:
                                    if drop_df_row(j,"is same length, but greater than or same chain ID",i):
                                        break
                # end of dup drops in models
    
                j += 1
            else: # Applies to outer while.  Do not increment i if a row was dropped (and break was called)
                i += 1
    
        # Important: This function drops by index, not the physical row number as the duplicate drop code above
        def drop_df_low_identity_rows(droplist,reason):
            global df_dropped
            if len(droplist) < 1:
                return
            for k in droplist:
                LOGGER.info("Dropping of row %d because %s %s"%(k,ci_df.loc[k]['structure_id'],reason))
            dropped_models = ci_df.loc[droplist].copy()
            dropped_models['drop_reason'] = reason 
            df_dropped = df_dropped.append(dropped_models)
            ci_df.drop(droplist,inplace=True) 
    
        models_sub10_seq_identity = ci_df.index[(((ci_df['label'] == 'modbase') | (ci_df['label'] == 'swiss')) & (ci_df['template_identity'] < 10.0))]
        if len(models_sub10_seq_identity) > 0:
            drop_df_low_identity_rows(models_sub10_seq_identity,"template_identity < 10.0 unacceptable")
        models_10to20_seq_identity = ci_df.index[(((ci_df['label'] == 'modbase') | (ci_df['label'] == 'swiss')) & ((10.0 <= ci_df['template_identity']) & (ci_df['template_identity'] < 20.0)))]
        models_20to30_seq_identity = ci_df.index[(((ci_df['label'] == 'modbase') | (ci_df['label'] == 'swiss')) & ((20.0 <= ci_df['template_identity']) & (ci_df['template_identity'] < 30.0)))]
        models_over30_seq_identity = ci_df.index[(((ci_df['label'] == 'modbase') | (ci_df['label'] == 'swiss')) & (30.0 <= ci_df['template_identity']))]
    
        # Do we have a solid group of anlyzable (mutation covering models?
        analyzable_over_30 = ci_df.index[(ci_df["analyzable"] == "Yes") & ((ci_df['label'] == 'pdb') | (((ci_df['label'] == 'modbase') | (ci_df['label'] == 'swiss')) & (30.0 <= ci_df['template_identity'])))]
        analyzable_under_30 = ci_df.index[(ci_df["analyzable"] == "Yes") & (((ci_df['label'] == 'modbase') | (ci_df['label'] == 'swiss')) & (ci_df['template_identity'] < 30.0))]
    
        # If we have 3 pretty good models, throw away the low identity ones
        # Unless we are hurting badly for analyzable models, and we have some
        if len(models_over30_seq_identity) > 3:
            if len(analyzable_over_30) == 0 and len(analyzable_under_30) > 0: # In this case (no 'good' structures/models with mutation coverage)
                models_not_analyzable_under_30 = ci_df.index[(((ci_df['label'] == 'modbase') | (ci_df['label'] == 'swiss')) & (ci_df['analyzable'] != 'Yes') & (ci_df['template_identity'] < 30.0))]
                drop_df_low_identity_rows(models_not_analyzable_under_30,"Retaining only low template_identity Models with mutation coverage")
            else:
                drop_df_low_identity_rows(models_10to20_seq_identity,"Low template_identity Models are not needed")
                drop_df_low_identity_rows(models_20to30_seq_identity,"Low template_identity Models are not needed")
    
        elif len(models_20to30_seq_identity) > 3:
            if len(analyzable_over_30) == 0 and len(analyzable_under_30) > 0: # In this case (no 'good' structures/models with mutation coverage)
                models_not_analyzable_under_20 = ci_df.index[(((ci_df['label'] == 'modbase') | (ci_df['label'] == 'swiss')) & (ci_df['analyzable'] != 'Yes') & (ci_df['template_identity'] < 30.0))]
                drop_df_low_identity_rows(models_not_analyzable_under_20,"Retaining only low template_identity Models with mutation coverage")
            else:
                drop_df_low_identity_rows(models_10to20_seq_identity,"Low template_identity Models are not needed")
    
        # Do we have a solid group of anlyzable (mutation covering models?
        # If we have analyzable models over 20 at all - then drop all models under_20
        analyzable_over_20 = ci_df.index[(ci_df["analyzable"] == "Yes") & ((ci_df['label'] == 'pdb') | (((ci_df['label'] == 'modbase') | (ci_df['label'] == 'swiss')) & (20.0 <= ci_df['template_identity'])))]
        if len(analyzable_over_20) > 0:
            analyzable_under_20 = ci_df.index[(ci_df["analyzable"] == "Yes") & (((ci_df['label'] == 'modbase') | (ci_df['label'] == 'swiss')) & (ci_df['template_identity'] < 20.0))]
            if len(analyzable_under_20) > 0:
                drop_df_low_identity_rows(analyzable_under_20,"Low template_identity covering models not needed")
    
        ci_df = ci_df.sort_values(by=["distance_from_aligned","method","nresidues",
                                              "resolution","template_identity"],
                                              ascending=[True,False,False,True,False])
    
        dropped_models_filename = "%s/%s_%s_%s_dropped_models.csv"%(mutation_dir,gene,refseq,mutation)
        LOGGER.info("Output the %d dropped  ModBase and Swiss models to: %s"%(len(df_dropped),dropped_models_filename))
        df_dropped.to_csv(dropped_models_filename,sep='\t',index=True)
        #
        # End if not df_empty

    ENST_transcript_ids = PDBMapProtein.unp2enst(unp)
    if ENST_transcript_ids:
        LOGGER.info("ENST transcript for %s will be %s",unp,ENST_transcript_ids[0])
        if len(ENST_transcript_ids) > 1:
            LOGGER.warning("%d additional ENST transcripts will be ignored",len(ENST_transcript_ids)-1)

    elif (config_pathprox_dict['disease1_variant_sql_label'] == 'tcga') or \
             (config_pathprox_dict['disease2_variant_sql_label'] == 'tcga'):
        LOGGER.error("The uniprot identifier must map to an ENST transcript ID in order use tcga variants",user_model,model_count)
        sys.exit(1)

    if user_model: # Then add it to our structure report
        if not ENST_transcript_ids:
            LOGGER.error("The uniprot identifier must map to an ENST transcript ID in order to incorporate a user model",user_model,model_count)
            sys.exit(1)

        from Bio.PDB.PDBParser import PDBParser 
        p = PDBParser()
        user_model_split = user_model.split(':')
        if len(user_model_split) > 2:
            LOGGER.error("I don't understand user model %s.  Use only one colon to delimit filename from chain please",user_model)
            sys.exit(1)

        ci = ChainInfo('usermodel')
        # If chain ID specified in missense.csv
        if len(user_model_split) == 2:
            ci.chain_id = user_model_split[1]

        ci.struct_filename = user_model_split[0]
        ci.structure_id = ci.struct_filename.split('.')[0] # Drop the .mmcif or .pdb from user model

        structure = p.get_structure('user_model',ci.struct_filename)

        model_count = 0
        chain_count = 0;

        for model in structure:
            model_count += 1
            if model_count > 1:
                LOGGER.critical("User model %s must contain one and only one model, not %d",user_model,model_count)
                sys.exit(1)
            for chain in model:
                # If no chain ID specified in missense.csv
                if len(user_model_split) == 1:
                    ci.chain_id = chain.get_id()
                chain_count += 1
                # If multi-chain model, and chain NOT specified by user
                # then we cannot continue
                if chain_count > 1 and len(user_model_split) == 1: 
                    LOGGER.critical("User_model %s has multiple chains.  You must specifiy a :ChainID for pathprox",user_model)
                    sys.exit(1)

        ci.mers = ci.mers_string(chain_count)

        ci.trans_mut_pos = int(mutation[1:-1])
        if ci.chain_id not in structure[0]:
            LOGGER.critical("Chain ID %s specified is not in user model %s",ci.chain_id.ci.structure_id)
            sys.exit(1)

        if (' ',ci.trans_mut_pos,' ') in structure[0][ci.chain_id]:
             ci.mut_pdb_res = ci.trans_mut_pos
             ci.mut_pdb_icode = ' ' 
        else:
             LOGGER.warning("Mutation POS %d is not in user model %s",ci.trans_mut_pos,ci.structure_id)


        """seq_start = -1 
        seq_end = -1
        for residue in structure[0][ci.chain_id]:
            resid = residue.id
            assert((resid[0] == None or resid[0] == ' ') and ('User model residues cannot be WAT or HET entries'))
            assert((resid[2] == None or resid[2] == ' ') and ('User model insertion codes must be blank'))
            if seq_start == -1:
                seq_start = residue.id[1]
            seq_end = residue.id[1]"""

        # ci_usermodel_df = ci_usermodel_df.append(vars(ci),ignore_index=True)
        alignment = PDBMapAlignment(structure,ci.chain_id,None,'trivial')
        assert ci.chain_id in structure[0],"Chain ID %s not found in file %s"%(ci.chain_id,structure_filename)
        success = False
        success,errormsg = alignment.align_trivial(transcript,structure,ci.chain_id)
        if not success:
            LOGGER.critical("Trivial alignment failed - which is bad for user model",errormsg)
    
            success,errormsg = alignment.align_biopython(transcript,structure,ci.chain_id)
            if not success:
                LOGGER.critical("Also Unable to align with biopython: %s",errormsg)

        if not success:
            LOGGER.critical("Halting on usermodel %s:%s because of alignment failure"%(ci.structure_id,ci.chain_id))
            sys.exit(1)

        ci.method = 'usermodel'
        ci.set_alignment_profile(alignment,structure)

        ci_usermodel_df = ci_usermodel_df.append(vars(ci),ignore_index=True)
      
        # Make the user model to go first in all processing and reports
        ci_df = pd.concat([ci_usermodel_df,ci_df],ignore_index = True).reset_index(drop = True)

 
    # df could be empty or not....
    structure_report_filename = "%s/%s_%s_%s_structure_report.csv"%(mutation_dir,gene,refseq,mutation)
    LOGGER.info("Output the %d retained PDB, ModBase and Swiss structures to: %s"%(len(ci_df),structure_report_filename))
    ci_df.to_csv(structure_report_filename,sep='\t',index=True,na_rep='nan')
    if mutation and len(ci_df) > 0:
        LOGGER.info("Total with coverage of %s: %d"%(mutation,(~ci_df["mut_pdb_res"].isnull()).sum()))
    
    #
    # With the list of structures complete, now create a schedule of work to be later launched
    #
    df_all_jobs = pd.DataFrame()
    
    # This dictionary used to feed patching of .slurm template scripts
    # It now carries info to create a .csv file of work to be performed
    params = {"config":args.config,"userconfig":args.userconfig,"collab":config_dict['collaboration'],"case":args.project,
                        "project":args.project,
                        "pdbid":"N/A",
                        "chain":"N/A",
                        "mers": "N/A",
                        "gene":gene,"refseq":refseq,"unp":unp,'mutation_dir':mutation_dir,
                        "mutation":mutation,
                        "transcript_mutation":mutation,
                        "pdb_mutation":"N/A" # Replace for ddG and Pathprox when we have coverage
             }
 
    if ENST_transcript_ids:
        params["transcript"] = ENST_transcript_ids[0]

    # We run one sequence analysis on each transcript and mutation point.. Get that out of the way
    df_all_jobs = df_all_jobs.append(makejob_udn_sequence(params),ignore_index=True)
    
    # LOGGER.info("%d Structures/models of %s %s (unp=%s)"%(len(df),refseq,gene,unp))
    # If no structures at all are available, then the work plan takes a very different course
    if ci_df.empty:
        if args.mutation:
            LOGGER.info("No structures with coverage of the mutation are available. No report will be generated.")
        else:
            LOGGER.info("No structures are available. No report will be generated.")
    
    
        # If the script made it to this point, set this mutation to "complete"
        # sql  = "UPDATE psb_collab.MutationSummary SET complete=True WHERE "
        # sql += "collab=%s AND project=%s AND gene=%s AND mutation=%s "
        # c = io._con.cursor()
        # c.execute(sql,(config_dict['collaboration'],args.project,gene,mutation))
        # c.close()
    else:
        # Launch ddG on all structures with coverage
        if mutation:
            for index,ci_row in ci_df.iterrows():  # For each chain_info record...
                # Don't run ddg unless this is a simple monomer
                if ci_row['mut_pdb_res'] and ci_row['mers'] == 'monomer':
                    _params = params.copy() # Dont affect the dictionary passed in  Shallow copy fine
                    _params['method'] = ci_row['method']
                    _params['structure_id'] = ci_row['structure_id']
                    _params['struct_filename'] = ci_row['struct_filename'] if ci_row['struct_filename'] else ci_row['structure_id']
                    _params['chain'] = "%s"%ci_row['chain_id']  # Quote it in case the IDentifier is a space, as in many models
                    _params['mers'] = ci_row['mers']
                    _params['label'] = ci_row['label']
                    # Only launchain_id the udn_structure (ddG) script if there is mutation coverage
                    _params['pdb_mutation'] = (
                         mutation[0] + 
                         str(int(ci_row['mut_pdb_res'])) + 
                         str(ci_row['mut_pdb_icode']).strip() + 
                         mutation[-1])
                    LOGGER.info("ddG will be run at %s of %s"%(_params['pdb_mutation'],_params['struct_filename']))
                    df_all_jobs = df_all_jobs.append(makejob_ddg_monomer(_params,),ignore_index = True)

        
        # Launch sequence and spatial jobs on a minimally overlapping subset of structures
        # These are the pathprox jobs, which take a lot longer
        dfSwiss = ci_df.loc[ci_df['label'] == 'swiss']
        dfModbase = ci_df.loc[ci_df['label'] == 'modbase']
        dfOther = ci_df.loc[(ci_df['label'] != 'modbase') & (ci_df['label'] != 'swiss')]
        ci_df    = pd.concat([repr_subset(dfOther),repr_subset(dfSwiss),repr_subset(dfModbase)])

        #Recalculate pdb_muts[i] positions because this new df is quite different
        if mutation:
            pdb_mut_tuples = ci_df[['mut_pdb_res','mut_pdb_icode']].apply(tuple,axis=1)
            pdb_muts = [(mutation[0] +str(int(mut_pdb_res)) +str(mut_pdb_icode.strip()) + mutation[-1]) if (mut_pdb_res and (not np.isnan(mut_pdb_res))) else None for mut_pdb_res,mut_pdb_icode in pdb_mut_tuples]
            transcript_muts = [mutation[0]+str(int(transcript_pos))+mutation[-1] if not np.isnan(transcript_pos) else None for transcript_pos in ci_df['trans_mut_pos'].values]
        else:
            pdb_muts = [None]*len(ci_df)
            transcript_muts = [None]*len(ci_df)
 
        # life easier when the dataframe knows everything in play
        ci_df['transcript_mut'] = transcript_muts
        ci_df['pdb_mutation'] = pdb_muts

        # Command line for  pathprox analyses on all structures

        # Where we can Send pathprox all the multimeric biounit complex as one flavor pathprox job
        # Certainly send pathprox all chains for single-chain processing, even if redundant to multimer
        df_pathprox_single_chain_jobs = makejobs_pathprox_df(params,ci_df[ci_df['mers'] == 'monomer'],multimer=False)
        df_all_jobs = df_all_jobs.append(df_pathprox_single_chain_jobs)
        del df_pathprox_single_chain_jobs

        df_pathprox_biounit_jobs = makejobs_pathprox_df(params,ci_df[ci_df['mers'] != 'monomer'],multimer=True) 
        df_all_jobs = df_all_jobs.append(df_pathprox_biounit_jobs)
        del df_pathprox_biounit_jobs

        
    # whether we had some structures or not, we have our completed workplan to save - and then exit  
    workplan_filename = "%s/%s_%s_%s_workplan.csv"%(mutation_dir,gene,refseq,mutation)

    workstatus_filename = "%s/%s_%s_%s_workstatus.csv"%(mutation_dir,gene,refseq,mutation)

    df_all_jobs.set_index('uniquekey',inplace=True);
    df_all_jobs.sort_index().to_csv(workplan_filename,sep='\t')
    LOGGER.info("Workplan written to %s"%workplan_filename)

    if os.path.exists(workstatus_filename):
        LOGGER.warning("Removing prior workstatus file: %s"%workstatus_filename)
        os.remove(workstatus_filename)

    # Close out the local log file for this mutation
    local_fh.flush()
    local_fh.close()
    LOGGER.removeHandler(local_fh)
    return df_all_jobs,workplan_filename,ci_df,df_dropped,log_filename
    
def makejob_DigenicAnalysis(params,options:str):
    # Secondary structure prediction and sequence annotation
    return makejob("DigenicAnalysis", # Flavor
        "python /dors/capra_lab/projects/psb_collab/UDN/souhrid_scripts/get_digenic_metrics.py", # Command
        params, # Parameters for dataframe
        options) # Command line options

def plan_casewide_work(original_case_xlsx_file):
    # Use the global database connection
    global args
    global io

    LOGGER.info("Planning case-wide work for %s"%args.project)
    # 2017-10-03 Chris Moth modified to key off NT_ refseq id

    casewideString = "casewide"
    casewide_dir = os.path.join(collaboration_dir,casewideString);
    psb_permissions.makedirs(casewide_dir)
    casewide_log_dir = casewide_dir # os.path.join(casewide_dir,"log")
    psb_permissions.makedirs(casewide_log_dir)

    # pw_name = pwd.getpwuid( os.getuid() ).pw_name # example jsheehaj or mothcw
    log_filename = os.path.join(casewide_log_dir,"psb_plan.log")

    if oneMutationOnly:
        sys.stderr.write("psb_plan log file is %s\n"%log_filename)
    else:
        LOGGER.info("Additionally logging to file %s"%log_filename)
    needRoll = os.path.isfile(log_filename)

    local_fh = RotatingFileHandler(log_filename, backupCount=7)
    formatter = logging.Formatter('%(asctime)s %(levelname)-4s [%(filename)20s:%(lineno)d] %(message)s',datefmt="%H:%M:%S")
    local_fh.setFormatter(formatter)
    local_fh.setLevel(logging.INFO)
    LOGGER.addHandler(local_fh)

    if needRoll:
        local_fh.doRollover()

    # whether we had some structures or not, we have our completed workplan to save - and then exit  
    workplan_filename = os.path.join(casewide_dir,"%s_workplan.csv"%casewideString)

    workstatus_filename = os.path.join(casewide_dir,"%s_workstatus.csv"%casewideString)

    df_all_jobs = pd.DataFrame()
    params = {"config":args.config,"userconfig":args.userconfig,"collab":config_dict['collaboration'],"case":args.project,
                        "project":args.project,
                        "pdbid":"N/A",
                        "chain":"N/A",
                        "mers":"N/A",
                        "gene":casewideString,"refseq":"N/A","unp":"N/A",'mutation_dir':casewide_dir,
                        "mutation":"N/A",
                        "transcript_mutation":"N/A",
                        "pdb_mutation":"N/A"}


    digenic_options = "-e %s"%original_case_xlsx_file

    
    # We run one sequence analysis on each transcript and mutation point.. Get that out of the way
    df_all_jobs = df_all_jobs.append(makejob_DigenicAnalysis(params,digenic_options),ignore_index=True)

    df_all_jobs.set_index('uniquekey',inplace=True);
    df_all_jobs.sort_index().to_csv(workplan_filename,sep='\t')
    LOGGER.info("Workplan written to %s"%workplan_filename)

    if os.path.exists(workstatus_filename):
        LOGGER.warning("Removing prior workstatus file: %s"%workstatus_filename)
        os.remove(workstatus_filename)

    # Close out the local log file for this mutation
    local_fh.flush()
    local_fh.close()
    LOGGER.removeHandler(local_fh)
    return df_all_jobs,workplan_filename,log_filename
        

if oneMutationOnly:
    print("Planning work for one mutation only: %s %s %s %s"%(args.project,args.entity,args.refseq,args.mutation))
    df_all_jobs,workplan_filename,df,df_dropped,log_filename = plan_one_mutation(args.entity,args.refseq,args.mutation)
    print("Workplan of %d jobs written to:"%len(df_all_jobs),workplan_filename)
    print("%3d structures/models will be processed."%len(ci_df))
    print("%3d structures/models were considered, but dropped."%len(df_dropped))
    print("Full details in %s",log_filename)
else:
    original_case_xlsx_file = os.path.join(collaboration_dir,"%s.xlsx"%args.project)

    # The exception to the rule is/are the case-wide jobs that the pipeline launches.
    # if not os.path.exists(phenotypes_filename):
    #     LOGGER.critical("File %s was not created from the UDN report."%phenotypes_filename)
    #     LOGGER.critical("Thus, psb_genedicts.py will NOT be part of the planned work")
    # else:
    df_all_jobs,workplan_filename,log_filename = plan_casewide_work(original_case_xlsx_file)

    fulldir, filename = os.path.split(log_filename)
    fulldir, casewide_dir = os.path.split(fulldir)
    fulldir, project_dir = os.path.split(fulldir)
    print(" %4d casewide jobs will run.  See: $UDN/%s"%(len(df_all_jobs),os.path.join(project_dir,casewide_dir,filename)))

    # Now plan the per-mutation jobs
    udn_csv_filename = os.path.join(collaboration_dir,"%s_missense.csv"%args.project)
    print("Retrieving project mutations from %s"%udn_csv_filename)
    df_all_mutations = pd.read_csv(udn_csv_filename,sep=',',index_col = None,keep_default_na=False,encoding='utf8',comment='#',skipinitialspace=True)
    print("Work for %d mutations will be planned"%len(df_all_mutations))

    if 'unp' not in df_all_mutations.columns:
        df_all_mutations['unp'] = None
        LOGGER.warning("Populating unp column of missense from gene and refseq")
        unps_from_refseq_gene = {}
        for index,row in df_all_mutations.iterrows():
            unpsForRefseq = PDBMapProtein.refseqNT2unp(row['refseq'])
            if len(unpsForRefseq):
                      unp    = ','.join(PDBMapProtein.refseqNT2unp(row['refseq']))

            if unp == None:
                if refseq:
                      LOGGER.warning("Could not map refseq input [%s] uniprot ID (or it was missing).  Gene_refseq input in missense.csv file is: %s"%(row['refseq'],row['gene']))
                else:
                      LOGGER.warning("Could not map refseq input to uniprot ID (or it was missing).  Gene_refseq input in missense.csv file is: %s"%row['gene'])
                refseq = "RefSeqNotFound_UsingGeneOnly"
                unp    = PDBMapProtein.hgnc2unp(gene)
                unps_from_refseq_gene[index] = unp

        for index in unps_from_refseq_gene:
                df_all_mutations.loc[index]['unp'] = unps_from_refseq_gene[index]
    else: # 'unp' is in the original columns - make sure the unps are known and queryable
        for unp in df_all_mutations['unp']:
            if not PDBMapProtein.unp2uniparc(unp):
                if '-' in unp and PDBMapProtein.unp2uniparc(unp.split('-')[0]):
                   msg = "unp %s appears to be an invalid isoform of %s as it references no UniParc sequence"%(unp,unp.split('-')[0])
                else:
                    msg = "unp %s was not found in the uniprot IDMapping file with a UniParc AA sequence"%unp
                LOGGER.critical(msg)
                sys.exit(msg)
                
    ui_final_table = pd.DataFrame()
    for index,row in df_all_mutations.iterrows():
        print("Planning %3d,%s,%s,%s,%s"%(index,row['gene'],row['refseq'],row['mutation'],row['unp'] if 'unp' in row else "???"))
        if ('unp' in row) and ('user_model' in row):
            print("....Including user_model %s"%row['user_model'])
        ui_final = {}
        for f in ['gene','refseq','mutation','unp']:
            ui_final[f] = row[f]
    
        df_all_jobs,workplan_filename,df_structures,df_dropped,log_filename = plan_one_mutation(index,row['gene'],row['refseq'],row['mutation'],row['user_model'] if 'user_model' in row and row['user_model'] else None,unp = row['unp'] if 'unp' in row else None)

        fulldir, filename = os.path.split(log_filename)
        fulldir, mutation_dir = os.path.split(fulldir)
        fulldir, project_dir = os.path.split(fulldir)
        print(" %4d structures retained  %4d dropped. %4d jobs will run.  See: $UDN/%s"%(len(df_structures),len(df_dropped),len(df_all_jobs),os.path.join(project_dir,mutation_dir,filename)))
        ui_final['retained'] = len(df_structures)
        ui_final['dropped'] = len(df_dropped)
        ui_final['jobs'] = len(df_all_jobs)
        ui_final['planfile'] = os.path.join(mutation_dir,filename)
        ui_final_table = ui_final_table.append(ui_final,ignore_index=True)

    myLeftJustifiedGene = lambda x: '%-8s'%x
    myLeftJustifiedRefseq = lambda x: '%-14s'%x
    myLeftJustifiedPlanfile = lambda x: '%-40s'%x
    myLeftJustifiedUNP = lambda x: '%-9s'%x
    final_structure_info_table = ui_final_table.to_string(columns=['gene','refseq','mutation','unp','retained','dropped','jobs','planfile'],float_format="%1.0f",justify='center',
                 formatters={'gene': myLeftJustifiedGene, 'refseq': myLeftJustifiedRefseq,'unp': myLeftJustifiedUNP, 'planfile': myLeftJustifiedPlanfile})
    print(("Structure Report\n%s"%final_structure_info_table))
    LOGGER.info("%s",final_structure_info_table)

    # It is so easy to forget to create this phenotypes file - so remind user again!
    # if not os.path.exists(phenotypes_filename):
    #     LOGGER.critical("Reminder: File %s was not created from the UDN report."%phenotypes_filename)
    #     LOGGER.critical("          Thus, psb_genedicts.py will NOT be part of the planned casewide work")
