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

import logging,os,pwd,sys,grp,stat
from logging.handlers import RotatingFileHandler
from logging import handlers

sh = logging.StreamHandler()
logger = logging.getLogger()
logger.addHandler(sh)

# Now that we've added streamHandler, basicConfig will not add another handler (important!)
log_format_string = '%(asctime)s %(levelname)-4s [%(filename)16s:%(lineno)d] %(message)s'
date_format_string = '%H:%M:%S'
log_formatter = logging.Formatter(log_format_string,date_format_string)

logger.setLevel(logging.DEBUG)
sh.setLevel(logging.WARNING)
sh.setFormatter(log_formatter)


import pandas as pd
pd.options.mode.chained_assignment = 'raise'
pd.set_option("display.max_columns",100)
pd.set_option("display.width",1000)
import numpy as np
import time
from get_structures import get_pdbs,get_modbase_swiss
import argparse,configparser
from glob import glob
from pdbmap import PDBMapIO,PDBMapProtein,PDBMapSwiss
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
cmdline_parser.add_argument("entity",nargs='?',type=str,help="Omit or Gene ID or SwissProt AC (ex. TTN)")
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


required_config_items = ["dbhost","dbname","dbuser","dbpass",
  "collaboration",
  "idmapping",
  "interpro_dir",
  "pdb_dir",
  "sec2prim",
  "modbase2013_dir",
  "modbase2013_summary",
  "modbase2016_dir",
  "modbase2016_summary",
  "output_rootdir",
  "swiss_dir",
  "swiss_summary"]

# parser.add_argument("udn_excel",type=str,help="Raw input UDN patient report (format: xls/xlsx)")
# parser.add_argument("udn_csv",type=str,help="Parsed output pipeline filename (e.g. filename.csv)")

config,config_dict = psb_config.read_config_files(args,required_config_items)

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

logger.info("Loading UniProt ID mapping...")

io = None

def detect_entity(io,entity):
  logger.info("\nAttempting to identify entity: %s..."%entity)
  etype = io.detect_entity_type(entity)
  if etype == "unp":
    unp  = entity
    hgnc = PDBMapProtein.unp2hgnc(unp)
    logger.info("UniProt: %s => HGNC: %s..."%(unp,hgnc))
  elif etype and etype!="structure" and etype!="model": # HGNC gene ID
    hgnc = entity
    unp  = etype
    logger.info("HGNC: %s => UniProt: %s..."%(hgnc,unp))
  else:
    logger.error("TERMINATING: Unrecognized ID %s in detect_entity(). Gene name or UniProt/Swiss-Prot AC required."%entity)
    sys.exit(1)
  return unp,hgnc


config_dict_reduced = {x:config_dict[x] for x in required_config_items}

config_dict = config_dict_reduced

config_dict_shroud_password = {x:config_dict[x] for x in required_config_items}
dbpass = config_dict.get('dbpass','?')
config_dict_shroud_password['dbpass'] = '*' * len(dbpass)

oneMutationOnly = args.entity or args.refseq or args.mutation

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
  logger.addHandler(global_fh)

  if needRoll:
    global_fh.doRollover()

# print "Command Line Arguments"
# pprint.pprint(vars(args))
logger.info("Command Line Arguments:\n%s"%pprint.pformat(vars(args)))
logger.info("Configuration File parameters:\n%s"%pprint.pformat(config_dict_shroud_password))

config_pathprox_dict = dict(config.items("PathProx"))  
logger.info("Pathprox config entries:\n%s"%pprint.pformat(config_pathprox_dict))

def pathprox_config_to_argument(disease1_or_2_or_neutral,default):
  variants_filename_key = "%s_variant_filename"%disease1_or_2_or_neutral
  if variants_filename_key in config_pathprox_dict:
    config_str = config_pathprox_dict[variants_filename_key]
    return "--pathogenic %s --pathogenic_label %s"%(config_str,config_pathprox_dict["%s_variant_sql_label"%disease1_or_2_or_neutral]) if 'neutral' not in variants_filename_key else "--neutral %s"%config_str
  else:
    variants_sql_label_key = "%s_variant_sql_label"%disease1_or_2_or_neutral
    if variants_sql_label_key not in config_pathprox_dict:
      logger.critical("A configuration file must define a variant set for %s or %s in the [PathProx] section"%(variants_sql_label_key,variants_filename_key))
    config_str = config_pathprox_dict[variants_sql_label_key]
    if config_str == "clinvar":
      return "--add_pathogenic"
    if config_str == "tcga":
      return "--add_tcga"
    if config_str == "cosmic":
      return "--add_cosmic"
    if config_str == "1kg3":
      return "--add_1kg3"
    if config_str == "gnomad":
      return "--add_gnomad"
    if config_str == "exac":
      return "--add_exac"
    if config_str == "ERROR":
      return "BAD VALUE IN PATHPROX CONFIG FILE SECTION OR CALLER"
    logger.critical("pathprox variant config of [%s] was not found. Applying default of [%s]",config_str,default)
    return pathprox_config_to_argument(default,"ERROR")

pathprox_arguments = {
  'disease1': pathprox_config_to_argument('disease1',"clinvar"),
  'disease2': pathprox_config_to_argument('disease2',"COSMIC"),
  'neutral' : pathprox_config_to_argument('neutral',"exac")
}

logger.info("Pathprox Disease 1 command line argument: %s"%pathprox_arguments['disease1'])
logger.info("Pathprox Disease 2 command line argument: %s"%pathprox_arguments['disease2'])
logger.info("Pathprox Neutral   command line argument: %s"%pathprox_arguments['neutral'])


logger.info("Results for patient case %s will be rooted in %s"%(args.project,collaboration_dir))

swiss_filename = os.path.join(config_dict['swiss_dir'],config_dict['swiss_summary'])
if not infoLogging:
  print("Loading swiss model JSON metadata from %s"%swiss_filename)
logger.info("Loading swiss model JSON metadata from %s"%swiss_filename)

# import cProfile
# cProfile.run("PDBMapSwiss.load_swiss_INDEX_JSON(config_dict['swiss_dir'],config_dict['swiss_summary'])")
# sys.exit(0)
PDBMapSwiss.load_swiss_INDEX_JSON(config_dict['swiss_dir'],config_dict['swiss_summary'])
# In order to set the directories, we need the gene name worted out, interestingly enough
if not infoLogging:
  print("Loading idmapping file from %s"%config_dict['idmapping'])

logger.info("Loading idmapping")
PDBMapProtein.load_idmapping(config_dict['idmapping'])
logger.info("Loading done")

# print "Establishing connection with the PDBMap database..."
if (not io):
  io = PDBMapIO(config_dict['dbhost'],config_dict['dbuser'],config_dict['dbpass'],config_dict['dbname']) # ,slabel=args.slabel)

def repr_subset(df):
  # Sort order: Mutation Coverage, PDBs, Length (rounded), Resolution/Identity
  tdf = df.copy()
  # Convert these sequence start/end fields to ints
  tdf['rndlen'] = (tdf['Seq End'] - tdf['Seq Start']).round(-1)
  tdf['mutcov'] = tdf['Distance to Boundary']==0
  tdf['ispdb']  = tdf['Resolution (PDB)'].notnull()
  tdf = tdf.sort_values(by=['mutcov','ispdb','rndlen','Resolution (PDB)','Seq Identity'],ascending=[False,False,False,True,False])
  nr_cov = pd.DataFrame()
  cov = set([])
  template_set = set([])
  for _,row in tdf.iterrows():
    s_cov = list(range(row['Seq Start'],row['Seq End']))
    # Test if <10% of residues in *this model* overlap current coverage
    if row['PDB Template'] is None:
      templatePDB = row['Structure ID'] 
    else:
      templatePDB = row['PDB Template']
    if len(templatePDB) > 4:
      templatePDB = templatePDB[:4]
    if ( (templatePDB and (not templatePDB in template_set)) or
         (not cov or len(set(s_cov).intersection(cov))/float(len(s_cov))<0.1)):
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
      logger.info("WARNING: Residue %s is missing in %s %s.%s."%(trans_pos,label,sid,chain))
      return np.nan

  return df["chain_seqid"].values[0]

def makejob(flavor,command,params,options,cwd=os.getcwd()):
  job = {}
  job['flavor'] = flavor
  job['config'] = args.config
  job['userconfig'] = args.userconfig
  job['command'] = command
  job['options'] = options
  job['gene'] = params['gene']
  job['refseq'] = params['refseq']
  job['method'] = params.get('method','')
  job['mutation'] = params['mutation']
  job['project'] = args.project
  job['unp'] = params['unp']
  
  job['pdbid'] = params.get('pdbid','')

  # If our structure name ends in .pdb, that's fine - but trim that out for the formation
  # of our output directory and unique key
  dot_pdb = job['pdbid'].find(".pdb")
  if dot_pdb > 2:
    job['pdbid'] = job['pdbid'][0:dot_pdb]
  job['chain'] = params.get('chain','')
  if "SequenceAnnotation" in flavor and job['pdbid'] == 'N/A':
    job['uniquekey'] = "%s_%s_%s_%s"%(job['gene'],job['refseq'],job['mutation'],job['flavor'])
    job['outdir'] = params['mutation_dir']
  elif "MakeGeneDictionaries" in flavor and job['pdbid'] == 'N/A':
    job['uniquekey'] = job['flavor']
    job['outdir'] = params['mutation_dir']
  else: # Most output directories have pdbid and chain suffixes
    if not job['chain'].strip() or job['chain'] == "''" or job['chain'] == "' '":
      pdb_chain_segment = job['pdbid']
    else:
      pdb_chain_segment = "%s_%s"%(job['pdbid'],job['chain'])
    job['outdir'] = os.path.join(params['mutation_dir'] , pdb_chain_segment )
    job['uniquekey'] = "%s_%s_%s_%s_%s"%(job['gene'],job['refseq'],job['mutation'],pdb_chain_segment,job['flavor'])
  job['cwd'] = cwd
  # too much detail logger.debug(pprint.pformat(job))
  return pd.Series(job)

# Jobs for PathProx analyses
def makejobs_pathprox(params,method,sid,ch,pdb_mut_i):
  params['method'] = method
  params['pdbid'] = sid
  params['chain'] = "%s"%ch  # Quote it in case the IDentifier is a space, as in many models
  if (not pdb_mut_i) or (pdb_mut_i == 'None'):
    params['pp_mut'] = ''
  else:
    params['pp_mut'] = pdb_mut_i

  params['sqlcache'] = os.path.join(params['mutation_dir'],"sqlcache")

  params_transcript = ''
  if 'transcript' in params and params['transcript']:
    params_transcript = " --isoform '%s'"%params['transcript']

  df_jobs_to_run = pd.DataFrame()
  # ClinVa params['disease1argument'] = pathprox_disease1argument

  #  Make jobs for Disease 1 (maybe ClinVar) PathProx and Disease 2 (maybe COSMIC) PathProx
  for disease_1or2 in [ 'disease1', 'disease2' ]:
      disease_variant_sql_label = config_pathprox_dict['%s_variant_sql_label'%disease_1or2]
      params['diseaseArgument'] = pathprox_arguments[disease_1or2]
      params['neutralArgument'] = pathprox_arguments['neutral']
      # import pdb; pdb.set_trace()
      j = makejob("PathProx_%s"%disease_variant_sql_label,"pathprox2.py",params,
        ("-c %(config)s -u %(userconfig)s %(pdbid)s %(refseq)s %(pp_mut)s " +
         "--chain=%(chain)s %(diseaseArgument)s %(neutralArgument)s --radius=D " + 
         "--sqlcache=%(sqlcache)s " + 
         "--overwrite")%params + params_transcript)
      df_jobs_to_run = df_jobs_to_run.append(j,ignore_index=True)

  return df_jobs_to_run

def makejob_udn_sequence(params):
  # Secondary structure prediction and sequence annotation
  return makejob("SequenceAnnotation","udn_pipeline2.py",params,   #command
          ("--config %(config)s --userconfig %(userconfig)s " +
           "--project %(collab)s --patient %(project)s --gene %(gene)s --transcript %(transcript_mutation)s")%params)

def makejob_udn_structure(params,method,sid,ch,pdb_mut_i):
  params['method'] = method
  params['pdbid'] = sid
  params['chain'] = ch 
  params['pdb_mutation'] = pdb_mut_i
  # ddg_monomer (plus redundant secondary structure/uniprot annotation)
  return makejob("ddG","udn_pipeline2.py",params,   #command
          ("--config %(config)s --userconfig %(userconfig)s " +
           "--project %(collab)s --patient %(project)s --gene %(gene)s --structure %(pdbid)s --chain %(chain)s --mutation %(pdb_mutation)s --transcript %(transcript_mutation)s")%params) # options

# Save ourselves a lot of trouble with scoping of df_dropped by declaring it global with apology
# Fix this in python 3 - where we have additional local scoping options
df_dropped = None

def plan_one_mutation(entity,refseq,mutation,user_model=None,unp=None):
  # Use the global database connection
  global io
  global args
  global oneMutationOnly
  global df_dropped

  logger.info("Planning mutation for %s %s %s %s",args.project,entity,refseq,mutation)
  if user_model:
    logger.info("Additional User model %s requested",user_model)

  # 2017-10-03 Chris Moth modified to key off NT_ refseq id
  unps = PDBMapProtein.refseqNT2unp(refseq)
  if unps:
    newunp = unps[0]
    if unp and (newunp != unp):
      logger.warning("unp determined from refseq is %s which does not match %s",newunp,unp)
    unp = newunp
    gene = PDBMapProtein.unp2hgnc(unp)
    del unps
  else: # Sort this out the old way
    if unp:
      logger.warning("Setting refseq to NA and using unp=%s",unp) 
      gene = PDBMapProtein.unp2hgnc(unp)
      refseq = 'NA'
    else:
      logger.warning("Refseq %s does not map to a uniprot identifier.  Attempting map of gene %s\n"%(refseq,entity))
      io = PDBMapIO(config_dict['dbhost'],config_dict['dbuser'],config_dict['dbpass'],config_dict['dbname']) # ,slabel=args.slabel)
      unp,gene = detect_entity(io,entity)

  if not gene:
    logger.critical("A gene name was not matched to the refseq or uniprot ID.  Gene will be set to %s",entity)
    gene=entity
  
  mutation_dir = os.path.join(collaboration_dir,"%s_%s_%s"%(gene,refseq,mutation))
  psb_permissions.makedirs(mutation_dir)
  mutation_log_dir = mutation_dir
  psb_permissions.makedirs(mutation_log_dir)

  # pw_name = pwd.getpwuid( os.getuid() ).pw_name # example jsheehaj or mothcw
  log_filename = os.path.join(mutation_log_dir,"psb_plan.log")


  if oneMutationOnly:
    sys.stderr.write("psb_plan log file is %s\n"%log_filename)
  else:
    logger.info("Additionally logging to file %s"%log_filename)
  needRoll = os.path.isfile(log_filename)

  local_fh = RotatingFileHandler(log_filename, backupCount=7)
  formatter = logging.Formatter('%(asctime)s %(levelname)-4s [%(filename)20s:%(lineno)d] %(message)s',datefmt="%H:%M:%S")
  local_fh.setFormatter(formatter)
  local_fh.setLevel(logging.INFO)
  logger.addHandler(local_fh)

  if needRoll:
    local_fh.doRollover()
  
  logger.info("MySQL: Marking psb_collab.MutationSummary incomplete for %s %s %s %s %s"%
      (args.project,gene,refseq,mutation,unp))
  
  # Ensure that this mutation is tracked in the database
  # Even though typically already added in Mtuation Summary
  c = io._con.cursor()
  # No-op if this mutation is already recorded in the database
  sql  = "INSERT IGNORE INTO psb_collab.MutationSummary "
  sql += "(collab,project,gene,refseq,mutation,uniprot) "
  sql += "VALUES (%s,%s,%s,%s,%s,%s) "
  c.execute(sql,(config_dict['collaboration'],args.project,gene,refseq,mutation,unp))
  c.close()

  """ 
  # Update the database to include the UniProt AC for this gene
  # Update to incomplete status, in case this is a rerun of a failed report
  c   = io._con.cursor()
  sql  = "UPDATE psb_collab.MutationSummary SET uniprot=%s, complete=0, comments='' WHERE "
  sql += "collab=%s AND project=%s AND gene=%s AND refseq=%s AND mutation=%s"
  c.execute(sql,(unp,config_dict['collaboration'],args.project,gene,refseq,mutation))
  c.close()
  """
  
  # Extract the mutation position
  if mutation:
    mut_pos = int(mutation[1:-1])
  else:
    mut_pos = None
  
  logger.info("Calling get_structures() to load PDB/ModBase/Swiss structures of %s %s unp=%s from MySQL..."%(gene,refseq,unp))
  logger.info("Identifying PDB/ModBase/Swiss structures of %s %s unp=%s..."%(gene,refseq,unp))
  df_pdbs = get_pdbs(io,unp,mut_pos,maybe_threshold = args.maybe)
  logger.info("%d pdbs identified."%len(df_pdbs))
  df_modbase_swiss = get_modbase_swiss(io,unp,mut_pos,maybe_threshold = args.maybe)
  logger.info("%d modbase swiss identified."%len(df_modbase_swiss))

  # Funny bug - but if df_pdbs is empty, the pd.concat turns all the second dF_modbase_swiss columns from ints to floats
  if len(df_pdbs) == 0:
    df = df_modbase_swiss.copy()
  elif len(df_modbase_swiss) == 0:
    df = df_pdbs.copy()
  else:
    df = pd.concat([df_pdbs,df_modbase_swiss],ignore_index = True)

  # df = get_structures(io,unp,mut_pos,maybe_threshold = args.maybe)
  # logger.info("%d structures identified."%len(df))
  
  # Use shorter header description for tighter htm output
  df.rename(columns={'Analysis Possible?': 'Analyzable?','Sequence Start': 'Seq Start','Sequence End': 'Seq End','Sequence Identity': 'Seq Identity'}, inplace=True)
  
  # Add in Swiss-model quality metrics from the REMARK 3 entries of the model file
  # This probably should be moved off to tables from load-time
  # logger.warning("DF columns = %s"%str(df.columns))

  for i in range(len(df)):
    row = df.iloc[i]
    if row['Label'] == 'pdb':
      df.iat[i,df.columns.get_loc('Seq Identity')] = float('nan')
    elif row['Label'] == 'swiss':
      remark3_metrics = PDBMapSwiss.load_REMARK3_metrics(row['Structure ID'])
      df.iat[i,df.columns.get_loc('Seq Identity')] = float(remark3_metrics['sid'])
      
  # Add the canonical sequence position for this mutation
  df["Transcript Pos"] = mut_pos
  
  
  # Query the PDB position for this position in the transcript
  if not df.empty and mutation:
    logger.info("Extracting the position of %s in each PDB/ModBase/Swiss structure..."%mutation)
    df["PDB Pos"] = df[["Label","Structure ID","Chain","Transcript","Transcript Pos",
                          "Distance to Boundary"]].apply(get_pdb_pos,axis=1)
  
    # This appears to be a kind of cleanup.  get_structures() above will already have marked
    # the structures as "Maybe" in case of "Distance to Boundary" < 20
    # This appears to be a situation where we have an experimental structure and the
    # PDB Pos of the mutation is missing.
    maybe_list = np.where((df["Distance to Boundary"] == 0) & (df["PDB Pos"].isnull()))
    if len(maybe_list) > 0:
      maybe_list = maybe_list[0]
      if len(maybe_list) > 0:
        df.iloc[maybe_list,df.columns.get_loc("Analyzable?")] = "Maybe"
    
    # df.iloc[(df["PDB Pos"].isnull()) & (df["Distance to Boundary"]==0),df.columns.get_loc("Analyzable?")] = "Maybe"
  elif not mutation:
    df["PDB Pos"] = np.nan
  
  if not df.empty: #ChrisMoth added this line and indented as the code below that makes no sense for empty df
    df["Gene"] = gene
    df = df[["Gene","Label","Structure ID","Chain","Method","Resolution (PDB)","Seq Identity",
            "PDB Template","Residues","Seq Start","Seq End",
            "Transcript Pos","PDB Pos","Distance to Boundary","Analyzable?"]]
  
    df_dropped = pd.DataFrame(columns=df.columns.values.tolist() + ["Reason"])
  
    # We now take pains to rid our candidate structures of models which are clearly redundant/degenerate
    # The algorithm retains all experimental structures
    def drop_df_row(k,reason,keeping = -1):
      global df_dropped
      if keeping >= 0: # Add the other structure ID being kept, to the reason text, if asked 
        reason += ' vs %s.%s'%(df.iloc[keeping]['Structure ID'],df.iloc[keeping]['Chain'])
      logger.info( "Dropping df duplicate row %d because %s %s"%(k,df.iloc[k]['Structure ID'],reason))
 
      # Record the removed model for the drop report
      dropped_model = df.iloc[k].copy()
  
      dropped_model['Reason'] = reason 
      df_dropped = df_dropped.append(dropped_model)
  
      # Remove model from set to process
      df.drop(df.index[k],inplace=True) 

    # Jonathan requested that all "No" structures be dropped - so let's do that first 
    i = 0
    while i < df.shape[0]:
      row_i = df.iloc[i]
      if row_i['Analyzable?'] == 'No':
          drop_df_row(i,"Not Analyzable")
      else:
        i += 1

    i = 0
    j = 0
    while i < df.shape[0]-1:
      row_i = df.iloc[i]
  
      row_i_template = None
      if row_i['PDB Template']: # This is an experimenal structure
        row_i_template = row_i['PDB Template'].split('.')[0].upper() 
  
      j = i+1
      while j < df.shape[0]:
        row_j = df.iloc[j]
  
  
        # If both are experimental structures, then unique care is warranted before delete
        if row_i['Label'] == 'pdb' and row_j['Label'] == 'pdb':
          # If we have a "yes",then drop the maybes and nos
          # If we have a "maybe",then drop the nos
          if row_i['Analyzable?'] != row_j['Analyzable?']:
            if row_i['Analyzable?'] == 'No':
              drop_df_row(i,"is not as analyzable",j)
              break
            elif row_j['Analyzable?'] == 'No':
              drop_df_row(j,"is not as analyzable",i)
              breaktest_SCN11A_psb_plan.log
            elif row_i['Analyzable?'] == 'Maybe':
              drop_df_row(i,"is not as analyzable",j)
              break
            elif row_j['Analyzable?'] == 'Maybe':
              drop_df_row(j,"is not as analyzable",i)
              break
          else: # Both experimental structures _are_ (or are not) analyzable - same level...
            row_i_len = int(row_i['Seq End'] or 0) - int(row_i['Seq Start'] or 0)
            row_j_len = int(row_j['Seq End'] or 0) - int(row_j['Seq Start'] or 0)
            row_i_resolution = float(row_i['Resolution (PDB)'])
            row_j_resolution = float(row_j['Resolution (PDB)'])
            if row_i['Analyzable?'] == 'Yes': # Generally this is good - but some peptides appear 100x in pdb and we do not need all that
              if (row_i['Seq Start'] == row_j['Seq Start']) and (row_i['Seq End'] == row_j['Seq End']):
                # Both penptides cover same region.. Delete lower resolution
                if row_i_resolution > row_j_resolution:
                  drop_df_row(i,"is lower resolution",j)
                  break
                elif row_j_resolution == row_i_resolution:
                  drop_df_row(i,"is same resolution",j)
                  break
                elif row_j_resolution > row_i_resolution:
                  drop_df_row(j,"is lower resolution",i)
                  break
              elif ((row_i_len > row_j_len) 
                     and  (row_i['Seq Start'] <= row_j['Seq Start']) 
                     and (row_i['Seq End'] >= row_j['Seq Start'])):
                drop_df_row(j,"is shorter",i)
                break
              elif ((row_j_len > row_i_len) 
                     and (row_j['Seq Start'] <= row_i['Seq Start']) 
                     and (row_j['Seq End'] >= row_i['Seq Start'])):
                drop_df_row(i,"is shorter",j)
                break
              # Otherwise, the structures might overlap - but these are experimental structures - do not throw away! 
  
            elif row_i['Analyzable?'] == 'No' or row_i['Analyzable?'] == 'Maybe': # Neither is very analyzable
              if int(row_i['Distance to Boundary']) > int(row_j['Distance to Boundary']):
                drop_df_row(i,"is farther from mutation",j)
                break
              elif int(row_i['Distance to Boundary']) < int(row_j['Distance to Boundary']):
                drop_df_row(j,"is farther from mutation",i)
                break
              if row_i_len > row_j_len:
                drop_df_row(j,"is shorter",i)
                break
              elif row_j_len > row_i_len:
                drop_df_row(i,"is shorter",j)
                break
              else: # Both 'Np' structures have same length.   Eliminate by resolution difference
                if row_i_resolution > row_j_resolution:
                  drop_df_row(i,"is higher resolution",j)
                  break
                elif row_j_resolution > row_i_resolution:
                  drop_df_row(j,"is higher resolution",i)
                  break
                else: # having tried all else - just drop the higher struct id
                  row_i_structid = row_i['Structure ID']
                  row_j_structid = row_j['Structure ID']
                  if (row_i_structid > row_j_structid):
                    drop_df_row(i,"is same length and resolution, but greater structid ID",j)
                  else:
                    drop_df_row(j,"is same length and resolution, but greater than or same structid ID",i)
                  break
            else: # Both structures are "yes" to analyzable - and much greater caution before dropping one
              pass
        # Don't remove anything if comparing an experimental structure to a model
        elif row_i['Label'] == 'pdb' or row_j['Label'] == 'pdb':
          pass
        else: # We are going to compare two models (not comparing experimental structures)
          row_j_template = row_j['PDB Template'].split('.')[0].upper()
          if (row_i_template == row_j_template): # If both models use same template, we need to drop one of them
            # If both are ENSP models, then we need to just get rid of the 2013 one, and keep the 2016 one containing a period
            if row_i['Label'] == 'modbase' and row_j['Label'] == 'modbase':
              if (('.' in row_i['Structure ID'] and '.' in row_j['Structure ID']) or # if Both structures modbase16
                  ('.' not in row_i['Structure ID'] and '.' not in row_j['Structure ID'])): # or if Both structures modbase13
                row_i_len = int(row_i['Seq End'] or 0) - int(row_i['Seq Start'] or 0)
                row_j_len = int(row_j['Seq End'] or 0) - int(row_j['Seq Start'] or 0)
                if row_i_len > row_j_len:
                  drop_df_row(j,"is shorter modbase",i)
                elif row_j_len > row_i_len:
                  drop_df_row(i,"is shorter modbase",j)
                else: # The two modbase models are on the same structure - but must be two different chains
                  # The Modbase Structure ID has the chain after the final . in the name
                  row_i_chain = row_i['Structure ID'].split('.')[-1].upper()
                  row_j_chain = row_j['Structure ID'].split('.')[-1].upper()
                  if (row_i_chain > row_j_chain):
                    drop_df_row(i,"is same length, but greater chain ID",j)
                  else:
                    drop_df_row(j,"is same length, but greater than or same chain ID",i)
                break;
              elif '.' in row_i['Structure ID']: # Then row_i is Modbase2016 - so drop row_j
                drop_df_row(j,"is modbase13",i)
                break
              elif '.' in row_j['Structure ID']: # Then row_j is Modbase2016 - so drop row_i
                drop_df_row(i,"is modbase13",j)
                break
            # If one of the same-template models is swiss, other modebase
            # Kepp swiss
            elif row_i['Label'] == 'swiss' and row_j['Label'] == 'modbase':
              drop_df_row(j,"is modbase (keeping swiss)",i)
              break
            elif row_i['Label'] == 'modbase' and row_j['Label'] == 'swiss':
              drop_df_row(i,"is modbase (keeping swiss)",j)
              break
            # If they are both swiss, then keep the longer chain... or lower chain ID if same length chains
            elif row_i['Label'] == 'swiss' and row_j['Label'] == 'swiss':
              row_i_len = int(row_i['Seq End'] or 0) - int(row_i['Seq Start'] or 0)
              row_j_len = int(row_j['Seq End'] or 0) - int(row_j['Seq Start'] or 0)
              if row_i_len > row_j_len:
                drop_df_row(j,"is shorter swissmodel",i)
              elif row_j_len > row_i_len:
                drop_df_row(i,"is shorter swissmodel",j)
              else: # The two swiss models are on the same structure - but must be two different chains
                # The Swiss Structure ID has the chain after the final . in the name
                row_i_chain = row_i['Structure ID'].split('.')[-1].upper()
                row_j_chain = row_j['Structure ID'].split('.')[-1].upper()
                if (row_i_chain > row_j_chain):
                  drop_df_row(i,"is same length, but greater chain ID",j)
                else:
                  drop_df_row(j,"is same length, but greater than or same chain ID",i)
              break;
        # end of dup drops in models
  
        j += 1
      else: # Do not increment i if a row was dropped (and break was called)
        i += 1
  
    # Important: This function drops by index, not the physical row number as the duplicate drop code above
    def drop_df_low_identity_rows(droplist,reason):
      global df_dropped
      if len(droplist) < 1:
        return
      for k in droplist:
        logger.info("Dropping of row %d because %s %s"%(k,df.loc[k]['Structure ID'],reason))
      dropped_models = df.loc[droplist].copy()
      dropped_models['Reason'] = reason 
      df_dropped = df_dropped.append(dropped_models)
      df.drop(droplist,inplace=True) 
  
    models_sub10_seq_identity = df.index[(((df['Label'] == 'modbase') | (df['Label'] == 'swiss')) & (df['Seq Identity'] < 10.0))]
    if len(models_sub10_seq_identity) > 0:
      drop_df_low_identity_rows(models_sub10_seq_identity,"Seq Identity < 10.0 unacceptable")
    models_10to20_seq_identity = df.index[(((df['Label'] == 'modbase') | (df['Label'] == 'swiss')) & ((10.0 <= df['Seq Identity']) & (df['Seq Identity'] < 20.0)))]
    models_20to30_seq_identity = df.index[(((df['Label'] == 'modbase') | (df['Label'] == 'swiss')) & ((20.0 <= df['Seq Identity']) & (df['Seq Identity'] < 30.0)))]
    models_over30_seq_identity = df.index[(((df['Label'] == 'modbase') | (df['Label'] == 'swiss')) & (30.0 <= df['Seq Identity']))]
  
    # Do we have a solid group of anlyzable (mutation covering models?
    analyzable_over_30 = df.index[(df["Analyzable?"] == "Yes") & ((df['Label'] == 'pdb') | (((df['Label'] == 'modbase') | (df['Label'] == 'swiss')) & (30.0 <= df['Seq Identity'])))]
    analyzable_under_30 = df.index[(df["Analyzable?"] == "Yes") & (((df['Label'] == 'modbase') | (df['Label'] == 'swiss')) & (df['Seq Identity'] < 30.0))]
  
    # If we have 3 pretty good models, throw away the low identity ones
    # Unless we are hurting badly for analyzable models, and we have some
    if len(models_over30_seq_identity) > 3:
      if len(analyzable_over_30) == 0 and len(analyzable_under_30) > 0: # In this case (no 'good' structures/models with mutation coverage)
        models_not_analyzable_under_30 = df.index[(((df['Label'] == 'modbase') | (df['Label'] == 'swiss')) & (df['Analyzable?'] != 'Yes') & (df['Seq Identity'] < 30.0))]
        drop_df_low_identity_rows(models_not_analyzable_under_30,"Retaining only low Seq Identity Models with mutation coverage")
      else:
        drop_df_low_identity_rows(models_10to20_seq_identity,"Low Seq Identity Models are not needed")
        drop_df_low_identity_rows(models_20to30_seq_identity,"Low Seq Identity Models are not needed")
  
    elif len(models_20to30_seq_identity) > 3:
      if len(analyzable_over_30) == 0 and len(analyzable_under_30) > 0: # In this case (no 'good' structures/models with mutation coverage)
        models_not_analyzable_under_20 = df.index[(((df['Label'] == 'modbase') | (df['Label'] == 'swiss')) & (df['Analyzable?'] != 'Yes') & (df['Seq Identity'] < 30.0))]
        drop_df_low_identity_rows(models_not_analyzable_under_20,"Retaining only low Seq Identity Models with mutation coverage")
      else:
        drop_df_low_identity_rows(models_10to20_seq_identity,"Low Seq Identity Models are not needed")
  
    # Do we have a solid group of anlyzable (mutation covering models?
    # If we have analyzable models over 20 at all - then drop all models under_20
    analyzable_over_20 = df.index[(df["Analyzable?"] == "Yes") & ((df['Label'] == 'pdb') | (((df['Label'] == 'modbase') | (df['Label'] == 'swiss')) & (20.0 <= df['Seq Identity'])))]
    if len(analyzable_over_20) > 0:
      analyzable_under_20 = df.index[(df["Analyzable?"] == "Yes") & (((df['Label'] == 'modbase') | (df['Label'] == 'swiss')) & (df['Seq Identity'] < 20.0))]
      if len(analyzable_under_20) > 0:
        drop_df_low_identity_rows(analyzable_under_20,"Low Seq Identity covering models not needed")
  
    df = df.sort_values(by=["Distance to Boundary","Method","Residues",
                          "Resolution (PDB)","Seq Identity"],
                          ascending=[True,False,False,True,False])
  
    dropped_models_filename = "%s/%s_%s_%s_dropped_models.csv"%(mutation_dir,gene,refseq,mutation)
    logger.info("Output the %d dropped  ModBase and Swiss models to: %s"%(len(df_dropped),dropped_models_filename))
    df_dropped.to_csv(dropped_models_filename,sep='\t',index=True)
    #
    # End if not df_empty

  ENST_transcript_ids = PDBMapProtein.unp2enst(unp)
  if ENST_transcript_ids:
    logger.info("ENST transcript for %s will be %s",unp,ENST_transcript_ids[0])
    if len(ENST_transcript_ids) > 1:
      logger.warning("%d additional ENST transcripts will be ignored",len(ENST_transcript_ids)-1)

  elif (config_pathprox_dict['disease1_variant_sql_label'] == 'tcga') or \
       (config_pathprox_dict['disease2_variant_sql_label'] == 'tcga'):
    logger.error("The uniprot identifier must map to an ENST transcript ID in order use tcga variants",user_model,model_count)
    sys.exit(1)

  if user_model: # Then add it to our structure report
    if not ENST_transcript_ids:
      logger.error("The uniprot identifier must map to an ENST transcript ID in order to incorporate a user model",user_model,model_count)
      sys.exit(1)

    from Bio.PDB.PDBParser import PDBParser 
    p = PDBParser()
    structure = p.get_structure('user_model',user_model)
    model_count = 0
    chain_count = 0;
    chain_id = ' '
    for model in structure:
      import pdb; pdb.set_trace()
      model_count += 1
      for chain in model:
        chain_id = chain.get_id()
        chain_count += 1
        seq_start = -1 
        seq_end = -1
        for residue in chain:
          resid = residue.id
          assert((resid[0] == None or resid[0] == ' ') and ('User model residues cannot be WAT or HET entries'))
          assert((resid[2] == None or resid[2] == ' ') and ('User model insertion codes must be blank'))
          if seq_start == -1:
            seq_start = residue.id[1]
          seq_end = residue.id[1]

    if model_count != 1:
      logger.error("User model %s must contain one and only one model, not %d",user_model,model_count)
      sys.exit(1)

    if chain_count != 1:
      logger.error("User model %s must contain one and only one chain, not %d",user_model,chain_count)
      sys.exit(1)

    for chain in structure: # This leaves chain referencing the first chain
      break

    df = df.append({'Gene': gene, 'Label': 'UserModel', 'Method': 'UserModel', 'transcript': ENST_transcript_ids[0], 'Structure ID': user_model, 'PDB Template': 'N/A', "Analyzable?": 'Yes', 'Chain': chain_id, "Transcript Pos": mut_pos, "PDB Pos": mut_pos, "Distance to Boundary": 0.0, "Seq Start":seq_start, "Seq End": seq_end},ignore_index=True);
 
 
  # df could be empty or not....
  structure_report_filename = "%s/%s_%s_%s_structure_report.csv"%(mutation_dir,gene,refseq,mutation)
  logger.info("Output the %d retained PDB, ModBase and Swiss structures to: %s"%(len(df),structure_report_filename))
  df.to_csv(structure_report_filename,sep='\t',index=True,na_rep='nan')
  if mutation and len(df) > 0:
    logger.info("Total with coverage of %s: %d"%(mutation,(~df["PDB Pos"].isnull()).sum()))
  
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
            "gene":gene,"refseq":refseq,"unp":unp,'mutation_dir':mutation_dir,
            "mutation":mutation,
            "transcript_mutation":mutation,
            "pdb_mutation":"Error: replace for each structure"}
 
  if ENST_transcript_ids:
    params["transcript"] = ENST_transcript_ids[0]
  # We run one sequence analysis on each transcript and mutation point.. Get that out of the way
  df_all_jobs = df_all_jobs.append(makejob_udn_sequence(params),ignore_index=True)
  
  # logger.info("%d Structures/models of %s %s (unp=%s)"%(len(df),refseq,gene,unp))
  # If no structures at all are available, then the work plan takes a very different course
  if df.empty:
    if args.mutation:
      logger.info("No structures with coverage of the mutation are available. No report will be generated.")
    else:
      logger.info("No structures are available. No report will be generated.")
  
  
    # If the script made it to this point, set this mutation to "complete"
    sql  = "UPDATE psb_collab.MutationSummary SET complete=True WHERE "
    sql += "collab=%s AND project=%s AND gene=%s AND mutation=%s "
    c = io._con.cursor()
    c.execute(sql,(config_dict['collaboration'],args.project,gene,mutation))
    c.close()
  else:
    # Generate a PDB-indexed mutation string for each chain
    # Chains that do not cover the mutation site are set to None in this array
    if mutation:
      pdb_muts = [mutation[0]+str(int(pdb_pos))+mutation[-1] if not np.isnan(pdb_pos) else None for pdb_pos in df["PDB Pos"].values]
      transcript_muts = [mutation[0]+str(int(transcript_pos))+mutation[-1] if not np.isnan(transcript_pos) else None for transcript_pos in df["Transcript Pos"].values]
    else:
      pdb_muts = [None]*len(df)
      transcript_muts = [None]*len(df)
    if mutation:
      logger.info("Evaluating mutation: %s"%(mutation))
      logger.info("Corresponding PDB/ModBase mutations:")
      for i,(sid,ch) in enumerate(df[["Structure ID","Chain"]].values):
        if transcript_muts[i]:
          logger.info("%s\t%s\t%s\tPDB pos: %s"%(sid,ch,transcript_muts[i],pdb_muts[i]))
        else:
          logger.info("%s\t%s\t%s"%(sid,ch,"Not covered by structure."))
  
    # Launch structural analyses on all structures
    for i,(method,sid,ch) in enumerate(df[["Method","Structure ID","Chain"]].values):
      if pdb_muts[i]: # Don't launch the udn_structure (ddG) script if there is no mutation coverage
        df_all_jobs = df_all_jobs.append(makejob_udn_structure(params,method,sid,ch,pdb_muts[i]),ignore_index = True)

    
    # Launch sequence and spatial jobs on a minimally overlapping subset of structures
    # These are the pathprox jobs, which take a lot longer
    dfSwiss = df.loc[df['Label'] == 'swiss']
    dfModbase = df.loc[df['Label'] == 'modbase']
    dfOther = df.loc[(df['Label'] != 'modbase') & (df['Label'] != 'swiss')]
    df    = pd.concat([repr_subset(dfOther),repr_subset(dfSwiss),repr_subset(dfModbase)])

    #Recalculate pdb_muts[i] positions because this new df is quite different
    if mutation:
      pdb_muts = [mutation[0]+str(int(pdb_pos))+mutation[-1] if not np.isnan(pdb_pos) else None for pdb_pos in df["PDB Pos"].values]
      transcript_muts = [mutation[0]+str(int(transcript_pos))+mutation[-1] if not np.isnan(transcript_pos) else None for transcript_pos in df["Transcript Pos"].values]
    else:
      pdb_muts = [None]*len(df)
      transcript_muts = [None]*len(df)
 
    # import pdb; pdb.set_trace() 
    # Launch pathprox analyses on all structures
    df_all_jobs = df_all_jobs.append([makejobs_pathprox(params,method,sid,ch,transcript_muts[i]) \
      for i,(method,sid,ch) in enumerate(df[["Method","Structure ID","Chain"]].values) ],ignore_index = True)
    
  # whether we had some structures or not, we have our completed workplan to save - and then exit  
  workplan_filename = "%s/%s_%s_%s_workplan.csv"%(mutation_dir,gene,refseq,mutation)

  workstatus_filename = "%s/%s_%s_%s_workstatus.csv"%(mutation_dir,gene,refseq,mutation)

  df_all_jobs.set_index('uniquekey',inplace=True);
  df_all_jobs.sort_index().to_csv(workplan_filename,sep='\t')
  logger.info("Workplan written to %s"%workplan_filename)

  if os.path.exists(workstatus_filename):
    logger.warning("Removing prior workstatus file: %s"%workstatus_filename)
    os.remove(workstatus_filename)

  # Close out the local log file for this mutation
  local_fh.flush()
  local_fh.close()
  logger.removeHandler(local_fh)
  return df_all_jobs,workplan_filename,df,df_dropped,log_filename
    
def makejob_MakeGeneDictionaries(params):
  # Secondary structure prediction and sequence annotation
  return makejob("MakeGeneDictionaries","psb_genedicts.py",params,   #command
          ("--config %(config)s --userconfig %(userconfig)s %(project)s")%params)

def plan_casewide_work():
  # Use the global database connection
  global args
  global io

  logger.info("Planning case-wide work for %s"%args.project)
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
    logger.info("Additionally logging to file %s"%log_filename)
  needRoll = os.path.isfile(log_filename)

  local_fh = RotatingFileHandler(log_filename, backupCount=7)
  formatter = logging.Formatter('%(asctime)s %(levelname)-4s [%(filename)20s:%(lineno)d] %(message)s',datefmt="%H:%M:%S")
  local_fh.setFormatter(formatter)
  local_fh.setLevel(logging.INFO)
  logger.addHandler(local_fh)

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
            "gene":casewideString,"refseq":"N/A","unp":"N/A",'mutation_dir':casewide_dir,
            "mutation":"N/A",
            "transcript_mutation":"N/A",
            "pdb_mutation":"N/A"}
  
  # We run one sequence analysis on each transcript and mutation point.. Get that out of the way
  df_all_jobs = df_all_jobs.append(makejob_MakeGeneDictionaries(params),ignore_index=True)
  df_all_jobs.set_index('uniquekey',inplace=True);
  df_all_jobs.sort_index().to_csv(workplan_filename,sep='\t')
  logger.info("Workplan written to %s"%workplan_filename)

  if os.path.exists(workstatus_filename):
    logger.warning("Removing prior workstatus file: %s"%workstatus_filename)
    os.remove(workstatus_filename)

  # Close out the local log file for this mutation
  local_fh.flush()
  local_fh.close()
  logger.removeHandler(local_fh)
  return df_all_jobs,workplan_filename,log_filename
    

if oneMutationOnly:
  print("Planning work for one mutation only: %s %s %s %s"%(args.project,args.entity,args.refseq,args.mutation))
  df_all_jobs,workplan_filename,df,df_dropped,log_filename = plan_one_mutation(args.entity,args.refseq,args.mutation)
  print("Workplan of %d jobs written to:"%len(df_all_jobs),workplan_filename)
  print("%3d structures/models will be processed."%len(df))
  print("%3d structures/models were considered, but dropped."%len(df_dropped))
  print("Full details in %s",log_filename)
else:
  phenotypes_filename = os.path.join(collaboration_dir,"Phenotypes","%s_phenotypes.txt"%args.project)

  # The exception to the rule is/are the case-wide jobs that the pipeline launches.
  if not os.path.exists(phenotypes_filename):
    logger.critical("File %s was not created from the UDN report."%phenotypes_filename)
    logger.critical("Thus, psb_genedicts.py will NOT be part of the planned work")
  else:
    df_all_jobs,workplan_filename,log_filename = plan_casewide_work()

    fulldir, filename = os.path.split(log_filename)
    fulldir, casewide_dir = os.path.split(fulldir)
    fulldir, project_dir = os.path.split(fulldir)
    print(" %4d casewide jobs will run.  See: $UDN/%s"%(len(df_all_jobs),os.path.join(project_dir,casewide_dir,filename)))

  # Now plan the per-mutation jobs
  udn_csv_filename = os.path.join(collaboration_dir,"%s_missense.csv"%args.project)
  print("Retrieving project mutations from %s"%udn_csv_filename)
  df_all_mutations = pd.read_csv(udn_csv_filename,sep=',',index_col = None,keep_default_na=False,encoding='utf8')
  print("Work for %d mutations will be planned"%len(df_all_mutations))

  ui_final_table = pd.DataFrame()
  for index,row in df_all_mutations.iterrows():
    print("Planning %3d,%s,%s,%s,%s"%(index,row['gene'],row['refseq'],row['mutation'],row['unp'] if 'unp' in row else "???"))
    if ('unp' in row) and ('user_model' in row):
      print("....Including user_model %s"%row['user_model'])
    ui_final = {}
    for f in ['gene','refseq','mutation','unp']:
      ui_final[f] = row[f]
   
    df_all_jobs,workplan_filename,df_structures,df_dropped,log_filename = plan_one_mutation(row['gene'],row['refseq'],row['mutation'],row['user_model'] if 'user_model' in row and row['user_model'] else None,unp = row['unp'] if 'unp' in row else None)

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
  logger.info("%s",final_structure_info_table)

  # It is so easy to forget to create this phenotypes file - so remind user again!
  if not os.path.exists(phenotypes_filename):
    logger.critical("Reminder: File %s was not created from the UDN report."%phenotypes_filename)
    logger.critical("          Thus, psb_genedicts.py will NOT be part of the planned casewide work")
