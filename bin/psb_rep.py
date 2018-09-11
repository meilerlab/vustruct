#!/usr/bin/env python2.7
#
# Project        : PSB Pipeline
# Filename       : psb_monitor.py
# Authors        : Chris Moth and R. Michael Sivley
#project Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : chris.moth@vanderbilt.edu
# Date           : 2018-01-22
# Description    : Output final html and pdf reports from completed jobs
#                : (output from psb_launch.py)
#=============================================================================#

"""\
Create final html and pdf reports for completed pipeline jobs from either
a project id (ex. UDN123456) or a single workstatus.csv file output from psb_monitor.py

Configuration options must be provided with
   -c global.config file
   -u user.config overrides
   structure_report.csv originally output by psb_plan.py
   workstatus.csv as previously updated with psb_monitor.py
"""

print "4 of 4  %s Pipeline report generator.  -h for detailed help"%__file__

import logging,os,pwd,sys,grp,stat,string
import time, datetime
import argparse,ConfigParser
import pprint
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

rootdir_log_filename = "psb_rep.root.log"
needRoll = os.path.isfile(rootdir_log_filename)
rootdir_fh = RotatingFileHandler(rootdir_log_filename, backupCount=7)
formatter = logging.Formatter('%(asctime)s %(levelname)-4s [%(filename)20s:%(lineno)d] %(message)s',datefmt="%H:%M:%S")
rootdir_fh.setFormatter(formatter)
rootdir_fh.setLevel(logging.INFO)
logger.addHandler(rootdir_fh)

if needRoll:
  rootdir_fh.doRollover()

sys.stderr.write("Root (case) directory log file is %s\n"%rootdir_log_filename)


import shutil
import numpy as np
import pandas as pd
from jinja2 import Environment, FileSystemLoader
from StringIO import StringIO
import pycurl
import json
from xml2json import Cxml2json

from subprocess import Popen, PIPE
from weasyprint import HTML,CSS

default_global_config=os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))),"config","global.config")

cmdline_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
cmdline_parser.add_argument("-c","--config",
help="PDBMap configuration profile for database access\n(default: %(default)s)", required=False,metavar="FILE",default=default_global_config)
cmdline_parser.add_argument("-u","--userconfig",
help="User specific settings and configuration profile overrides", required=True,metavar="FILE")
cmdline_parser.add_argument("-v","--verbose",
help="Include routine info log entries on stderr", action = "store_true")
cmdline_parser.add_argument("-s","--slurm",
help="Create a slurm file to launch all the reports", action = "store_true")
cmdline_parser.add_argument("-d","--debug",
help="Include routine info AND 'debug' log entries on stderr", action = "store_true")
cmdline_parser.add_argument("projectORstructures",type=str,help="The project ID UDN123456 to report on all mutations.  Else, a specific structures file from psb_plan.py  Example: ....../UDN/UDN123456/GeneName_NM_12345.1_S123A_structure_report.csv")
args,remaining_argv = cmdline_parser.parse_known_args()
cmdline_parser.add_argument("workstatus",nargs='?',type=str,help="Blank for all mutations, or specific output file from psb_monitor.py  Example: ....../UDN/UDN123456/GeneName_NM_12345.1_S123A_workstatus.csv")
args,remaining_argv = cmdline_parser.parse_known_args()

# A period in the argument means the user wants to monitor one mutation only,
# directly from a single mutation output file of psb_launch.py
oneMutationOnly = ('.' in args.projectORstructures) and args.workstatus != None
infoLogging = False

if args.debug:
  infoLogging = True
  sh.setLevel(logging.DEBUG)
elif args.verbose:
  infoLogging = True
  sh.setLevel(logging.INFO)

# Load the structure details of all jobs that were not dropped

# print "Command Line Arguments"
logger.info("Command Line Arguments:\n%s"%pprint.pformat(vars(args)))

if not os.path.exists(args.config):
  logger.critical("Global config file not found: " + args.config)
  sys,exit(1)

config = ConfigParser.SafeConfigParser(allow_no_value=True)
configFilesRead = config.read([args.config,args.userconfig])
if len(configFilesRead) == 0:
  logger.critical("Unable to open config files: %s or %s  Exiting"%(args.config,args.userconfig))
  sys.exit(1)
logger.info("Successful read of config files: %s"%str(configFilesRead))
config_dict = dict(config.items("Genome_PDB_Mapper")) # item() returns a list of (name, value) pairs
# Init all the user config dictionaries to empty
SlurmParametersReport = {}  
userSlurmParametersAll = userSlurmParametersReport = {}

if args.slurm:
  try:
    globalSlurmParametersAll = dict(config.items("SlurmParametersAll"))
  except Exception as ex:
    globalSlurmParametersAll = {}

  try:
    globalSlurmParametersReport = dict(config.items("SlurmParametersReport"))
  except Exception as ex:
    globalSlurmParametersReport = {}

if args.userconfig:
  userconfig = ConfigParser.SafeConfigParser() 
  userconfig.read([args.userconfig])
  config_dict.update(dict(userconfig.items("UserSpecific"))) # item() returns a list of (name, value) pairs
  
  try:
    userSlurmParametersReport = dict(userconfig.items("SlurmParametersReport"))
  except:
    pass

if args.slurm:
  SlurmParametersReport.update(globalSlurmParametersAll)
  SlurmParametersReport.update(userSlurmParametersAll)
  SlurmParametersReport.update(globalSlurmParametersReport)
  SlurmParametersReport.update(userSlurmParametersReport)




# The collaboration_dir is the master directory for the case, i.e. for one patient
# Example: /dors/capra_lab/projects/psb_collab/UDN/UDN532183
udn_root_directory = os.path.join(config_dict['output_rootdir'],config_dict['collaboration'])
if oneMutationOnly:
  collaboration_dir = os.path.dirname(os.path.dirname(args.projectORstructures))
else:
  collaboration_dir = os.path.join(udn_root_directory,args.projectORstructures)

# Return a dictionary of "generic Interactions" that can flow into the summary report
# Return an html table that can be placed at end of summary report
# Return html long report
# Return text that can be dumped to the log file
def GeneInteractionReport(case_root,case,CheckInheritance):
    genepair_dict_file = os.path.join(case_root,'casewide','MakeGeneDictionaries','%s.json'%case)
    try:
        with open(genepair_dict_file,'r') as fd:
    	    genepair_dict = json.load(fd)
    except:
        logger.exception("Cannot open %s.  Did you run Souhrid's analysis scripts?"%genepair_dict_file)
        return None,None,None,None
       
    genes_filename = os.path.join(case_root,"%s_genes.txt"%case)
    try:
        with open(genes_filename) as fd:
    	      unstripped_genes = fd.readlines()
    except:
        logger.exception("Cannot open %s.  Did you run Souhrid's analysis scripts?"%genes_filename)
        return None,None,None,None

    geneInteractions = {}
    html_table =  StringIO()
    html_report = StringIO()
    text_report = StringIO()

    structure_positive_genes = [gene.strip() for gene in unstripped_genes]
    
    pairs = {}
    pairs['direct'] = []
    pairs['pathways'] = []
    pairs['phenotypes'] = []
    pairs['other'] = []

    # We let the template html open the <table> with formatting as it prefers
    # print >>html_table, "<table>" 
    print >>html_table, "<thead><tr>" 
    print >>html_table, "<th></th>"# top left is blank cell in header
    for gene1 in structure_positive_genes:
      print >>html_table, "<th>" + gene1 + "</th>" 
    print >>html_table, "</tr></thead>" 

    for gene1 in structure_positive_genes:
        print >>text_report, "%s:"%gene1
        print >>html_report, '%s:'%gene1
        print >>html_table, "<tr><td>%s</td>"%gene1 
        if gene1 in genepair_dict:
          print >>html_report, '<ul style="list-style: none;">'
        for gene2 in structure_positive_genes:
          print >>html_table, "<td>",
          if gene1 in genepair_dict and gene2 in genepair_dict[gene1] and ((not CheckInheritance) or (genepair_dict[gene1][gene2]['Inheritance'])):
              pairs_hits = []
              if genepair_dict[gene1][gene2]['Direct_Interaction']:
                  pairs_hits.append('direct')
              if genepair_dict[gene1][gene2]['Common_pathways']:
                  pairs_hits.append('pathways')
              if genepair_dict[gene1][gene2]['Common_phenotypes'] and genepair_dict[gene1][gene2]['Patient_Phenotype_Overlap'] > 0:
                  pairs_hits.append('phenotypes')
              if genepair_dict[gene1][gene2]['Distance'] < 3 and genepair_dict[gene1][gene2]['Coexpression'] < 50:
                  pairs_hits.append('other')
              if pairs_hits:
                  for pairs_hit in pairs_hits:
                      pairs[pairs_hit].append((gene1,gene2))
                      if pairs_hit != 'phenotypes': # We don't get too excited for phenotypes at this point
                          if not gene1 in geneInteractions:
                              geneInteractions[gene1] = []
                          if not gene2 in geneInteractions[gene1]:
                              geneInteractions[gene1].append(gene2)
                  logging.debug("Gene Interactions: %s %s %s"%(gene1, gene2, str(pairs_hits)))
                  print >>text_report, gene2, str(pairs_hits) 
                  print >>html_report, "<li>",gene2 + ':', ", ".join(pairs_hits),"</li>"
                  print >>html_table, '<br>'.join(pairs_hits) 
          print >>html_table, "</td>" 
        print >>html_report, '</ul>'
        print >>html_table, "</tr>"
    # print >>html_table, "</table>" 
    return geneInteractions,html_table.getvalue(),html_report.getvalue(),text_report.getvalue()

# This complex routine gather PathProx, and ddG results from the pipeline run, where available
# It creates a sigle .html file representing the mutation, and returns a dictionary that is
# input to the overall case report in gathered_info
def report_one_mutation(structure_report,workstatus):
  gathered_info = {} # Dictionary info to return, of great use by caller
  # Load the status of jobs that was created by psb_launch.py (or prior run of psb_monitor.py)
  df_all_jobs_status = pd.read_csv(workstatus,'\t',keep_default_na = False,na_filter=False,dtype=str)
  msg = "%d rows read from work status file %s"%(len(df_all_jobs_status),workstatus)
  if len(df_all_jobs_status) < 1:
    logger.critical(msg)
    return None;
  else:
    if not infoLogging:
      print msg
    logger.info(msg)

  row0 = df_all_jobs_status.iloc[0]
  gathered_info = {col: row0[col] for col in ['project','unp','gene','refseq','mutation']}
  template_vars = gathered_info.copy()

  mutation = row0['mutation'] # Sorry for alias - but this is all over the code

  gathered_info['Error'] = '' # < This is False for truth but looks fine on the report
  gathered_info['ddG Max'] = gathered_info['ddG Min'] = None
  gathered_info['csm_pp'] = gathered_info['cln_pp'] = None


  # import pdb; pdb.set_trace()
  
  udn_sequence_annotation_jobs = ((df_all_jobs_status['flavor'] == 'SequenceAnnotation'))
  if udn_sequence_annotation_jobs.sum() > 0:
    logger.info("Dropping %d rows where jobs are UDN SequenceAnnotation"%udn_sequence_annotation_jobs.sum())
    df_all_jobs_status = df_all_jobs_status[~udn_sequence_annotation_jobs]


  # ExitCode column arrives when monitor program finds a clear exit, but may not arrive
  # if programs exit without leaving a trail of info 
  if not 'ExitCode' in df_all_jobs_status.columns: 
    gathered_info['Error'] = "'ExitCode' column missing."
    logger.error(gathered_info['Error'] + "in " + workstatus)
    return gathered_info

  failed_jobs = df_all_jobs_status['ExitCode'] > '0'
  if failed_jobs.sum() > 0:
    logger.info("Dropping %d rows for jobs that failed (ExitCode > '0')."%failed_jobs.sum())
    df_all_jobs_status = df_all_jobs_status[~failed_jobs]

  incomplete_jobs = df_all_jobs_status['ExitCode'] != '0'
  if incomplete_jobs.sum() > 0:
    logger.info("Dropping %d rows for incomplete jobs (ExitCode missing)."%incomplete_jobs.sum())
    df_all_jobs_status = df_all_jobs_status[~incomplete_jobs]
    if len(df_all_jobs_status) < 1:
      gathered_info['Error'] = "All jobs are marked incomplete."
      logger.error(gathered_info['Error'] + "  Be sure to run pdb_monitr.py before this script")
      return gathered_info

  if len(df_all_jobs_status) < 1:
      gathered_info['Error'] = "No completed PathProx or ddG jobs."
      logger.warning(gathered_info['Error'] + " in " + workstatus)
      return gathered_info

  # Sanity check - make sure mutation is same in all rows!
  test_count = len(df_all_jobs_status.groupby('mutation'))
  if test_count != 1:
    gathered_info['Error'] = "%s  may be coreupted, has some rows which lack the same mutation"%workstatus
    logger.error(gathered_info['Error'])
    gathered_info['Error'] = "corrupted = see log file"
    return gathered_info

  # Sanity check - make sure refseq is same in all rows!
  assert len(df_all_jobs_status.groupby('refseq')) == 1
  
  # Sanity check - make sure gene is same in all rows!
  assert len(df_all_jobs_status.groupby('gene')) == 1




  mutation_dir = os.path.join(collaboration_dir,"%s_%s_%s"%(gathered_info['gene'],gathered_info['refseq'],mutation))

  if not os.path.exists(mutation_dir):  # python 3 has exist_ok parameter...
    os.makedirs(mutation_dir)
  mutation_log_dir = mutation_dir # os.path.join(mutation_dir,"log")
  if not os.path.exists(mutation_log_dir):  # python 3 has exist_ok parameter...
    errorMsg = "%s log directory should have been created by psb_plan.py. Terminating"%mutation_log_dir
    logger.critical(errorMsg)
    sys.stderr.write(errorMsg)
    sys.exit(1)

  # pw_name = pwd.getpwuid( os.getuid() ).pw_name # example jsheehaj or mothcw
  log_filename = os.path.join(mutation_log_dir,"psb_rep.log")

  if oneMutationOnly:
    sys.stderr.write("psb_rep log file is %s\n"%log_filename)
  else:
    logger.info("additionally logging to file %s"%log_filename)
  needRoll = os.path.isfile(log_filename)
  local_fh = RotatingFileHandler(log_filename, backupCount=7)
  formatter = logging.Formatter('%(asctime)s %(levelname)-4s [%(filename)20s:%(lineno)d] %(message)s',datefmt="%H:%M:%S")
  local_fh.setFormatter(formatter)
  local_fh.setLevel(logging.INFO)
  logger.addHandler(local_fh)

  if needRoll:
    local_fh.doRollover()
# Load the structure details of all jobs that were not dropped

# Load the structure details of all jobs that were not dropped
  nan_values = ['-1.#IND', '1.#QNAN', '1.#IND', '-1.#QNAN', '#N/A','N/A', '#NA', 'NULL', 'NaN', '-NaN', 'nan', '-nan']
  df_structure_report = pd.read_csv(structure_report,'\t',keep_default_na=False,na_values = nan_values,dtype={'Resolution': float, 'Seq Identity': float})
  # import pdb; pdb.set_trace()
  logger.info("%d rows read from structure_report file %s"%(len(df_structure_report),structure_report))
  df_structure_report.set_index(["Method","Structure ID","Chain"],inplace = True)
  
  # For final report directory, trim that last directory from the structure-specific one
  final_report_directory = df_all_jobs_status.iloc[0]['outdir']
  # Now trim final directory name
  final_report_directory = os.path.dirname(final_report_directory)
  
  # project = 'TestJED' # MUST FIX THIS ERROR
  # unp = 'Q9NY15-1' # MUST FIX THIS ERROR
  # import pdb; pdb.set_trace() 
  df_all_jobs_status.set_index('uniquekey',inplace=True)
  
  df_complete_jobs = df_all_jobs_status[df_all_jobs_status['ExitCode'] == '0'].copy()
  
  
  print "Generating html and pdf final reports for %d of %d complete jobs:"%(len(df_complete_jobs),len(df_all_jobs_status))
  
  def structure_html_vars(df_onerow,df_structure_report_onerow,thestruct):
    # Subset the dataframe to the results of this chain
    # This could need 2 rows - one for each pathprox - we just have to keep an eye onthis
    structure  = {"structid":thestruct['pdbid'],
  		"chain":thestruct['chain'],
  		"html_id":(thestruct['pdbid'] + thestruct['chain']).translate(None, string. punctuation),
  		"mutation":gathered_info['mutation']}

  
    structure.update(df_structure_report_onerow.to_dict())
  
    # one_job_status = df_onerow.to_frame()
    # one_structure_report = structure_report_onerow.toframe()
  
    # Make a one row dataframe just to convert to html
    # tempdf = df_structure_report_onerow
    struct_detail_columns = (["Analyzable?","Distance to Boundary","Seq Identity","PDB Pos","PDB Template","Residues",
             "Resolution (PDB)","Seq Start","Seq End","Transcript Pos"])
    dfT = df_structure_report_onerow.to_frame().transpose()[struct_detail_columns]
    # import pdb;pdb.set_trace()
    if df_structure_report_onerow['Label'] == 'pdb':
      dfT['Seq Identity'] = 100.0
    # dfT[['Seq Identity']] = dfT[['Seq Identity']].astype('float')
   
    structure["details"] = dfT.to_html(float_format=lambda x: "%.2f"%x)
    return structure
  
  def ddG_html_vars(structure,thestruct):
    if not 'ddG' in thestruct:
      return structure,None

    df_row = thestruct['ddG']
    output_flavor_directory = os.path.join(df_row['outdir'],'ddG','ddg')

    ddGPrefix = "%s_%s_%s"%(thestruct['pdbid'],thestruct['chain'],thestruct['mutation'])
    ddG_summary_file =  "%s/%s_ddg.results"%(output_flavor_directory,ddGPrefix)
    if os.path.exists(ddG_summary_file):
      summary = pd.read_csv(ddG_summary_file,sep='\t')
      logger.debug('We got the file %s'%ddG_summary_file)
      structure['ddG'] = summary['ddG']
      return structure,summary

    logger.debug('Unable to open %s'%ddG_summary_file)

    return structure,None

  
  # Compile PathProx results into report (one page per structure+chain+mutation)
  def PathProx_html_vars(structure,thestruct,Clinvar_or_COSMIC):
    """ Generate the HTML for a single page structure report """
  
    # If there is no thestruct['Clinvar'] or thestruct['COSMIC'] then return right away - no big deal
    if not Clinvar_or_COSMIC in thestruct:
      return structure,None
  
  
    if (Clinvar_or_COSMIC == 'Clinvar'):
      clinvar_or_cosmic = 'clinvar'
      cln_or_csm = 'cln'
    elif (Clinvar_or_COSMIC == 'COSMIC'):
      clinvar_or_cosmic = 'cosmic'
      cln_or_csm = 'csm'
    else:
      assert Clinvar_or_COSMIC == 'Clinvar' or Clinvar_or_COSMIC == 'COSMIC'
  
    df_row = thestruct[Clinvar_or_COSMIC]
  
    # Define directory prefixes
    # This is where we'll search for those final pathprox .pdb files, etc 
    output_flavor_directory = os.path.join(df_row['outdir'],df_row['flavor'])
  
    PathProxPrefix = "%s_%s_%s_exac_D"%(thestruct['pdbid'],thestruct['chain'],clinvar_or_cosmic)
  
    # Rename the primary columns for each results set
    coldir = {"Kz_path":"Pathogenic Kz","Kt_path":"Pathogenic Kt",
                 "Kz_neut":"Neutral Kz","Kt_neut":"Neutral Kt",
                 "roc_auc":"ROC AUC","pr_auc":"PR AUC",
                 "%s_pathprox"%mutation:"PathProx",
                 "%s_neutcon"%mutation:"Neutral Constraint",
                 "%s_pathcon"%mutation:"Pathogenic Constraint"}
    colord = ["Neutral Kz","Neutral Kt","Pathogenic Kz","Pathogenic Kt",
                "ROC AUC","PR AUC","Mutation","PathProx","Neutral Constraint",
                "Pathogenic Constraint"]
    # Load the Pathprox summary of results for this mutation
    pathprox_summary_file =  "%s/%s_summary.csv"%(output_flavor_directory,PathProxPrefix)
    if os.path.exists(pathprox_summary_file):
      # For this read below, we expect occasional meaningful Nans
      summary = pd.read_csv(pathprox_summary_file,sep='\t',header=0)
      summary["Mutation"] = gathered_info['mutation']
      # Dummy code coverage-dependent columns if missing
      if "%s_pathprox"%mutation not in summary.columns:
        summary["%s_pathprox"%mutation] = np.nan
      if "%s_neutcon"%mutation not in summary.columns:
        summary["%s_neutcon"%mutation] = np.nan
      if "%s_pathcon"%mutation not in summary.columns:
        summary["%s_pathcon"%mutation] = np.nan
      summary = summary.rename(columns=coldir)
      summary = summary[colord]
      structure["%s_results"%cln_or_csm] = summary.to_html(index=False,float_format=lambda x: "%.2f"%x)
    else:
      print "PathProx summary file not found: %s"%pathprox_summary_file
      structure["%s_results"%cln_or_csm] = ""
      summary = pd.DataFrame()
    """
    # Load the COSMIC summary of results for this mutation
    if os.path.exists("%s/%s_summary.txt"%(dors_ppdir,csm_prefix)):
      csm_sum = pd.read_csv("%s/%s_summary.txt"%(dors_ppdir,csm_prefix),sep='\t',header=0)
      csm_sum["Mutation"] = mutation
      # Dummy code coverage-dependent columns if missing
      if "%s_pathprox"%mutation not in csm_sum.columns:
        csm_sum["%s_pathprox"%mutation] = np.nan
      if "%s_neutcon"%mutation not in csm_sum.columns:
        csm_sum["%s_neutcon"%mutation] = np.nan
      if "%s_pathcon"%mutation not in csm_sum.columns:
        csm_sum["%s_pathcon"%mutation] = np.nan
      csm_sum = csm_sum.rename(columns=coldir)
      csm_sum = csm_sum[colord]
      structure["csm_results"] = csm_sum.to_html(index=False,float_format=lambda x: "%.2f"%x)
    else:
      print "COSMIC PathProx summary file not found: %s/%s_summary.txt"%(dors_ppdir,csm_prefix)
      structure["csm_results"] = ""
      csm_sum = pd.DataFrame()
    """
    # Load the variant plots
    # Use  root-relative pathnames as inputs to the structure[] diction that flows into the html creation
    # vus_var and exac_var should be populated from either Clinvar or COSMIC calculations.  Get them we don't
    # have them yet
    if not "vus_var" in structure and os.path.exists("%s/%s_structure.png"%(output_flavor_directory,PathProxPrefix)):
      structure["vus_var"]    = "%s/%s_structure.png"%(output_flavor_directory,PathProxPrefix)
  
    if not "exac_var" in structure and os.path.exists("%s/%s_neutral.png"%(output_flavor_directory,PathProxPrefix)):
      structure["exac_var"] = "%s/%s_neutral.png"%(output_flavor_directory,PathProxPrefix)
 
    cln_or_csmPathogenicPNGfile = "%s/%s_pathogenic.png"%(output_flavor_directory,PathProxPrefix)
    if os.path.exists(cln_or_csmPathogenicPNGfile):
      structure["%s_var"%cln_or_csm]    = cln_or_csmPathogenicPNGfile
    else:
      structure["%s_var"%cln_or_csm]    = 0

    # Where pathprox has left us ngl viewer variant setups, incorporate those into the final report
    # As with the static .png im

    def residue_list_to_ngl(residues):
      if len(residues) == 0:
        return None
      retval = '('
      for r in residues:
        retval += "%d "%r[0]
      return retval + ')'

    def residue_list_to_ngl_CAs(residues):
      if len(residues) == 0:
        return None
      return residue_list_to_ngl(residues) + "and .CA"
      
    json_filename = os.path.join(output_flavor_directory,PathProxPrefix + "_ResiduesOfInterest.json")
    if os.path.exists(json_filename):
      with open(json_filename) as f:
        residuesOfInterest = json.load(f)
        logger.info('Loaded %d variants, %d neutrals, %d pathogenic from %s'%(
          len(residuesOfInterest['variants']),
          len(residuesOfInterest['neutrals']),
          len(residuesOfInterest['pathogenics']),
          json_filename))

      # For final .html - we need to be specific regarding type of pathogenic.  We weren't when json file was written
      residuesOfInterest['%s_pathogenics'%cln_or_csm] = residuesOfInterest['pathogenics']
      del residuesOfInterest['pathogenics']
      structure.update(residuesOfInterest)
      structure['ngl_variant_residue_count'] = len(residuesOfInterest['variants'])
      structure['ngl_variant_residues'] = residue_list_to_ngl(residuesOfInterest['variants'])
      structure['ngl_neutral_residue_count'] = len(residuesOfInterest['neutrals'])
      structure['ngl_neutral_residues'] = residue_list_to_ngl_CAs(residuesOfInterest['neutrals'])
      structure['ngl_%s_residue_count'%cln_or_csm] = len(residuesOfInterest['%s_pathogenics'%cln_or_csm])
      structure['ngl_%s_residues'%cln_or_csm] = residue_list_to_ngl_CAs(residuesOfInterest['%s_pathogenics'%cln_or_csm])

    
    # if not "vus_ngl_html" in structure and os.path.exists("%s/%s_structure.ngl.html"%(output_flavor_directory,PathProxPrefix)):
      # with open("%s/%s_structure.ngl.html"%(output_flavor_directory,PathProxPrefix),"rb") as ngl_html_file:
        # structure["vus_ngl_html"] =  ngl_html_file.read()
 
    if not "exac_ngl_html" in structure:
      html_filename =  "%s/%s_neutral.ngl.html"%(output_flavor_directory,PathProxPrefix)
      if os.path.exists(html_filename):
        with open("%s/%s_neutral.ngl.html"%(output_flavor_directory,PathProxPrefix),"rb") as ngl_html_file:
          structure["exac_ngl_html"] = ngl_html_file.read()

    if os.path.exists("%s/%s_pathogenic.ngl.html"%(output_flavor_directory,PathProxPrefix)):
      with open("%s/%s_pathogenic.ngl.html"%(output_flavor_directory,PathProxPrefix),"rb") as ngl_html_file:
        structure["%s_ngl_html"%cln_or_csm]    = ngl_html_file.read()

    # Load the Ripley's K results for ClinVar and COSMIC

    # If an image file is in the filesystem, great, else return empty string  
    def checkf(plot_png):
      if os.path.exists(plot_png):
        return plot_png
      return ""

    # Don't re-do exac_K graphic if already loaded from the other directory
    if (not "exac_K" in structure) or len(structure["exac_K"]) == 0:
      structure["exac_K"]   = checkf("%s/%s_neutral_K_plot.png"%(output_flavor_directory,PathProxPrefix))

    structure["%s_K"%cln_or_csm]      = checkf("%s/%s_pathogenic_K_plot.png"%(output_flavor_directory,PathProxPrefix))
    # Load the Ripley's D results for ClinVar and COSMIC
    structure["%s_D"%cln_or_csm]      = checkf("%s/%s_D_plot.png"%(output_flavor_directory,PathProxPrefix))
    # Load the PathProx Mapping and Performance for ClinVar
    structure["%s_pp"%cln_or_csm]     = checkf("%s/%s_pathprox.png"%(output_flavor_directory,PathProxPrefix))
    structure["%s_roc"%cln_or_csm]    = checkf("%s/%s_pathprox_roc.png"%(output_flavor_directory,PathProxPrefix))
    structure["%s_pr"%cln_or_csm]     = checkf("%s/%s_pathprox_pr.png"%(output_flavor_directory,PathProxPrefix))
    
    # Return the structure and the ClinVar/COSMIC results
    return structure,summary
  
  def unp2PfamDomainGraphicString(unp,timeoutSeconds):
    url = 'http://pfam.xfam.org/protein/{0}/graphic'.format(unp.split(' ')[0].split('-')[0])
  
    buffer = StringIO()
  
    try:
      timeout = 60 #NOTE: using pycurl assures that this timeout is honored
      c = pycurl.Curl()
      c.setopt(pycurl.CONNECTTIMEOUT, timeoutSeconds)
      c.setopt(pycurl.TIMEOUT, timeoutSeconds)
      c.setopt(pycurl.WRITEFUNCTION, lambda x: None)
      c.setopt(pycurl.WRITEFUNCTION, buffer.write)
      c.setopt(pycurl.HEADERFUNCTION, lambda x: None)
      c.setopt(pycurl.NOSIGNAL, 1)
      c.setopt(pycurl.URL, url)
      c.setopt(pycurl.HTTPGET, 1)
      graphicJSON = str(c.perform())
      c.close()
  
    except Exception as ex:
     logger.exception('Failed to get JSON string from pfam.xfam.org')
     return None
  
    JSONstr= buffer.getvalue()
    if len(JSONstr) > 3 and JSONstr[0] == '[':
      JSONstr = JSONstr[1:-1]
    return JSONstr
  
  
  df_complete_jobs['Structure ID'] = df_complete_jobs['pdbid']
  df_complete_jobs['Method'] = df_complete_jobs['method']
  df_complete_jobs['Chain'] = df_complete_jobs['chain']
  
  tdf        = df_complete_jobs.set_index(["Method","Structure ID","Chain"])
  tdf_unique_structures = tdf[~tdf.index.duplicated(keep='first')]
  duplicate_structure_indices = ~tdf.index.duplicated(keep='first')
  tdf_unique_structures = tdf[duplicate_structure_indices]
  
 
   
  # The final .html (.pdf) output is organized structure by structure.
  # To stay sane, we somewhat need a dictionary with one element per structure
  # Here, we'll have struct_dict[Method/pdbid/chain] with sub hash elements:
  # mutation/pdbid/chain/
  # COSMIC dict/Clinvar dict/ddg Dict are dictionaries which have the entire rows
  # From the completed job status
  struct_dict = {}
  for index,row in df_complete_jobs.iterrows():
    key_tuple = (row['Method'],row['pdbid'],row['chain'])
    thestruct = struct_dict.get(key_tuple,{})
    thestruct['mutation'] = row['mutation']
    thestruct['pdbid'] = row['pdbid']
    thestruct['chain'] = row['chain']
    if 'COSMIC' in row['flavor']:
      thestruct['COSMIC'] = row.to_dict()
    elif 'Clinvar' in row['flavor']:
      thestruct['Clinvar'] = row.to_dict()
    elif 'SequenceAnnotation' == row['flavor']: 
      continue # UDN Sequence annotations are NOT a part of generated reports
    elif 'ddG' == row['flavor']: 
      thestruct['ddG'] = row.to_dict()
    else:
      logger.critical("flavor in row is neither COSMIC nor Clinvar nor ddG - cannot continue:\n%s",str(row))
      sys.exit(1)
    # This will often reassign over prior assignments - That's OK
    struct_dict[key_tuple] = thestruct
  
  print "%d structures have at least one complete report"%len(struct_dict)
  assert len(struct_dict) == len(tdf_unique_structures)
  print "*** THINK ABOUT THIS ASSERT BELOW ***"
  # assert len(struct_dict) == len(df_structure_report)
  
  # Generate vars for individual structure pages
  print "\nLoading the results from each structural analysis..."
  
  structures = []
  cln_list = []
  csm_list = []
  ddG_list = []
  # import pdb; pdb.set_trace()
  # Load up variables for each structure in our set
  for key_tuple in struct_dict:
    # import pdb; pdb.set_trace()
    s = structure_html_vars(tdf_unique_structures.loc[key_tuple],df_structure_report.loc[key_tuple],struct_dict[key_tuple])
    s,Clinvar_sum = PathProx_html_vars(s,struct_dict[key_tuple],'Clinvar')
    if Clinvar_sum is not None:
      cln_list.append(Clinvar_sum)
    s,COSMIC_sum = PathProx_html_vars(s,struct_dict[key_tuple],'COSMIC')
    if COSMIC_sum is not None:
      csm_list.append(COSMIC_sum)
    s,ddG_summary = ddG_html_vars(s,struct_dict[key_tuple])
    if ddG_summary is not None:
      ddG_list.append(ddG_summary)
 
    structures.append(s)
  if cln_list:
    cln_sum = pd.concat(cln_list)
  else:
    cln_sum = pd.DataFrame() 
  if csm_list:
    csm_sum = pd.concat(csm_list)
  else:
    csm_sum = pd.DataFrame() 
  if ddG_list:
    ddG_sum = pd.concat(ddG_list)
    gathered_info['ddG Max'] = ddG_sum['ddG'].max()
    gathered_info['ddG Min'] = ddG_sum['ddG'].min()
  else:
    ddG_sum = pd.DataFrame() 
    gathered_info['ddG Max'] = None
    gathered_info['ddG Min'] = None
    
  # Summarize the results for this gene mutation
  aggcols = ["Mutation","Neutral Kz","Neutral Kt","Pathogenic Kz","Pathogenic Kt","ROC AUC",
              "PR AUC","PathProx","Neutral Constraint","Pathogenic Constraint"]
  if not cln_sum.empty:
    msg = "Summarizing results from %d ClinVar PathProx analyses..."%len(cln_sum)
    if not infoLogging:
      print msg
    logger.info(msg)
    # Fancy way to compete the average for all the aggcols 
    cln_sum = cln_sum.groupby("Mutation")[aggcols].aggregate(np.nanmedian)
    cln_auc = cln_sum["ROC AUC"].values[0]
    gathered_info['cln_pp']  = cln_sum["PathProx"].values[0]
    cln_auc = cln_auc if not np.isnan(cln_auc) else None
    gathered_info['cln_pp']  = gathered_info['cln_pp'] if not np.isnan(gathered_info['cln_pp']) else None
  else:
    cln_auc = gathered_info['cln_pp'] = None
  
  if not csm_sum.empty:
    msg = "Summarizing results from %d COSMIC PathProx analyses..."%len(csm_sum)
    if not infoLogging:
      print msg
    logger.info(msg)
    # Fancy way to compete the average for all the aggcols 
    csm_sum = csm_sum.groupby("Mutation")[aggcols].aggregate(np.nanmedian)
    csm_auc = csm_sum["ROC AUC"].values[0]
    gathered_info['csm_pp']  = csm_sum["PathProx"].values[0]
    csm_auc = csm_auc if not np.isnan(csm_auc) else None
    gathered_info['csm_pp']  = gathered_info['csm_pp'] if not np.isnan(gathered_info['csm_pp']) else None
  else:
    csm_auc = gathered_info['csm_pp'] = None


  if not ddG_sum.empty:
    msg = "Summarizing results from %d ddg runs..."%len(ddG_sum)
    if not infoLogging:
      print msg
    logger.info(msg)

  # import pdb;pdb.set_trace()
  logger.warn("REminder: DATABASE IS NOT BEING UPDATED")
 
  """ 
  # Upload results to the database
  from pdbmap import PDBMapIO
  io = PDBMapIO(config_dict['dbhost'],config_dict['dbuser'],config_dict['dbpass'],config_dict['dbname']) # ,slabel=args.slabel)
  c = io._con.cursor()
  sql  = "UPDATE psb_collab.MutationSummary SET "
  sql += "clinvar_roc_med=%s,  cosmic_roc_med=%s,"
  sql += "clinvar_pathprox=%s, cosmic_pathprox=%s "
  sql += "WHERE collab=%s AND project=%s AND gene=%s AND refseq=%s AND mutation=%s"
  c.execute(sql,(cln_auc,csm_auc,gathered_info['gathered_info['cln_pp']'],gathered_info['csm_pp'],
            config_dict['collaboration'],gathered_info['project'],gathered_info['gene'],gathered_info['refseq'],mutation)
  c.close()
  """
  
  # Generate vars for overall html template
  pdb_html = "No structures or homology models were identified for this gene."
  cln_html = "Insufficient variation was available for PathProx analysis."
  csm_html = cln_html
  ddG_html = "No ddG calculations were completed for this mutation."
 
  columns_for_html = ['Analyzable?','Label','Distance to Boundary','Seq Start','Seq End','Seq Identity','PDB Pos','PDB Template','Residues','Resolution (PDB)','Transcript Pos']
  
  if not tdf.empty:
    # pdb_html = tdf.to_html(float_format=lambda x: "%.2f"%x,formatters={'Seq Start': '${:,1.0f}'.format,'Seq End': '${:,1.0f}'.format})
    df_structure_report.style.set_properties(**{'text-align': 'center'})
    # dfTemp = df_structure_report[columns_for_html].copy()
    # dfTemp[['Seq Identity']] = dfTemp[['Seq Identity']].astype('float')
    # import pdb; pdb.set_trace()
    df_copy =  df_structure_report[columns_for_html].copy() # .to_html(justify='center') # ,formatters={'Seq Identity': '{:,.2f}'.format})
    # pd.set_option('float_format', '%.2f')
    pdb_html = df_copy.to_html(justify='center',float_format=lambda x: "%.2f"%x) 
  if not cln_sum.empty:
    cln_html = cln_sum.to_html(float_format=lambda x: "%.2f"%x)
  if not csm_sum.empty:
    csm_html = csm_sum.to_html(float_format=lambda x: "%.2f"%x)
  if not ddG_sum.empty:
    ddG_html = ddG_sum.to_html(float_format=lambda x: "%.2f"%x)
  
  
  graphicsLegend = "Pfam Domain Graphic unavailable"
  # Get Uniprot Graphics string or put in a blank
  domainGraphicsJSONstr = None
  if gathered_info['unp']: # Fetch the JSON - wait up to 20 seconds
    xmlHandle = Cxml2json(config_dict['interpro_dir'] + 'match_humanonly.xml')
    if not xmlHandle.isCanonical(gathered_info['unp']): # Then we need to manually create the xml ourselves
      graphicsLegend = "Pfam Domain Graphic for Non-canonical Isoform %s"%gathered_info['unp']
      nonCanonicalGraphicsJSON = xmlHandle.PFAMgraphicsJSONfromUnp(gathered_info['unp'])
      domainGraphicsJSONstr = json.dumps(nonCanonicalGraphicsJSON)
    else:
      graphicsLegend = "Pfam Domain Graphic for Canonical Isoform %s"%gathered_info['unp']
      domainGraphicsJSONstr = unp2PfamDomainGraphicString(gathered_info['unp'],20)
      # Not sure why - but graphics string from web is a list (inside outer [])
      # Removed this - did not work: domainGraphicsJSONstr = domainGraphicsJSONstr[0]
  
    if (domainGraphicsJSONstr != None):
      domainGraphicsDict = json.loads(domainGraphicsJSONstr)
      # Add our mutation point of interest to this
      mutationSiteDict = {u'colour': u'#e469fe',\
       u'display': True,\
       u'headStyle': u'diamond',\
       u'lineColour': u'#333333',\
       u'metadata': {u'database': u'UDN Case',\
                u'description': u'%s'%mutation,\
                u'start': mutation[1:-1],\
                u'type': u'Mutation'} ,\
       u'residue': u'X',\
       u'start': mutation[1:-1],\
       u'type': u'UDN Mutation site',\
       u'v_align': u'top'}
      markupDictList = domainGraphicsDict.get('markups',None)
      if (markupDictList == None):
        markupDictList = []
      markupDictList.append(mutationSiteDict)
      domainGraphicsDict['markups'] = markupDictList
      for region in domainGraphicsDict.get('regions',[]):
        if 'href' in region and region['href'].startswith('/'):
          region['href'] = 'http://pfam.xfam.org' + region['href']
      domainGraphicsDictList = [domainGraphicsDict]
      domainGraphicsJSONstr = json.dumps(domainGraphicsDictList)
  
  # Build a table to explaine the pfam graphics
  PfamResultTableHTML = '''<table class="resultTable details"
           id="imageKey"
           summary="Key for the Pfam domain image">
      <thead>
        <tr>
          <th class="dh" rowspan="2">Source</th>
          <th class="dh" rowspan="2">Domain</th>
          <th class="dh" rowspan="2">Start</th>
          <th class="dh" rowspan="2">End</th>
          <th class="sh" style="display: none" colspan="2">Gathering threshold (bits)</th>
          <th class="sh" style="display: none" colspan="2">Score (bits)</th>
          <th class="sh" style="display: none" colspan="2">E-value</th>
        </tr>
        <tr>
          <th class="sh" style="display: none">Sequence</th>
          <th class="sh" style="display: none">Domain</th>
          <th class="sh" style="display: none">Sequence</th>
          <th class="sh" style="display: none">Domain</th>
          <th class="sh" style="display: none">Sequence</th>
          <th class="sh" style="display: none">Domain</th>
        </tr>
      </thead>'''
  
  PfamResultTableHTML += '<tbody>'
  
  table_rows = []
  for motif in domainGraphicsDict.get('motifs',[]):
    table_rows.append( {
      'type': motif['metadata']['type'],
      'domain': 'n/a',
      'start': motif['start'],
      'end': motif['end']} )
  
  for region in domainGraphicsDict.get('regions',[]):
    table_rows.append( {
      'type': region['metadata']['type'],
      'domain': "<a href=" + region['href'] + ">" + region['text'] + "</a>",
      'start': region['start'],
      'end': region['end']} )
  
  table_rows.append( {
      'type': 'Mutation',
      'domain': 'n/a',
      'start': mutation[1:-1],
      'end': mutation[1:-1]} )
  
  sorted_table_rows = sorted(table_rows, key=lambda d: int(d['start']))
  
  for sr in sorted_table_rows:
    table_row = '<tr class="odd">'
    table_row += '<td class="domain">%s</td>'%sr['type']
    table_row += '<td><span class="inactive">%s</span></td>'%sr['domain']
    table_row += '<td>%s</td>'%sr['start']
    table_row += '<td>%s</td>'%sr['end']
    table_row += 6* '<td class="sh" style="display: none"><span class="inactive">n/a</span></td>'
    table_row += '</tr>\n'
    PfamResultTableHTML += table_row
  
  
  PfamResultTableHTML +=  '</tbody></table>'
  pfamGraphicsIframe_fname = "html/%s_%s_PfamGraphicsIframe.html"%(gathered_info['gene'],gathered_info['mutation'])
  
 
  # This data structure directly feeds the complex items in the psb_report.html template 
  # elements for project/case/gene/mutation/refseq were populated at top of this function
  template_vars.update({"title" : "PSB Structure Report - %s %s"%(gathered_info['unp'],gathered_info['mutation']),
                   "date"  : time.strftime("%Y-%m-%d"),
                   "pfamGraphicsIframe_fname" : pfamGraphicsIframe_fname,
                   "graphicsLegend" : graphicsLegend,
                   "psb_structure_results" : pdb_html,
                   "clinvar_summary"   : cln_html,
                   "cosmic_summary"    : csm_html,
                   "ddG_summary"    : ddG_html,
                   "DomainGraphicsJSON"  : domainGraphicsJSONstr,
                   "PfamResultTableHTML"  : PfamResultTableHTML,
                   "Structure_Detail"  : structures})
  if mutation:
    base_fname = "%s_%s_structure_final_report.%%s"%(gathered_info['gene'],mutation)
  else:
    base_fname = "%s_structure_final_report.%%s"%gathered_info['gene']
 
  # The html directory, and it's large set of supporting files, is located beneath the directory
  # where this file, psb_rep.py, is located 
  src =  os.path.join(os.path.dirname(os.path.realpath(__file__)),"html")
  dest = "%s/html"%final_report_directory
  logger.info("Copying supporting javascript and html from %s to %s"%(src,dest))

  
  try:
    shutil.rmtree(dest,ignore_errors=True)
    shutil.copytree(src,dest)
  except Exception as err:
    print err
    print "REMINDER: Need to exit"
    # sys.exit(1)


  
  src_file_list = [os.path.join(root,name) for root,dirs,files in os.walk(src) for name in files]
  
  
  dest_file_list = [os.path.join(root,name) for root,dirs,files in os.walk(dest) for name in files]
  if len(src_file_list) == len(dest_file_list):
    infostr = "Copy of %d html support files succeeded"%len(src_file_list)
    print infostr
    logger.info(infostr)
    # print '\n'.join(src_file_list)
  else:
    logger.critical("FAILURE to cp html support files (%d source files,%d dest files)"%(len(src_file_list),len(dest_file_list)))
    sys.exit(1)
  
  # Now give read access to everyone, read-write to group, execute on directories
  capra_group = grp.getgrnam('capra_lab').gr_gid

  os.chmod(dest,0o775) 
  os.chown(dest, -1, capra_group)
  
  for root,dirs, files in os.walk(dest):
    for d in dirs:
      os.chmod(os.path.join(root,d),0o775) 
      os.chown(os.path.join(root, d), -1, capra_group)
    for f in files:
      fname = os.path.join(root, f)
      os.chmod(fname,0o664)
      os.chown(os.path.join(fname), -1, capra_group)
  
  env        = Environment(loader=FileSystemLoader(os.path.join(os.path.dirname(os.path.realpath(__file__)))))
  template   = env.get_template("psb_report.html")
  html_out = template.render(template_vars)

  templatePfamGraphics = env.get_template("html/pfamGraphicsIframeTemplate.html")
  htmlPfamGraphics = templatePfamGraphics.render(template_vars)

  save_cwd = os.getcwd();
  os.chdir(final_report_directory);
  
  # WE ARE NOW OPERATING FROM ..../UDN/CaseName target directory
  logger.info("Now working in %s",final_report_directory)

  html_fname = base_fname%"html"
  with open(html_fname,"w") as html_f:
    html_f.write(html_out)

  with open(pfamGraphicsIframe_fname,"w") as html_f:
    html_f.write(htmlPfamGraphics)

  pdf_fname = base_fname%"pdf"
  wkhtmltopdf_fname = base_fname%"wkhtml.pdf"
  print "\nWriting final reports to %s directory:\n  %s\n  %s\n  %s\n  %s\n"%(final_report_directory,pdf_fname,wkhtmltopdf_fname,html_fname,pfamGraphicsIframe_fname)

  logger.warning("Temporarily disabling all logging prior to calling weasyprint write_pdf()")

  logging.disable(sys.maxint)
  HTML(string=html_out).write_pdf(pdf_fname,stylesheets=["./html/css/typography.css"])#,CSS(string="@page { size: A3 landscape; }")])
  logging.disable(logging.NOTSET)
  logger.warning("weasyprint write_pdf() has returned.  Restoring logging")
     
  
  # Write out another .pdf using the wkhtmltopdf tool

  wkhtml_command_list = ["wkhtmltopdf","--print-media-type","--no-outline","--minimum-font-size","12",html_fname,wkhtmltopdf_fname]
  process = Popen(wkhtml_command_list, stdout = PIPE,stderr = PIPE)
  (output,err) = process.communicate()

  # So unfortunately, the stderr from wkhtmltopdf is full of ==================\r things, which are a mess for the logger.  Get rid of stuff that ends in \r  Keep all other lines and ending in \r
  # import pdb; pdb.set_trace()
  err_nl_array = err.split('\n')
  err_legit_list = []
  for ena in err_nl_array:
    crs = ena.split('\r')
    if len(crs) > 0: ## Pick off the last string, which means all ending with Carriage Return are ignored
      cr_last = crs[-1].strip()
      if cr_last:
        err_legit_list.append(cr_last)
      
  err_legit = '\n'.join(err_legit_list)
  
 
  # ContentNotFound really just info - not warning 
  if (process.returncode != 0) and ("Exit with code 1 due to network error: ContentNotFoundError" not in err):
    print "Unable to complete %s exitstatus=%d due to error %s\n  output: %s\n"%(str(wkhtml_command_list),process.returncode,err_legit,output)
    logger.warning("wkhtmltopdf stderr:\n%s",err_legit)
  else:
    logger.info("wkhtmltopdf stderr:\n%s",err_legit)
  
  # WE HAVE NOW RETURNED TO the psb_pipeline/bin directory
  os.chdir(save_cwd);
  # Close out the local log file for this mutation
  local_fh.flush()
  local_fh.close()
  logger.removeHandler(local_fh)
  gathered_info['final_report_directory'] = final_report_directory
  gathered_info['html_fname'] = html_fname
  return gathered_info # End of function report_one_mutation()
  
# Main logic here.  A period in the argument means the user wants to launch one mutation only,
# directly from a single mutation output file of psb_plan.py
if oneMutationOnly:
  if args.slurm:
    logger.warning("--slurm option is ignored when only one report is requested")
  template_vars = report_one_mutation(args.projectORstructures,args.workstatus)
  if template_vars: # The usual case, we had some Pathprox outputs
    print "Single mutation report saved to %s/.pdf/.wkhtml.pdf"%template_vars['html_fname']
  else:
    print "Due to lack of pathprox outputs, no html (or pdf) reports were created from %s"%args.workstatus
else:
  udn_csv_filename = os.path.join(collaboration_dir,"%s_missense.csv"%args.projectORstructures) # The argument is an entire project UDN124356
  msg = "Retrieving project mutations from %s"%udn_csv_filename
  if not infoLogging:
    print msg
  logger.info(msg)
  df_all_mutations = pd.DataFrame.from_csv(udn_csv_filename,sep=',')

  if args.slurm:
    slurm_directory = os.path.join(collaboration_dir,"slurm")
    slurm_file = os.path.join(slurm_directory,"psb_reps.slurm")
    slurm_stdout_directory = os.path.join(slurm_directory,"stdout")
    if not os.path.exists(slurm_directory):  # make the containing directory if needed
      os.makedirs(slurm_directory)
    if not os.path.exists(slurm_stdout_directory):  # make the stdout directory if needed
      os.makedirs(slurm_stdout_directory)
    jobCount = len(df_all_mutations) 
    msg = "Slurm script to run %d reports for %s is %s"%(jobCount,args.projectORstructures,slurm_file)
  else: # not a slurm run - so just do one report after another- normal case 
    msg = "Reporting on all %d project %s mutations"%(len(df_all_mutations),args.projectORstructures)
  if not infoLogging:
    print msg
  logger.info(msg)

  mutation_summaries = []
  slurm_array = []
  for index,row in df_all_mutations.iterrows():
    msg = "%s %-10s %-10s %-6s"%("Generating slurm entry for" if args.slurm else "Reporting on", row['gene'],row['refseq'],row['mutation'])
    if not infoLogging:
      print msg
    logger.info(msg)
    mutation_dir = os.path.join(collaboration_dir,"%s_%s_%s"%(row['gene'],row['refseq'],row['mutation']))
    if not os.path.exists(mutation_dir):  # python 3 has exist_ok parameter... 
      logger.critical("The specific mutation directory %s should have been created by psb_status.py.  Fatal problem.")
      sys.exit(1)

    gene_refseq_mutation = "%s_%s_%s"%(row['gene'],row['refseq'],row['mutation'])
       
    structure_report_filename = "%s/%s_structure_report.csv"%(mutation_dir,gene_refseq_mutation)
    workstatus_filename = "%s/%s_workstatus.csv"%(mutation_dir,gene_refseq_mutation)

    if args.slurm:
      slurm_array.append("./psb_rep.py -c %s -u %s %s %s"%(args.config,args.userconfig,structure_report_filename,workstatus_filename))
    else:
      gathered_info = report_one_mutation(structure_report_filename,workstatus_filename)
      if gathered_info:
        gathered_info['#'] = index
        if not gathered_info['Error'] and gathered_info['final_report_directory']:
          gathered_info['html_fname'] = os.path.join(gathered_info['final_report_directory'],gathered_info['html_fname'])
          msg = "%s %s %s report saved to %s/.pdf/.wkhtml.pdf"%(row['gene'],row['refseq'],row['mutation'],gathered_info['html_fname'])
          if not infoLogging:
            print msg
          logger.info(msg)
        gathered_info['Gene Mutation'] = "%-9.9s %-10.10s %-10.10s %-7.7s"%(row['gene'],row['refseq'],row['unp'],row['mutation'])
        gathered_info['Gene Hit Generic'] = ''
        mutation_summaries.append(gathered_info)
      else:
        msg = "Due to lack of pathprox or ddG outputs, %s %s %s has no html (or pdf) report"%(row['gene'],row['refseq'],row['mutation'])
        if not infoLogging:
          print msg
        logger.info(msg)

  if args.slurm:
    with open(slurm_file,"w") as slurmf:
      slurmf.write("""\
#!/bin/sh
#
# Project        : PSB Pipeline
# Filename       : %s
# Generated By   : %s
# For case       : %s
# Organization   : Vanderbilt Genetics Institute,
#                : Program in Personalized Structural Biology,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University
# Generated on   : %s
# Description    : Runs a python script on a cluster to launch psb_rep.py for 
#                : each and all mutations in parallel
#===============================================================================
# Slurm Parameters
"""%(slurm_file,__file__,args.projectORstructures,str(datetime.datetime.now())))

      slurmDict = SlurmParametersReport
      slurmDict['output'] = os.path.join(slurm_stdout_directory,"%s_%%A_%%a_psb_rep.out"%args.projectORstructures)
      slurmDict['job-name'] = "%s_psb_reps"%args.projectORstructures

      jobCount = len(slurm_array)
      if jobCount > 1:
        slurmDict['array'] = "0-%d"%(jobCount-1)
             
      # We can set a max on the number of simultaneously running jobs in the array.  For now, don't do this
      # if jobCount > 10:
      #   slurmDict['array'] += "%%10"
      
      for slurm_field in ['job-name','mail-user', 'mail-type', 'ntasks', 'time', 'mem', 'account','output','array']:
        if slurm_field in slurmDict:
          slurmf.write("#SBATCH --%s=%s\n"%(slurm_field,slurmDict[slurm_field]))

      if jobCount > 1:
        slurmf.write("""
echo "SLURM_ARRAY_TASKID="$SLURM_ARRAY_TASKID
""")
      slurmf.write("""
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
# echo "SLURMTMPDIR="$SLURMTMPDIR
echo "SLURM_SUBMIT_DIR = "$SLURM_SUBMIT_DIR

""")

      cwd = os.getcwd()

      slurmf.write("""
cd %s
if [ $? != 0 ]; then
echo Failure at script launch: Unable to change to directory %s
exit 1
fi
"""%(cwd,cwd))


      if jobCount == 1:
        # No need to fiddle with the slurm case statement
        slurmf.write(slurm_array[0])
      else:
        slurmf.write("\ncase $SLURM_ARRAY_TASK_ID in\n")
        slurm_array_id = 0
        for launchstring in slurm_array:
          # Write out the "case" tag 
          slurmf.write("%d)\n"%slurm_array_id)
          # Save this to put back later
          slurm_array_id += 1
          slurmf.write("echo Command: '%s'\n"%launchstring)
          slurmf.write(launchstring)
          # end the case tag
          slurmf.write("\n;;\n\n")
        slurmf.write("esac\n")
    msg = "Created slurm script to launch all psb_rep.py processes for this case: %s"%slurm_file   
    if args.verbose:
      logger.info(msg);
    else:
      print(msg)
       
  if mutation_summaries:
    # pprint.pformat(mutation_summaries)
    # Grab Souhrids gene interaction information
    geneInteractions_generic,html_table_generic,html_report_generic,text_report_generic = GeneInteractionReport(collaboration_dir,args.projectORstructures,False)
    geneInteractions_familial,html_table_familial,html_report_familial,text_report_familial = GeneInteractionReport(collaboration_dir,args.projectORstructures,True)
    if geneInteractions_generic:
      for mutation_summary in mutation_summaries:
         if mutation_summary['gene'] in geneInteractions_generic:
            mutation_summary['Gene Interactions Generic'] = ' '.join(geneInteractions_generic[mutation_summary['gene']])
    if geneInteractions_familial:
      for mutation_summary in mutation_summaries:
         if mutation_summary['gene'] in geneInteractions_familial:
            mutation_summary['Gene Interactions Familial'] = ' '.join(geneInteractions_familial[mutation_summary['gene']])
    env        = Environment(loader=FileSystemLoader(os.path.join(os.path.dirname(os.path.realpath(__file__)))))
    template   = env.get_template("case_report_template.html")
    # print html_table
    # import pdb; pdb.set_trace()
    final_gathered_info = {'mutation_summaries': mutation_summaries, 
             'firstGeneTable': html_table_generic,
             'firstGeneReport': html_report_generic,
             'secondGeneTable': html_table_familial,
             'secondGeneReport': html_report_familial,
             'case': args.projectORstructures, 
             "date"  : time.strftime("%Y-%m-%d")
             }
    html_out = template.render(final_gathered_info)
    case_summary_filename = os.path.join(collaboration_dir,"%s.html"%args.projectORstructures) # The argument is an entire project UDN124356
    with open(case_summary_filename,"w") as f:
      f.write(html_out)
    lastmsg = "The case summary report is: " + case_summary_filename
  else:
    lastmsg = "No mutation summaries - bummer"

  print ""
  print ""

  if args.verbose:
    logger.info(lastmsg);
  else:
    print(lastmsg)

  # The html directory, and it's large set of supporting files, is located beneath the directory
  # where this file, psb_rep.py, is located 
  src =  os.path.join(os.path.dirname(os.path.realpath(__file__)),"html")
  dest = "%s/html"%collaboration_dir
  logger.info("Copying supporting javascript and html from %s to %s"%(src,dest))
  
  try:
    shutil.rmtree(dest,ignore_errors=True)
    shutil.copytree(src,dest)
  except Exception as err:
    print err
    print "REMINDER: Need to exit"
 
sys.exit(0)
