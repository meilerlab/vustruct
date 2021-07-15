#!/usr/bin/env python3
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

print ("4 of 4  %s Pipeline report generator.  -h for detailed help"%__file__)

import logging,os,pwd,sys,grp,stat,string
import time, datetime
import argparse,configparser
import pprint
from logging.handlers import RotatingFileHandler
from logging import handlers

from psb_shared.ddg_repo import DDG_repo
from psb_shared.ddg_monomer import DDG_monomer

sh = logging.StreamHandler()
LOGGER = logging.getLogger()
LOGGER.addHandler(sh)

# Now that we've added streamHandler, basicConfig will not add another handler (important!)
log_format_string = '%(asctime)s %(levelname)-4s [%(filename)16s:%(lineno)d] %(message)s'
date_format_string = '%H:%M:%S'
log_formatter = logging.Formatter(log_format_string,date_format_string)

LOGGER.setLevel(logging.DEBUG)
sh.setLevel(logging.INFO)
sh.setFormatter(log_formatter)

rootdir_log_filename = "psb_rep.root.log"
needRoll = os.path.isfile(rootdir_log_filename)
rootdir_fh = RotatingFileHandler(rootdir_log_filename, backupCount=7)
formatter = logging.Formatter('%(asctime)s %(levelname)-4s [%(filename)20s:%(lineno)d] %(message)s',datefmt="%H:%M:%S")
rootdir_fh.setFormatter(formatter)
rootdir_fh.setLevel(logging.INFO)
LOGGER.addHandler(rootdir_fh)

if needRoll:
  rootdir_fh.doRollover()

sys.stderr.write("Root (case) directory log file is %s\n"%rootdir_log_filename)

# As report is generated, build up a list of just the files needed to make a portable website
# This list omits the intermediate files left behind in the calculation
website_filelist=['html/']


import shutil
import numpy as np
import pandas as pd
from jinja2 import Environment, FileSystemLoader
from io import StringIO, BytesIO
import pycurl
import json
from xml2json import Cxml2json

import subprocess
# from subprocess import Popen, PIPE
from weasyprint import HTML,CSS

from psb_shared import psb_config

cmdline_parser = psb_config.create_default_argument_parser(__doc__,os.path.dirname(os.path.dirname(__file__)))


cmdline_parser.add_argument("-s","--slurm",
                            help="Create a slurm file to launch all the reports", action = "store_true")
cmdline_parser.add_argument("projectORstructures",type=str,
                             help="The project ID UDN123456 to report on all mutations.  Else, a specific structures file from psb_plan.py  Example: ....../UDN/UDN123456/GeneName_NM_12345.1_S123A_structure_report.csv",
                             default= os.path.basename(os.getcwd()),nargs='?')
cmdline_parser.add_argument("workstatus",nargs='?',type=str,help="Blank for all mutations, or specific output file from psb_monitor.py  Example: ....../UDN/UDN123456/GeneName_NM_12345.1_S123A_workstatus.csv")
args,remaining_argv = cmdline_parser.parse_known_args()

# A period in the argument means the user wants to monitor one mutation only,
# directly from a single mutation output file of psb_launch.py
oneMutationOnly = ('.' in args.projectORstructures and os.path.isfile(args.projectORstructures) and args.workstatus != None)
infoLogging = False

if args.debug:
  infoLogging = True
  sh.setLevel(logging.DEBUG)
elif args.verbose:
  infoLogging = True
  sh.setLevel(logging.INFO)

# Load the structure details of all jobs that were not dropped

# print "Command Line Arguments"
LOGGER.info("Command Line Arguments:\n%s"%pprint.pformat(vars(args)))

if not os.path.exists(args.config):
  LOGGER.critical("Global config file not found: " + args.config)
  sys,exit(1)

required_config_items = ['output_rootdir','collaboration','ddg_config']

config,config_dict = psb_config.read_config_files(args,required_config_items)
config_dict_shroud_password = {x:config_dict[x] for x in required_config_items}
dbpass = config_dict.get('dbpass','?')
config_dict_shroud_password['dbpass'] = '*' * len(dbpass)
LOGGER.info("Configuration File parameters:\n%s"%pprint.pformat(config_dict_shroud_password))

try:
  slurmParametersAll = dict(config.items("SlurmParametersAll"))
except Exception as ex:
  slurmParametersAll = {}

import copy
slurmParametersReport = copy.deepcopy(slurmParametersAll)

try:
  slurmParametersReport.update(dict(config.items("SlurmParametersReport")))
except Exception as ex:
  pass

config_pathprox_dict = dict(config.items("PathProx"))  
LOGGER.info("Pathprox config entries:\n%s"%pprint.pformat(config_pathprox_dict))

# The collaboration_dir is the master directory for the case, i.e. for one patient
# Example: /dors/capra_lab/projects/psb_collab/UDN/UDN532183
udn_root_directory = os.path.join(config_dict['output_rootdir'],config_dict['collaboration'])

if oneMutationOnly:
  collaboration_absolute_dir = os.path.dirname(os.path.dirname(args.projectORstructures))
else:
  collaboration_absolute_dir = os.path.join(udn_root_directory,args.projectORstructures)

collaboration_absolute_dir = os.path.realpath(collaboration_absolute_dir)
LOGGER.info("Collaboration_absolute_dir = %s",collaboration_absolute_dir)

collaboration_dir = os.path.relpath(os.getcwd(),collaboration_absolute_dir)
LOGGER.info("Relative to cwd=%s, collaboration_directory= %s",os.getcwd(),collaboration_dir)

def copy_html_css_javascript():
  # The source html directory, and it's large set of supporting files, is located beneath the directory
  # where this file, psb_rep.py, is located.  All mutations share this one large html/css/javascript
  # sourcee directory
  src =  os.path.join(os.path.dirname(os.path.realpath(__file__)),"html")
  dest = "%s/html"%collaboration_dir
  LOGGER.info("Copying supporting javascript and html from %s to %s"%(src,dest))
  
  try:
    shutil.rmtree(dest,ignore_errors=True)
    shutil.copytree(src,dest)
  except Exception as e:
    LOGGER.exception(e)
    sys.exit(1)

  # Make sure 
  src_file_list = [os.path.join(root,name) for root,dirs,files in os.walk(src) for name in files]
  
  dest_file_list = [os.path.join(root,name) for root,dirs,files in os.walk(dest) for name in files]
  if len(src_file_list) == len(dest_file_list):
    infostr = "Copy of %d html support files succeeded"%len(src_file_list)
    print (infostr)
    LOGGER.info(infostr)
    # print '\n'.join(src_file_list)
  else:
    LOGGER.critical("FAILURE to cp html support files (%d source files,%d dest files)"%(len(src_file_list),len(dest_file_list)))
    sys.exit(1)
  
  # Now give read access to everyone, read-write to group, execute on directories
  try:
      capra_group = grp.getgrnam('capra_lab').gr_gid
  except KeyError:
      capra_group = os.getegid()

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

# Return a dictionary of "generic Interactions" that can flow into the summary report
# Return an html table that can be placed at end of summary report
# Return html long report
# Return text that can be dumped to the log file
def GeneInteractionReport(case_root,case,CheckInheritance):
    genepair_dict_file = os.path.join(case_root,'casewide','DigenicAnalysis','%s_all_gene_pairs_summary.json'%case)

    try:
        with open(genepair_dict_file,'r') as fd:
          genepair_dict = json.load(fd)
    except:
        LOGGER.exception("Cannot open %s.  Did you run Souhrid's analysis scripts?"%genepair_dict_file)
        return None,None,None,None
       
    # genes_filename = os.path.join(case_root,"%s_genes.txt"%case)
    # try:
    #    with open(genes_filename) as fd:
    #       unstripped_genes = fd.readlines()
    # except:
    #    LOGGER.exception("Cannot open %s.  Did you run Souhrid's analysis scripts?"%genes_filename)
    #     return None,None,None,None

    geneInteractions = {}
    html_table =  StringIO()
    html_report = StringIO()
    text_report = StringIO()

    structure_positive_genes = [gene.strip() for gene in list(sorted(genepair_dict.keys()))]
    
    pairs = {}
    pairs['direct'] = []
    pairs['pathways'] = []
    pairs['phenotypes'] = []
    pairs['other'] = []

    # We let the template html open the <table> with formatting as it prefers
    # print >>html_table, "<table>" 
    print("<thead><tr>", file=html_table)
    print("<th></th>", file=html_table) # Top left is blank cell in header
    for gene1 in structure_positive_genes:
      print("<th>" + gene1 + "</th>" , file=html_table)
    print("</tr></thead>", file=html_table)

    for gene1 in structure_positive_genes:
        print("%s:"%gene1, file=text_report)
        print('%s:'%gene1, file=html_report)
        print("<tr><td>%s</td>"%gene1, file=html_table)
        if gene1 in genepair_dict:
          print('<ul style="list-style: none;">', file=html_report)
        for gene2 in structure_positive_genes:
          print("<td>", file=html_table)
          if gene1 in genepair_dict and gene2 in genepair_dict[gene1] and ((not CheckInheritance) or ('nherit' in genepair_dict[gene1][gene2])):
              pairs_hits = []
              if 'direct' in genepair_dict[gene1][gene2]:
                  pairs_hits.append('direct')
              if 'pwy' in genepair_dict[gene1][gene2]:
                  pairs_hits.append('pathways')
              if 'pheno' in genepair_dict[gene1][gene2]:#  and genepair_dict[gene1][gene2]['Patient_Phenotype_Overlap'] > 0:
                  pairs_hits.append('phenotypes')
              # if genepair_dict[gene1][gene2]['Distance'] < 3 and genepair_dict[gene1][gene2]['Coexpression'] < 50:
              #    pairs_hits.append('other')
              if pairs_hits:
                  for pairs_hit in pairs_hits:
                      pairs[pairs_hit].append((gene1,gene2))
                      if pairs_hit != 'phenotypes': # We don't get too excited for phenotypes at this point
                          if not gene1 in geneInteractions:
                              geneInteractions[gene1] = []
                          if not gene2 in geneInteractions[gene1]:
                              geneInteractions[gene1].append(gene2)
                  logging.debug("Gene Interactions: %s %s %s"%(gene1, gene2, str(pairs_hits)))
                  print(gene2, str(pairs_hits) , file=text_report)
                  print("<li>",gene2 + ':', ", ".join(pairs_hits),"</li>", file=html_report)
                  print('<br>'.join(pairs_hits), file=html_table)
          print("</td>", file=html_table)
        print('</ul>', file=html_report)
        print("</tr>", file=html_table)
    # print("</table>", file=html_table)
    return geneInteractions,html_table.getvalue(),html_report.getvalue(),text_report.getvalue()

# This complex routine gather PathProx, and ddG results from the pipeline run, where available
# It creates a sigle .html file representing the mutation, and returns a dictionary that is
# input to the overall case report in gathered_info
def report_one_mutation(structure_report,workstatus_filename):
    gathered_info = {} # Dictionary info to return, of great use by caller
    # Load the status of jobs that was created by psb_launch.py (or prior run of psb_monitor.py)
    df_all_jobs_status = pd.read_csv(workstatus_filename,'\t',keep_default_na = False,na_filter=False,dtype=str)
    msg = "%d rows read from work status file %s"%(len(df_all_jobs_status),workstatus_filename)
    if len(df_all_jobs_status) < 1:
      LOGGER.critical(msg)
      return None
    else:
      if not infoLogging:
        print (msg)
      LOGGER.info(msg)

    web_dir=os.path.dirname(structure_report)

    row0 = df_all_jobs_status.iloc[0]
    gathered_info = {col: row0[col] for col in ['project','unp','gene','refseq','mutation']}
    if 'unp' in row0:
        gathered_info['unp'] = row0['unp']
    else: # Get this from PDBMapProtein
        LOGGER.critical('FIX THIS DURING PLAN PHASE BY PUTTING UNP COLUMN OUT!!')
        gathered_info['unp'] = PDBMapProtein.refseqNT2unp(gathered_info['refseq'])[0]
    
    template_vars = gathered_info.copy()

    mutation = row0['mutation'] # Sorry for alias - but this is all over the code

    gathered_info['Error'] = '' # < This is False for truth but looks fine on the report
    gathered_info['ddG Max'] = gathered_info['ddG Min'] = None
    gathered_info['disease2_pp'] = gathered_info['disease1_pp'] = None


    udn_sequence_annotation_jobs = ((df_all_jobs_status['flavor'] == 'SequenceAnnotation'))
    if udn_sequence_annotation_jobs.sum() > 0:
      LOGGER.info("Dropping %d rows where jobs are UDN SequenceAnnotation"%udn_sequence_annotation_jobs.sum())
      df_all_jobs_status = df_all_jobs_status[~udn_sequence_annotation_jobs]


    # ExitCode column arrives when monitor program finds a clear exit, but may not arrive
    # if programs exit without leaving a trail of info 
    if not 'ExitCode' in df_all_jobs_status.columns: 
      gathered_info['Error'] = "'ExitCode' column missing."
      LOGGER.error(gathered_info['Error'] + "in " + workstatus_filename)
      return gathered_info

    failed_jobs = df_all_jobs_status['ExitCode'] > '0'
    if failed_jobs.sum() > 0:
      LOGGER.info("Dropping %d rows for jobs that failed (ExitCode > '0')."%failed_jobs.sum())
      df_all_jobs_status = df_all_jobs_status[~failed_jobs]

    incomplete_jobs = df_all_jobs_status['ExitCode'] != '0'
    if incomplete_jobs.sum() > 0:
      LOGGER.info("Dropping %d rows for incomplete jobs (ExitCode missing).\n%s",incomplete_jobs.sum(),df_all_jobs_status[incomplete_jobs])
      df_all_jobs_status = df_all_jobs_status[~incomplete_jobs]
      if len(df_all_jobs_status) < 1:
        gathered_info['Error'] = "All jobs are marked incomplete."
        LOGGER.error(gathered_info['Error'] + "  Be sure to run pdb_monitor.py before this script")
        return gathered_info

    if len(df_all_jobs_status) < 1:
        gathered_info['Error'] = "No completed PathProx or ddG jobs."
        LOGGER.warning(gathered_info['Error'] + " in " + workstatus_filename)
        return gathered_info

    # Sanity check - make sure mutation is same in all rows!
    test_count = len(df_all_jobs_status.groupby('mutation'))
    if test_count != 1:
      gathered_info['Error'] = "%s  may be corrupted, has some rows which lack the same mutation"%workstatus_filename
      LOGGER.error(gathered_info['Error'])
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
      LOGGER.critical(errorMsg)
      sys.stderr.write(errorMsg)
      sys.exit(1)

    # pw_name = pwd.getpwuid( os.getuid() ).pw_name # example jsheehaj or mothcw
    log_filename = os.path.join(mutation_log_dir,"psb_rep.log")

    if oneMutationOnly:
      sys.stderr.write("psb_rep log file is %s\n"%log_filename)
    else:
      LOGGER.info("additionally logging to file %s"%log_filename)
    needRoll = os.path.isfile(log_filename)
    local_fh = RotatingFileHandler(log_filename, backupCount=7)
    formatter = logging.Formatter('%(asctime)s %(levelname)-4s [%(filename)20s:%(lineno)d] %(message)s',datefmt="%H:%M:%S")
    local_fh.setFormatter(formatter)
    local_fh.setLevel(logging.INFO)
    LOGGER.addHandler(local_fh)

    if needRoll:
      local_fh.doRollover()
# Load the structure details of all jobs that were not dropped

# Load the structure details of all jobs that were not dropped
    nan_values = ['-1.#IND', '1.#QNAN', '1.#IND', '-1.#QNAN', '#N/A','N/A', '#NA', 'NULL', 'NaN', '-NaN', 'nan', '-nan']
    df_structure_report = pd.read_csv(structure_report,'\t',keep_default_na=False,na_values = nan_values,
         # EM structures often have integer chain_ids that confuse the pipeline
         # Resolution and Seq Identites are floats, but nanS for some kinds of models and structures
         dtype={'Resolution': float, 'Template Identity': float,'Trans Identity': float, 'chain_id': str})

    df_structure_report['structure_id'].replace("^(.*)\.pdb$",r"\1", regex=True,inplace=True)

    # Create a column for PDB Pos in format pdb_res+icode+a chain id iff there are multi-chains
    df_structure_report['PDB Pos'] = df_structure_report[['mut_pdb_res','mut_pdb_icode','chain_id','biounit_chains']].apply(lambda x: str(x[0])+str(x[1]).strip()+".%s"%x[2] if x[3] > 1 else '', axis=1)

    LOGGER.info("%d rows read from structure_report file %s"%(len(df_structure_report),structure_report))
    df_structure_report.set_index(["method","structure_id","chain_id","mers"],inplace = True)
  
    # for the .html report for each variant, use the direction in which
    # The workstation_filename was retrieved
    # import pdb; pdb.set_trace()
    variant_report_directory = os.path.realpath(os.path.dirname(workstatus_filename))
    
    # project = 'TestJED' # MUST FIX THIS ERROR
    # unp = 'Q9NY15-1' # MUST FIX THIS ERROR
    df_all_jobs_status.set_index('uniquekey',inplace=True)
    
    df_complete_jobs = df_all_jobs_status[df_all_jobs_status['ExitCode'] == '0'].copy()
    
    
    print ("Generating html and pdf final reports for %d of %d complete jobs:"%(len(df_complete_jobs),len(df_all_jobs_status)))
  
    # Rename the columns on a duplicate copy, to make things
    # More Human readable, and honor the history of this report:
    xlate_columns_dict={'label':"Label",
        'analyzable':"Analyzable?",
        'distance_from_aligned':"Distance to Boundary",
        'perc_identity':'Trans Identity',
        'pdb_template':'PDB Template',
        'template_identity':'Template Identity',
        'biounit_chains': 'Chains',
        'nresidues':'Residues',
        'resolution':'Resolution (PDB)',
        'trans_first_aligned': 'Seq Start',
        'trans_last_aligned': 'Seq End',
        'structure_url': 'URL'}

    def structure_html_vars(method_pdbid_chain_mers_tuple,df_unique_structure:pd.DataFrame,df_structure_report_onerow:pd.Series,thestruct):
      """Subset the dataframe to the results of this chain
      This could need 2 rows - one for each pathprox - we just have to keep an eye on this"""
      map_punct_to_None = dict.fromkeys(string.punctuation) # for each punctuation character as key, the value is None
      structure  = {"structid":thestruct['pdbid'],
        "chain":thestruct['chain'],
        "mers":thestruct['mers'],
        "html_div_id":(thestruct['pdbid'] + thestruct['chain'] + (thestruct['mers'] if thestruct['mers'] != 'monomer' else '')).translate(str.maketrans(map_punct_to_None)),
        "mutation":gathered_info['mutation']}

    
      structure.update(df_structure_report_onerow.to_dict())
    
      # one_job_status = df_onerow.to_frame()
      # one_structure_report = structure_report_onerow.toframe()
    
      # Make a one row dataframe just to convert to html
      # tempdf = df_structure_report_onerow
      
      human_df = df_structure_report_onerow.copy(deep=True).rename(xlate_columns_dict)
      dfT = human_df.to_frame().transpose()[ list(xlate_columns_dict.values()) ]
      # Not sure this is quite right
      # if df_structure_report_onerow['label'] == 'pdb':
      #     dfT['Seq Id'] = 100.0

      # Blank out transcript identity for models
      if df_structure_report_onerow['label'] not in ['pdb','biounit']:
          dfT['Trans Identity'] = None

     
      structure["details"] = dfT.to_html(float_format=lambda x: "%.2f"%x)
      return structure
    
    def ddG_html_vars(structure,thestruct):
        if not 'ddG_monomer' in thestruct:
            return structure,None

        df_row = thestruct['ddG_monomer']

        ddg_repo = DDG_repo(config_dict['ddg_config'],
                    calculation_flavor='ddg_monomer')
      
        output_flavor_directory = os.path.join(df_row['outdir'],'ddG','ddg')

        structure_id = df_row['Structure ID']
        chain_id = df_row['Chain']
        if df_row['method'] == 'swiss':
            ddg_structure_dir = ddg_repo.set_swiss(structure_id,chain_id)
        elif df_row['method'] == 'modbase':
            ddg_structure_dir = ddg_repo.set_modbase(structure_id,chain_id)
        elif df_row['method'] == 'usermodel':
            ddg_structure_dir = ddg_repo.set_usermodel(structure_id,chain_id)
        else: # Has to be a pdb
            ddg_structure_dir = ddg_repo.set_pdb(structure_id.lower(),chain_id)

        variant = df_row['pdbmut']
        ddg_repo.set_variant(variant)
        ddg_monomer = DDG_monomer(ddg_repo,variant)

        ddg_results_df = ddg_monomer.retrieve_result()
        if ddg_results_df is None:
            LOGGER.info('No results found for %s.%s:%s'%(structure_id,chain_id,variant))
            return structure, None

        return structure,ddg_results_df

  
    # Compile PathProx results into report (one page per structure+chain+mutation)
    def PathProx_html_vars(structure,thestruct,disease1_or_disease2):
      """ Generate the HTML for a single page structure report """
    
      variant_sql_label = config_pathprox_dict['%s_variant_sql_label'%disease1_or_disease2]
      variant_short_description = config_pathprox_dict['%s_variant_short_description'%disease1_or_disease2]
      # If there is no thestruct['Clinvar'] or thestruct['COSMIC'] then return right away - no big deal

      if not disease1_or_disease2 in thestruct:
        return structure,None
    
      df_row = thestruct[disease1_or_disease2] # was Clinvar_or_COSMIC]
    
      save_cwd = os.getcwd()
      os.chdir(os.path.dirname(df_row['outdir']))


      # Define directory prefixes
      # This is where we'll search for those final pathprox .pdb files, etc 
      output_flavor_directory = os.path.join(os.path.basename(df_row['outdir']),df_row['flavor'])

      PathProxPrefix = "%s_%s_%s_%s_D"%(thestruct['pdbid'],thestruct['chain'],variant_sql_label,config_pathprox_dict['neutral_variant_sql_label'])
    
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
        summary["Gene"] = gathered_info['gene']
        summary["Protein"] = gathered_info['unp']
        # Dummy code coverage-dependent columns if missing
        if "%s_pathprox"%mutation not in summary.columns:
          summary["%s_pathprox"%mutation] = np.nan
        if "%s_neutcon"%mutation not in summary.columns:
          summary["%s_neutcon"%mutation] = np.nan
        if "%s_pathcon"%mutation not in summary.columns:
          summary["%s_pathcon"%mutation] = np.nan
        summary = summary.rename(columns=coldir)
        summary = summary[colord]
        structure["%s_results"%disease1_or_disease2] = summary.to_html(index=False,float_format=lambda x: "%.2f"%x)
      else:
        print ("PathProx summary file not found: %s"%pathprox_summary_file)
        structure["%s_results"%disease1_or_disease2] = ""
        summary = pd.DataFrame()
      """
      # Load the COSMIC summary of results for this mutation
      if os.path.exists("%s/%s_summary.txt"%(dors_ppdir,disease2_prefix)):
        disease2_sum = pd.read_csv("%s/%s_summary.txt"%(dors_ppdir,disease2_prefix),sep='\t',header=0)
        disease2_sum["Mutation"] = mutation
        # Dummy code coverage-dependent columns if missing
        if "%s_pathprox"%mutation not in disease2_sum.columns:
          disease2_sum["%s_pathprox"%mutation] = np.nan
        if "%s_neutcon"%mutation not in disease2_sum.columns:
          disease2_sum["%s_neutcon"%mutation] = np.nan
        if "%s_pathcon"%mutation not in disease2_sum.columns:
          disease2_sum["%s_pathcon"%mutation] = np.nan
        disease2_sum = disease2_sum.rename(columns=coldir)
        disease2_sum = disease2_sum[colord]
        structure["disease2_results"] = disease2_sum.to_html(index=False,float_format=lambda x: "%.2f"%x)
      else:
        print "COSMIC PathProx summary file not found: %s/%s_summary.txt"%(dors_ppdir,disease2_prefix)
        structure["disease2_results"] = ""
        disease2_sum = pd.DataFrame()
      """
      # Load the variant plots
      # Use  root-relative pathnames as inputs to the structure[] diction that flows into the html creation
      # vus_var and exac_var should be populated from either Clinvar or COSMIC calculations.  Get them we don't
      # have them yet

      # If an image file is in the filesystem, great, else return empty string  
      def check_existence_add_to_website(plot_png,structure_dict,structure_key,value_if_missing=None):
        if os.path.exists(plot_png):
          structure_dict[structure_key] = plot_png
          website_filelist.append(os.path.join(web_dir,plot_png))
          return plot_png
        LOGGER.info("Graphic png file %s not found",plot_png)
        if value_if_missing is not None:
          structure_dict[structure_key] = value_if_missing

        return 0 

      if not "vus_var" in structure:
        check_existence_add_to_website("%s/%s_structure.png"%(output_flavor_directory,PathProxPrefix),structure,"vus_var",None)
    
      if not "neutral_var" in structure:
        check_existence_add_to_website("%s/%s_neutral.png"%(output_flavor_directory,PathProxPrefix),structure,"neutral_var",None)
 
      check_existence_add_to_website("%s/%s_pathogenic.png"%(output_flavor_directory,PathProxPrefix),structure,"%s_var"%disease1_or_disease2,0)

      # Where pathprox has left us ngl viewer variant setups, incorporate those into the final report
      # As with the static .png im

      def residue_list_to_ngl(residues_onechain,chain_id):
        if len(residues_onechain) == 0:
            return 0,None
        ngl_res_string = ''
        for transcript_resno in residues_onechain:
            if ngl_res_string:
                ngl_res_string += ' '
            ngl_res_string += "%d"%transcript_resno
            if chain_id:
                ngl_res_string += ":%s"%chain_id
        return len(residues_onechain),ngl_res_string

      # def residue_list_to_ngl_CAs(residues_onechain):
      #    ngl_res_count,ngl_res_string = residue_list_to_ngl(residues_onechain)
      #    return ngl_res_count,"(" + ngl_res_string + ") and .CA" if ngl_res_string else ngl_res_string

      def residues_to_ngl(residues):
         ngl_res_count = 0
         ngl_res_string = None
         for chain_id in residues:
             ngl_res_chain_count,ngl_res_chain_string =  residue_list_to_ngl(residues[chain_id],chain_id)
             if ngl_res_chain_count and ngl_res_chain_string:
                 if ngl_res_string:
                     ngl_res_string += ' '
                 else:
                     ngl_res_string = ''
  
                 ngl_res_count += ngl_res_chain_count
                 ngl_res_string += ngl_res_chain_string
         return ngl_res_count,ngl_res_string

      def residues_to_ngl_CAs(residues):
          ngl_res_count,ngl_res_string = residues_to_ngl(residues)
          return ngl_res_count,"(" + ngl_res_string + ") and .CA" if ngl_res_string else ngl_res_string
       
      json_filename = os.path.join(output_flavor_directory,PathProxPrefix + "_ResiduesOfInterest.json")
      if not os.path.exists(json_filename):
          LOGGER.warning("Can't open %s for processing",json_filename)
      if os.path.exists(json_filename):
        with open(json_filename) as f:
          residuesOfInterest = json.load(f)
          # if 'cifSSfilename' in residuesOfInterest:
          #    cifSSbasename = os.path.basename(residuesOfInterest['cifSSfilename'])
          #    cifSSdirname = os.path.realpath(os.path.dirname(residuesOfInterest['cifSSfilename']))
          #    residuesOfInterest['cifSSfilename'] = os.path.join(os.path.relpath(cifSSdirname,variant_report_directory),cifSSbasename)
          if 'pdbSSfilename' in residuesOfInterest:
              pdbSSbasename = os.path.basename(residuesOfInterest['pdbSSfilename'])
              pdbSSdirname = os.path.realpath(os.path.dirname(residuesOfInterest['pdbSSfilename']))
              residuesOfInterest['pdbSSfilename'] = os.path.join(os.path.relpath(pdbSSdirname,variant_report_directory),pdbSSbasename)
          LOGGER.info('Loaded %d variants, %d neutrals, %d pathogenic from %s'%(
            len(residuesOfInterest['variants']),
            len(residuesOfInterest['neutrals']),
            len(residuesOfInterest['pathogenics']),
            json_filename))

        # For final .html - we need to be specific regarding type of pathogenic.  We weren't when json file was written
        residuesOfInterest['%s_pathogenics'%disease1_or_disease2] = residuesOfInterest['pathogenics']
        del residuesOfInterest['pathogenics']

        structure.update(residuesOfInterest)

        # import pdb; pdb.set_trace()
        website_filelist.append(os.path.join(web_dir,structure['pdbSSfilename']))
        # website_filelist.append(os.path.join(web_dir,structure['cifSSfilename']))

        structure['ngl_variant_residue_count'],structure['ngl_variant_residues'] = residues_to_ngl(residuesOfInterest['variants'])
        structure['ngl_neutral_residue_count'],structure['ngl_neutral_residues'] = residues_to_ngl_CAs(residuesOfInterest['neutrals'])
        structure['ngl_%s_residue_count'%disease1_or_disease2], \
          structure['ngl_%s_residues'%disease1_or_disease2] = residues_to_ngl_CAs(residuesOfInterest['%s_pathogenics'%disease1_or_disease2])

      
      # if not "vus_ngl_html" in structure and os.path.exists("%s/%s_structure.ngl.html"%(output_flavor_directory,PathProxPrefix)):
        # with open("%s/%s_structure.ngl.html"%(output_flavor_directory,PathProxPrefix),"rb") as ngl_html_file:
          # structure["vus_ngl_html"] =  ngl_html_file.read()
 
      if not "neutral_ngl_html" in structure:
        html_filename =  "%s/%s_neutral.ngl.html"%(output_flavor_directory,PathProxPrefix)
        if os.path.exists(html_filename):
          website_filelist.append(os.path.join(web_dir,html_filename))
          with open("%s/%s_neutral.ngl.html"%(output_flavor_directory,PathProxPrefix),"rb") as ngl_html_file:
            structure["neutral_ngl_html"] = ngl_html_file.read()

      pathogenic_ngl_html_filename = "%s/%s_pathogenic.ngl.html"%(output_flavor_directory,PathProxPrefix)
      if os.path.exists(pathogenic_ngl_html_filename):
        website_filelist.append(os.path.join(web_dir,pathogenic_ngl_html_filename))
        with open( pathogenic_ngl_html_filename ,"rb") as ngl_html_file:
          structure["%s_ngl_html"%disease1_or_disease2]    = ngl_html_file.read()

      # Load the Ripley's K results for Disease1(default ClinVar) and Disease2(default COSMIC)

      # Don't re-do neutral_K (exac) graphic if already loaded from the other directory
      if (not "neutral_K" in structure) or len(structure["neutral_K"]) == 0:
        check_existence_add_to_website("%s/%s_neutral_K_plot.png"%(output_flavor_directory,PathProxPrefix),structure,"neutral_K","")

      check_existence_add_to_website("%s/%s_pathogenic_K_plot.png"%(output_flavor_directory,PathProxPrefix),structure,"%s_K"%disease1_or_disease2,"")
      # Load the Ripley's D results for Disease1(ClinVar) and Disease2(COSMIC)
      check_existence_add_to_website("%s/%s_D_plot.png"%(output_flavor_directory,PathProxPrefix),structure,"%s_D"%disease1_or_disease2,"")
      # Load the PathProx Mapping and Performance for Disease1(ClinVar)
      check_existence_add_to_website("%s/%s_pathprox.png"%(output_flavor_directory,PathProxPrefix),structure,"%s_pp"%disease1_or_disease2,"")
      check_existence_add_to_website("%s/%s_pathprox_roc.png"%(output_flavor_directory,PathProxPrefix),structure,"%s_roc"%disease1_or_disease2,"")
      check_existence_add_to_website("%s/%s_pathprox_pr.png"%(output_flavor_directory,PathProxPrefix),structure,"%s_pr"%disease1_or_disease2,"")

      os.chdir(save_cwd)
      
      # Return the structure and the Disease1(ClinVar)/Disease2(COSMIC) results
      return structure,summary
    
    def unp2PfamDomainGraphicString(unp,timeoutSeconds):
      url = 'http://pfam.xfam.org/protein/{0}/graphic'.format(unp.split(' ')[0].split('-')[0])
    
      buffer = BytesIO()
    
      try:
        timeout = 60 #NOTE: using pycurl assures that this timeout is honored
        c = pycurl.Curl()
        c.setopt(pycurl.CONNECTTIMEOUT, timeoutSeconds)
        c.setopt(pycurl.TIMEOUT, timeoutSeconds)
        # c.setopt(pycurl.WRITEFUNCTION, lambda x: None)
        c.setopt(pycurl.WRITEFUNCTION, buffer.write)
        c.setopt(pycurl.HEADERFUNCTION, lambda x: None)
        c.setopt(pycurl.NOSIGNAL, 1)
        c.setopt(pycurl.URL, url)
        c.setopt(pycurl.HTTPGET, 1)

        c.perform()
        c.close()
    
      except Exception as ex:
       LOGGER.exception('Failed to get JSON string from %s: %s'%(url,str(ex)))
       return None
    
      JSONstr= buffer.getvalue().decode('latin')
      if len(JSONstr) > 3 and JSONstr[0] == '[':
        JSONstr = JSONstr[1:-1]
      return JSONstr
  
  
    df_complete_jobs['Structure ID'] = df_complete_jobs['pdbid']
    df_complete_jobs['Method'] = df_complete_jobs['method']
    df_complete_jobs['Chain'] = df_complete_jobs['chain']
    df_complete_jobs['Complex'] = df_complete_jobs['mers']
    
    tdf        = df_complete_jobs.set_index(["Method","Structure ID","Chain","Complex"])
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
      method_pdbid_chain_mers_tuple = (row['Method'],row['pdbid'],row['chain'],row['mers'])
      thestruct = struct_dict.get(method_pdbid_chain_mers_tuple,{})
      thestruct['mutation'] = row['mutation']
      thestruct['pdbid'] = row['pdbid']
      thestruct['chain'] = row['chain']
      thestruct['mers'] = row['mers']
      if row['flavor'].endswith(config_pathprox_dict['disease1_variant_sql_label']):
        thestruct['disease1'] = row.to_dict()
      elif row['flavor'].endswith(config_pathprox_dict['disease2_variant_sql_label']):
        thestruct['disease2'] = row.to_dict()
      elif 'SequenceAnnotation' == row['flavor']: 
        continue # UDN Sequence annotations are NOT a part of generated reports
      elif 'ddG_monomer' == row['flavor']: 
        thestruct['ddG_monomer'] = row.to_dict()
      else:
        LOGGER.critical("flavor in row is neither %s nor %s nor ddG_monomer - cannot continue:\n%s",config_pathprox_dict['disease1_variant_sql_label'],config_pathprox_dict['disease2_variant_sql_label'],str(row))
        sys.exit(1)
      # This will often reassign over prior assignments - That's OK
      struct_dict[method_pdbid_chain_mers_tuple] = thestruct
   
    print("%d structures have at least one complete report"%len(struct_dict))
    assert len(struct_dict) == len(tdf_unique_structures)
    # print("*** THINK ABOUT THIS ASSERT BELOW ***")
    # assert len(struct_dict) == len(df_structure_report)
    
    # Generate vars for individual structure pages
    print("\nLoading the results from each structural analysis...")
    
    structures = []
    disease1_list = []
    disease2_list = []
    ddG_list = []

    # Load up variables for each structure in our set
    for method_pdbid_chain_mers_tuple in struct_dict:
      s = structure_html_vars(
          method_pdbid_chain_mers_tuple,
          tdf_unique_structures.loc[method_pdbid_chain_mers_tuple],
          df_structure_report.loc[method_pdbid_chain_mers_tuple],
          struct_dict[method_pdbid_chain_mers_tuple])

      s,temp_disease_sum = PathProx_html_vars(s,struct_dict[method_pdbid_chain_mers_tuple],'disease1')
      if temp_disease_sum is not None:
        disease1_list.append(temp_disease_sum)

      s,temp_disease_sum = PathProx_html_vars(s,struct_dict[method_pdbid_chain_mers_tuple],'disease2')
      if temp_disease_sum is not None:
        disease2_list.append(temp_disease_sum)

      del temp_disease_sum

      # import pdb; pdb.set_trace()
      s,ddG_summary = ddG_html_vars(s,struct_dict[method_pdbid_chain_mers_tuple])
      if ddG_summary is not None:
        ddG_list.append(ddG_summary)

      structures.append(s)

    if disease1_list:
      disease1_sum = pd.concat(disease1_list)
    else:
      disease1_sum = pd.DataFrame() 

    if disease2_list:
      disease2_sum = pd.concat(disease2_list)
    else:
      disease2_sum = pd.DataFrame() 

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
    if not disease1_sum.empty:
      msg = "Summarizing results from %d %s PathProx analyses..."%(len(disease1_sum),config_pathprox_dict['disease1_variant_short_description'])
      if not infoLogging:
        print (msg)
      LOGGER.info(msg)
      # Fancy way to compete the average for all the aggcols 
      disease1_sum = disease1_sum.groupby("Mutation")[aggcols].aggregate(np.nanmedian)
      disease1_auc = disease1_sum["ROC AUC"].values[0]
      gathered_info['disease1_pp']  = disease1_sum["PathProx"].values[0]
      disease1_auc = disease1_auc if not np.isnan(disease1_auc) else None
      gathered_info['disease1_pp']  = gathered_info['disease1_pp'] if not np.isnan(gathered_info['disease1_pp']) else None
    else:
      disease1_auc = gathered_info['disease1_pp'] = None
    
    if not disease2_sum.empty:
      msg = "Summarizing results from %d %s PathProx analyses..."%(len(disease2_sum),config_pathprox_dict['disease2_variant_short_description'])

      if not infoLogging:
        print (msg)
      LOGGER.info(msg)
      # Fancy way to compete the average for all the aggcols 
      disease2_sum = disease2_sum.groupby("Mutation")[aggcols].aggregate(np.nanmedian)
      disease2_auc = disease2_sum["ROC AUC"].values[0]
      gathered_info['disease2_pp']  = disease2_sum["PathProx"].values[0]
      disease2_auc = disease2_auc if not np.isnan(disease2_auc) else None
      gathered_info['disease2_pp']  = gathered_info['disease2_pp'] if not np.isnan(gathered_info['disease2_pp']) else None
    else:
      disease2_auc = gathered_info['disease2_pp'] = None


    if not ddG_sum.empty:
      msg = "Summarizing results from %d ddg runs..."%len(ddG_sum)
      if not infoLogging:
        print (msg)
      LOGGER.info(msg)

    LOGGER.warning("Reminder: DATABASE IS NOT BEING UPDATED")
 
    # Generate vars for overall html template
    pdb_html = "No structures or homology models were identified for this gene."
    disease1_html = "Insufficient variation was available for PathProx analysis."
    disease2_html = disease1_html
    ddG_html = "No ddG calculations were completed for this mutation."
 
    columns_for_html = ['Analyzable?','PDB Template','Template Identity','Distance to Boundary','Seq Start','Seq End','Trans Identity','PDB Pos','Residues','Resolution (PDB)','URL']
    
    if not tdf.empty:
        # pdb_html = tdf.to_html(float_format=lambda x: "%.2f"%x,formatters={'Seq Start': '${:,1.0f}'.format,'Seq End': '${:,1.0f}'.format})

        df_structure_report_human = df_structure_report.rename(columns=xlate_columns_dict)
        df_structure_report_human.style.set_properties(**{'text-align': 'center'})
        df_copy =  df_structure_report_human[columns_for_html].reset_index() # .to_html(justify='center') # ,formatters={'Seq Identity': '{:,.2f}'.format})
        # Really we should let django deal with the href
        # For now  we sin, and output the links unescaped...
        # First replace null links with a space
        # import pdb; pdb.set_trace()
        df_copy.URL.fillna('',inplace=True)
        df_copy['URL'] = df_copy['URL'].apply(lambda x: x if x else '')
        df_copy['structure_id'] = '<a href=\"' + df_copy['URL']+'">' + df_copy['structure_id'] + '</a>'  
        df_copy.drop(['URL'],axis=1,inplace=True)
        # pd.set_option('float_format', '%.2f')
        pdb_html = df_copy.to_html(justify='center',index=False,escape=False,float_format=lambda x: "%.2f"%x) 
    if not disease1_sum.empty:
      disease1_html = disease1_sum.to_html(float_format=lambda x: "%.2f"%x)
    if not disease2_sum.empty:
      disease2_html = disease2_sum.to_html(float_format=lambda x: "%.2f"%x)
    if not ddG_sum.empty:
      ddG_html = ddG_sum.to_html(float_format=lambda x: "%.2f"%x)
    
    
    graphicsLegend = "Pfam Domain Graphic unavailable"
    # Get Uniprot Graphics string or put in a blank
    domainGraphicsJSONstr = None
    if gathered_info['unp']: # Fetch the JSON - wait up to 30 seconds
      xmlHandle = Cxml2json(os.path.join(config_dict['interpro_dir'],'match_humanonly.xml'))
      if xmlHandle.isCanonical(gathered_info['unp']): # Then we need to manually create the xml ourselves
        # Go for the web link because this is a canonical UNP - however that _can_ fail.
        graphicsLegend = "Downloaded Pfam Domain Graphic for Canonical Isoform %s"%gathered_info['unp']
        LOGGER.info("Attempting download of Domain Graphics for %s from xfam.pfam"%gathered_info['unp'])
        domainGraphicsJSONstr = unp2PfamDomainGraphicString(gathered_info['unp'],30)
        if not domainGraphicsJSONstr:
          LOGGER.info("Download returned nothing")
          graphicsLegend = "Pfam Domain Graphic for Canonical Isoform %s"%gathered_info['unp']
      else:
        graphicsLegend = "Pfam Domain Graphic for Non-canonical Isoform %s"%gathered_info['unp']

      # _Either_ communications failure OR non-canonical isoform
      # So we create our own graphic from our local xml database
      if not domainGraphicsJSONstr: # _Either_ communications failure OR non-canonical (no communication)
        LOGGER.info("Creating Domain Graphic for %s from xml",gathered_info['unp'])
        nonCanonicalGraphicsJSON = xmlHandle.PFAMgraphicsJSONfromUnp(gathered_info['unp'])
        domainGraphicsJSONstr = json.dumps(nonCanonicalGraphicsJSON)
        # Not sure why - but graphics string from web is a list (inside outer [])
        # Removed this - did not work: domainGraphicsJSONstr = domainGraphicsJSONstr[0]
    
      if domainGraphicsJSONstr:
        domainGraphicsDict = json.loads(domainGraphicsJSONstr)
        # Add our mutation point of interest to this
        mutationSiteDict = {'colour': '#e469fe',
         'display': True,
         'headStyle': 'diamond',
         'lineColour': '#333333',
         'metadata': {'database': 'UDN Case',
                  'description': '%s'%mutation,
                  'start': mutation[1:-1],
                  'type': 'Mutation'} ,
         'residue': 'X',
         'start': mutation[1:-1],
         'type': 'UDN Mutation site',
         'v_align': 'top'}
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
    pfamGraphicsIframe_fname = "%s_%s_PfamGraphicsIframe.html"%(gathered_info['gene'],gathered_info['mutation'])
    website_filelist.append(os.path.join(web_dir,pfamGraphicsIframe_fname))
 
    # This data structure directly feeds the complex items in the psb_report.html template 
    # elements for project/case/gene/mutation/refseq were populated at top of this function
    template_vars.update({"title" : "PSB Structure Report - %s %s"%(gathered_info['unp'],gathered_info['mutation']),
                     "date"  : time.strftime("%Y-%m-%d"),
                     "pfamGraphicsIframe_fname" : pfamGraphicsIframe_fname,
                     "graphicsLegend" : graphicsLegend,
                     "psb_structure_results" : pdb_html,
                     "config_pathprox_dict"   : config_pathprox_dict,
                     "disease1_summary"   : disease1_html,
                     "disease2_summary"    : disease2_html,
                     "ddG_summary"    : ddG_html,
                     "DomainGraphicsJSON"  : domainGraphicsJSONstr,
                     "PfamResultTableHTML"  : PfamResultTableHTML,
                     "Structure_Detail"  : structures})
    if mutation:
      base_fname = "%s_%s_structure_final_report.%%s"%(gathered_info['gene'],mutation)
    else:
      base_fname = "%s_structure_final_report.%%s"%gathered_info['gene']

 
    """ Copying of the html system for each mutation was replaced with single case-wide copy
    # The html directory, and it's large set of supporting files, is located beneath the directory
    # where this file, psb_rep.py, is located 
    src =  os.path.join(os.path.dirname(os.path.realpath(__file__)),"html")
    dest = "%s/html"%variant_report_directory
    LOGGER.info("Copying supporting javascript and html from %s to %s"%(src,dest))

    
    try:
      shutil.rmtree(dest,ignore_errors=True)
      shutil.copytree(src,dest)
    except Exception as err:
      print err
      print("REMINDER: Need to exit")
      # sys.exit(1)


    
    src_file_list = [os.path.join(root,name) for root,dirs,files in os.walk(src) for name in files]
    
    
    dest_file_list = [os.path.join(root,name) for root,dirs,files in os.walk(dest) for name in files]
    if len(src_file_list) == len(dest_file_list):
      infostr = "Copy of %d html support files succeeded"%len(src_file_list)
      print infostr
      LOGGER.info(infostr)
      # print '\n'.join(src_file_list)
    else:
      LOGGER.critical("FAILURE to cp html support files (%d source files,%d dest files)"%(len(src_file_list),len(dest_file_list)))
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
    """
    
    env        = Environment(loader=FileSystemLoader(os.path.join(os.path.dirname(os.path.realpath(__file__)))))
    template   = env.get_template("psb_report.html")
    html_out = template.render(template_vars)

    templatePfamGraphics = env.get_template("html/pfamGraphicsIframeTemplate.html")
    htmlPfamGraphics = templatePfamGraphics.render(template_vars)

    save_cwd = os.getcwd()
    os.chdir(variant_report_directory)
    
    # WE ARE NOW OPERATING FROM ..../UDN/CaseName target directory
    LOGGER.info("Now working in %s",variant_report_directory)

    html_fname = base_fname%"html"
    website_filelist.append(os.path.join(web_dir,html_fname))
    with open(html_fname,"w") as html_f:
      html_f.write(html_out)

    with open(pfamGraphicsIframe_fname,"w") as html_f:
      html_f.write(htmlPfamGraphics)

    write_pdfs = None
    if write_pdfs:
      pdf_fname = base_fname%"pdf"
      wkhtmltopdf_fname = base_fname%"wkhtml.pdf"
      print("\nWriting final reports to %s directory:\n  %s\n  %s\n  %s\n  %s\n"%(variant_report_directory,pdf_fname,wkhtmltopdf_fname,html_fname,pfamGraphicsIframe_fname))
    
      LOGGER.warning("Temporarily disabling all logging prior to calling weasyprint write_pdf()")
    
      logging.disable(sys.maxsize)
      HTML(string=html_out).write_pdf(pdf_fname,stylesheets=["../html/css/typography.css"])#,CSS(string="@page { size: A3 landscape; }")])
      logging.disable(logging.NOTSET)
      LOGGER.warning("weasyprint write_pdf() has returned.  Restoring logging")
         
      
      # Write out another .pdf using the wkhtmltopdf tool
    
      wkhtml_command_list = ["wkhtmltopdf","--print-media-type","--no-outline","--minimum-font-size","12",html_fname,wkhtmltopdf_fname]
      process = subprocess.Popen(wkhtml_command_list, stdout = subprocess.PIPE,stderr = subprocess.PIPE)
      (output,err) = process.communicate()
    
      # So unfortunately, the stderr from wkhtmltopdf is full of ==================\r things, which are a mess for the LOGGER.  Get rid of stuff that ends in \r  Keep all other lines and ending in \r
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
        print("Unable to complete %s exitstatus=%d due to error %s\n  output: %s\n"%(str(wkhtml_command_list),process.returncode,err_legit,output))
        LOGGER.warning("wkhtmltopdf stderr:\n%s",err_legit)
      else:
        LOGGER.info("wkhtmltopdf stderr:\n%s",err_legit)
      
    # WE HAVE NOW RETURNED TO the psb_pipeline/bin directory
    os.chdir(save_cwd);
    # Close out the local log file for this mutation
    local_fh.flush()
    local_fh.close()
    LOGGER.removeHandler(local_fh)
    gathered_info['variant_report_directory'] = os.path.basename(variant_report_directory)
    gathered_info['html_fname'] = html_fname
    return gathered_info # End of function report_one_mutation()
  
# Main logic here.  A period in the argument means the user wants to launch one mutation only,
# Whether one mutation, or the more typical set, we must have html/css/javascript statuc
# resources property installed in the destination tree
copy_html_css_javascript()

# directly from a single mutation output file of psb_plan.py
if oneMutationOnly:
  if args.slurm:
    LOGGER.warning("--slurm option is ignored when only one report is requested")
  template_vars = report_one_mutation(args.projectORstructures,args.workstatus)
  if template_vars: # The usual case, we had some Pathprox outputs
    print("Single mutation report saved to %s/.pdf/.wkhtml.pdf"%template_vars['html_fname'])
  else:
    print("Due to lack of pathprox outputs, no html (or pdf) reports were created from %s"%args.workstatus)
else:
  udn_csv_filename = os.path.join(collaboration_dir,"%s_missense.csv"%args.projectORstructures) # The argument is an entire project UDN124356
  msg = "Retrieving project mutations from %s"%udn_csv_filename
  if not infoLogging:
    print (msg)
  LOGGER.info(msg)
  df_all_mutations = pd.read_csv(udn_csv_filename,sep=',',index_col = None,keep_default_na=False,encoding='utf8',comment='#',skipinitialspace=True)
  df_all_mutations.fillna('NA',inplace=True)

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
    print (msg)
  LOGGER.info(msg)

  mutation_summaries = []
  slurm_array = []
  for index,row in df_all_mutations.iterrows():
    msg = "%s %-10s %-10s %-6s"%("Generating slurm entry for" if args.slurm else "Reporting on", row['gene'],row['refseq'],row['mutation'])
    if not infoLogging:
      print (msg)
    LOGGER.info(msg)
    if 'RefSeqNotFound_UsingGeneOnly' in row['refseq']:
      row['refseq'] = 'NA'
    mutation_dir = os.path.join(collaboration_dir,"%s_%s_%s"%(row['gene'],row['refseq'],row['mutation']))
    if not os.path.exists(mutation_dir):  # python 3 has exist_ok parameter... 
      LOGGER.critical("The specific mutation directory %s should have been created by psb_status.py.  Fatal problem.",mutation_dir)
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
        if not gathered_info['Error'] and gathered_info['variant_report_directory']:
          gathered_info['html_fname'] = os.path.join(gathered_info['variant_report_directory'],gathered_info['html_fname'])
          msg = "%s %s %s report saved to %s/.pdf/.wkhtml.pdf"%(row['gene'],row['refseq'],row['mutation'],gathered_info['html_fname'])
          if not infoLogging:
            print (msg)
          LOGGER.info(msg)
        gathered_info['Gene Mutation'] = "%-9.9s %-10.10s %-10.10s %-7.7s"%(row['gene'],row['refseq'],gathered_info['unp'],row['mutation'])
        gathered_info['Gene Hit Generic'] = ''
        mutation_summaries.append(gathered_info)
      else:
        msg = "Due to lack of pathprox or ddG outputs, %s %s %s has no html (or pdf) report"%(row['gene'],row['refseq'],row['mutation'])
        if not infoLogging:
          print (msg)
        LOGGER.info(msg)

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

      slurmDict = slurmParametersReport
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
      LOGGER.info(msg);
    else:
      print(msg)
       
  if mutation_summaries: # Excellent, we have a home page report for many variants
    # If there are DigenicInteraction .svg files in the right place, integrate them
    digenic_graphics_types = ['digenic_score_only','digenic_details','sys_bio_details']
    digenic_graphics_score_only_filename = None
    found_graphics_files = {}
    for digenic_graphic_type in digenic_graphics_types:
      digenic_graphic_filename = os.path.join('casewide','DigenicAnalysis',
        '%s_all_gene_pairs_%s.svg'%(args.projectORstructures,digenic_graphic_type))
      
      if os.path.exists(os.path.join(collaboration_dir,digenic_graphic_filename)):
        found_graphics_files[digenic_graphic_type] = digenic_graphic_filename
        website_filelist.append(os.path.join('.',digenic_graphic_filename))

    env = Environment(loader=FileSystemLoader(os.path.join(os.path.dirname(os.path.realpath(__file__)))))
    if len(found_graphics_files) == 3:
      LOGGER.info("Integrating DigenicAnalysis svg graphics")
      digenic_graphics_score_only_filename = os.path.join('.',found_graphics_files['digenic_score_only'])
      template   = env.get_template("html/DigenicInteractionsReportTemplate.html")
      # Probably a goof - but the template file restates the full graphics filenames
      html_out = template.render({'case': args.projectORstructures})
      digenic_graphics_html_filename = os.path.join(collaboration_dir,'casewide','DigenicGraphics.html') # The argument is an entire project UDN124356
      with open(digenic_graphics_html_filename,"w") as f:
        f.write(html_out)
      website_filelist.append(digenic_graphics_html_filename)
    else:
      LOGGER.warning("Digenic Analysis svg files are missing.  Did you run DigenicAnalysis")


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
    template   = env.get_template("case_report_template.html")
    # print html_table
    final_gathered_info = {'mutation_summaries': mutation_summaries, 
             'firstGeneTable': html_table_generic,
             'firstGeneReport': html_report_generic,
             'secondGeneTable': html_table_familial,
             'secondGeneReport': html_report_familial,
             'case': args.projectORstructures, 
             "date"  : time.strftime("%Y-%m-%d"),
             'disease1_variant_short_description':  config_pathprox_dict['disease1_variant_short_description'],
             'disease2_variant_short_description':  config_pathprox_dict['disease2_variant_short_description'],
             'digenic_graphics_score_only_filename': digenic_graphics_score_only_filename
             }
    html_out = template.render(final_gathered_info)
    case_summary_filename = os.path.join(collaboration_dir,"%s.html"%args.projectORstructures) # The argument is an entire project UDN124356
    website_filelist.append(case_summary_filename)
    with open(case_summary_filename,"w") as f:
      f.write(html_out)
    lastmsg = "The case summary report is: " + case_summary_filename
    case_summary_json = os.path.join(collaboration_dir,"%s.json"%args.projectORstructures) # The argument is an entire project UDN124356
    with open(case_summary_json,'w') as fp:
        json.dump(final_gathered_info,fp)
  else:
    lastmsg = "No mutation summaries - bummer"

  print ("")
  print ("")

  if args.verbose:
    LOGGER.info(lastmsg);
  else:
    print(lastmsg)

  if website_filelist:
    try:
      os.remove('index.html') # First time through will gen an exception.  That's A-OK
    except:
      pass
    os.symlink(case_summary_filename,'index.html')
    website_filelist.append('index.html')

    website_filelist_filename = os.path.join(collaboration_dir,"%s_website_files.list"%args.projectORstructures) # The argument is an entire project UDN124356
    with open(website_filelist_filename,"w") as f:
      f.write('\n'.join((os.path.relpath(
              website_file, # os.path.realpath(website_file),
                (os.path.realpath(os.path.join(config_dict['output_rootdir'],config_dict['collaboration']))))
                  for website_file in website_filelist)))
    LOGGER.info("A filelist for creating a website is in %s",website_filelist_filename)
    website_zip_filename = "%s.zip"%args.projectORstructures
    pkzip_maker = 'rm -f %s; cd ..; cat %s | zip -r@ %s > %s.stdout; cd -'%(website_zip_filename,
         os.path.join(args.projectORstructures,website_filelist_filename),
         os.path.join(args.projectORstructures,website_zip_filename),
         os.path.join(args.projectORstructures,website_zip_filename))
    LOGGER.info("Executing: %s",pkzip_maker)
    subprocess.call(pkzip_maker,shell=True)

    website_tar_filename = "%s.tar.gz"%args.projectORstructures
    tar_maker = 'rm -f %s; cd ..; tar cvzf %s --files-from %s --mode=\'a+rX,go-w,u+w\' > %s.stdout; cd -'%(website_tar_filename,
         os.path.join(args.projectORstructures,website_tar_filename),
         os.path.join(args.projectORstructures,website_filelist_filename),
         os.path.join(args.projectORstructures,website_tar_filename))
    LOGGER.info("Executing: %s",tar_maker)
    subprocess.call(tar_maker,shell=True)

    print ("Compressed website files are in %s and %s"%(website_zip_filename,website_tar_filename))
 
sys.exit(0)
