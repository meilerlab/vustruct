#!/usr/bin/env python
#
# Project        : PSB Pipeline
# Filename       : psb_monitor.py
# Authors        : Chris Moth and R. Michael Sivley
#project Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : chris.moth@vanderbilt.edu
# Date           : 2018-01-22
# Description    : Report job progress on a SLURM cluster 
#                : (output from psb_launch.py or prior run of psb_monitor.py)
#=============================================================================#

"""\
Given either a project name on the command line, or a single mutation workstatus.csv
file previously output from psb_launch.py (or prior run of psb_monitor.py),
this script will report on the progress of the launched jobs for each mutation.

psb_monitor.py interrogates the .../status and .../info directories of each

It could also be updated to  inquire status of slurm's "scontrol show job" feature
HOWEVER - this is very slow - and ACCRE staff have requested we not do this

The workstatus.csv file for each mutation is updated.

Configuration supplied by:
   -c global.config file
   -u user.config overrides
"""

print("%s: Pipeline monitor for launched jobs.  -h for detailed help."%__file__)

import logging,os,pwd,sys
import time, datetime
import argparse,configparser
import pprint
from logging.handlers import RotatingFileHandler
from logging import handlers
sh = logging.StreamHandler()
LOGGER = logging.getLogger()
LOGGER.addHandler(sh)

# Now that we've added streamHandler, basicConfig will not add another handler (important!)
log_format_string = '%(asctime)s %(levelname)-4s [%(filename)16s:%(lineno)d] %(message)s'
date_format_string = '%H:%M:%S'
log_formatter = logging.Formatter(log_format_string,date_format_string)

LOGGER.setLevel(logging.DEBUG)
sh.setLevel(logging.WARNING)
sh.setFormatter(log_formatter)

import pandas as pd
import numpy as np
import hashlib
from jinja2 import Environment, FileSystemLoader

from slurm import SlurmJob,slurm_scontrol_show_job,slurm_jobstate_isfinished

cmdline_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)

from psb_shared import psb_config

cmdline_parser = psb_config.create_default_argument_parser(__doc__,os.path.dirname(os.path.dirname(__file__)))
cmdline_parser.add_argument("projectORworkstatus",type=str,
                   help="Project ID (ex. UDN123456) or single workstatus.csv output file from psb_launch.py (or previous psb_monitor.py)  Example: ......$UDN/UDN123456/GeneName_NM_12345.1_S123A_workstatus.csv",
                   default = os.path.basename(os.getcwd()),nargs='?')

args,remaining_argv = cmdline_parser.parse_known_args()

infoLogging = False

if args.debug:
  infoLogging = True
  sh.setLevel(logging.DEBUG)
elif args.verbose:
  infoLogging = True
  sh.setLevel(logging.INFO)
# If neither option, then logging.WARNING (set at top of this code) will prevail for stdout


# A period in the argument means the user wants to monitor one mutation only,
# directly from a single mutation output file of psb_launch.py
oneMutationOnly = ('.' in args.projectORworkstatus and os.path.isfile(args.projectORworkstatus))

required_config_items = ['output_rootdir','collaboration']

config,config_dict = psb_config.read_config_files(args,required_config_items)

# print "Command Line Arguments"
LOGGER.info("Command Line Arguments:\n%s"%pprint.pformat(vars(args)))

# The collaboration_dir is the master directory for the case, i.e. for one patient
# Example: /dors/capra_lab/projects/psb_collab/UDN/UDN532183
udn_root_directory = os.path.join(config_dict['output_rootdir'],config_dict['collaboration'])
if oneMutationOnly:
  collaboration_dir = os.path.dirname(os.path.dirname(os.path.abspath(args.projectORworkstatus)))
else:
  collaboration_dir = os.path.join(udn_root_directory,args.projectORworkstatus)

if not os.path.exists(collaboration_dir):  # python 3 has exist_ok parameter...
  logging.error("%s not found.  It should have been created by psb_plan.py"%collaboration_dir)
  sys.exit(1)

# collaboration_log_dir = os.path.join(collaboration_dir,"log")
# if not os.path.exists(collaboration_log_dir):  # python 3 has exist_ok parameter...
#  logging.error("%s not found.  It should have been created by psb_plan.py"%collaboration_log_dir)
#  sys.exit(1)

def monitor_one_mutation(workstatus):
  slurmInfoColumns = (['scontrolTimestamp','JobState','ExitCode','RunTime','TimeLimit',
                     'SubmitTime','EligibleTime','StartTime','EndTime',
                     'NodeList','BatchHost','NumNodes','NumCPUs','NumTasks','StdErr','StdOut','StdIn'])
  # set the type of all the slurmInfoColumns as string
  type_dict = dict((c,str) for c in slurmInfoColumns)
  # Load the status of jobs that was created by psb_launch.py (or prior run of psb_monitor.py)
  # Load the schedule of jobs that was created by psb_launch.py (or previous run of psb_monitor.py)
  # keep_default_na causes empty strings to come is as such, and non pesky nan floats
  df_all_jobs_status = pd.read_csv(workstatus,'\t',dtype=str,keep_default_na=False)
  LOGGER.info("%d rows read from workstatus file %s"%(len(df_all_jobs_status),workstatus))
 
  df_all_jobs_status.set_index('uniquekey',inplace=True)
  if len(df_all_jobs_status) < 1:
    LOGGER.warning("No rows in work status file %s.  No jobs to monitor."%workstatus)
    return pd.DataFrame()

  df_all_jobs_original_status = df_all_jobs_status.copy()


  """
  These variables are not relevant to pdb_monitor.py: monitor_one_mutation
  row0 = df_all_jobs_status.iloc[0]
  gene = row0['gene']
  refseq = row0['refseq']
  mutation = row0['mutation']
  mutation_dir = os.path.join(collaboration_dir,"%s_%s_%s"%(gene,refseq,mutation))
  if not os.path.exists(mutation_dir):  # python 3 has exist_ok parameter...
    logging.error("%s not found.  It should have been created by psb_plan.py"%mutation_dir)
    sys.exit(1)
  mutation_log_dir = mutation_dir # os.path.join(mutation_dir,"log")
  if not os.path.exists(mutation_log_dir):  # python 3 has exist_ok parameter...
    logging.error("%s not found.  It should have been created by psb_plan.py"%mutation__log_dir)
    sys.exit(1)
  """


  # Life is easiest if we do not set the index until we are ready to update the df_all_jobs_status at
  # the end of the monitor check-in.
  # df_job updates is populated by not only slurm's scontrol command, but also the status directory
  # that is filled in with files 'complete', 'progress', and 'info'
  # When we run 'scontrol show job' we get a _LOT_ of things back.
  # For now grab lots of these things in a new dataframe
  df_job_updates = pd.DataFrame(columns=['jobid','uniquekey'] + slurmInfoColumns)
  
  df_job_updates[slurmInfoColumns] = df_job_updates[slurmInfoColumns].astype(str)
  
  
  for index,row in df_all_jobs_status.iterrows():
    # If we already have the ExitCode of the slurm job from prior run of this module
    # DO NOT attempt to further update.  There can be no new news after exit code arrives
    uniquekey = index
    jobid = row['jobid']
    interesting_info = {}
    interesting_info['uniquekey'] = uniquekey
    interesting_info['jobid'] = jobid
    if ('ExitCode' in row) and row['ExitCode'] and (len(row['ExitCode']) >= 1) and (int(row['ExitCode']) == 0):
      LOGGER.info("%15s:%-20s Exit Code %s recorded previously"%(jobid,uniquekey,str(row['ExitCode'])))
      continue
    # The most reliable source of Exit=0 is the arrival of a completed file in the status director
    progress_file_found = False
    statusdir = "%s/%s/status"%(row['outdir'],row['flavor'])
    try:
      with open("%s/progress"%statusdir) as progressfile:
        s = progressfile.read().replace('\n','')
        if s and len(s):
          interesting_info['jobprogress'] = s
          progress_file_found = True
    except Exception as ex:
      pass
 
    info_file_found = False
    try:
      with open("%s/info"%statusdir) as infofile:
        s = infofile.read().replace('\n','')
        if s and len(s):
          interesting_info['jobinfo'] = s
          info_file_found = True
    except Exception as ex:
      pass


    if os.path.exists("%s/complete"%statusdir):
      LOGGER.info("%15s:%-20s Found `complete` file.  Recording success"%(jobid,uniquekey))
      interesting_info['ExitCode'] = '0'
      interesting_info['jobprogress'] = 'Completed' # For now, we will overwrite memory on happy exit - need to think
      interesting_info['jobinfo'] = 'Completed'
    elif os.path.exists("%s/FAILED"%statusdir):
      LOGGER.info("%15s:%-20s Found `FAILED` file.  Recording failure"%(jobid,uniquekey))
      interesting_info['ExitCode'] = '1'
    else: # Try to grab updated progress info from the output files.  No worries if nothing there
      if not (progress_file_found or info_file_found):
        interesting_info['jobprogress'] = "No status updates.  Inspect %s"%statusdir
        LOGGER.info("%s",interesting_info['jobprogress'])
        interesting_info['jobinfo'] = 'No status updates'
        
      pass # We no longer call "scontrol show job" - it just takes way too long      
      """scontrol_info = slurm_scontrol_show_job(jobid)
      if (not scontrol_info) or (scontrol_info.get('JobState','') == 'UNKNOWN') :
        LOGGER.info("%15s:%-20s Nothing returned from scontrol"%(jobid,uniquekey))
        continue
      # We don't want _everything_ from the scontrol call - but as much as we can get
      else:
        interesting_info.update({k: scontrol_info.get(k,'') for k in slurmInfoColumns})
        LOGGER.info("%15s:%-20s JobState=%s"%(jobid,uniquekey,interesting_info['JobState']))
        interesting_info['scontrolTimestamp'] = datetime.datetime.now()
        # Don't update the ExitCode _unless_ the job is clearly finished
        if not slurm_jobstate_isfinished(interesting_info['JobState']):
          interesting_info['ExitCode'] = ''
        elif interesting_info['JobState'] == 'FAILED':
          interesting_info['ExitCode'] = '1'
        elif ':' in interesting_info.get('ExitCode',''):  # Slurm ExitCode format is exit:signal - get rid of signal
          interesting_info['ExitCode'] = interesting_info['ExitCode'].split(':')[0]
        info_as_series = pd.Series(data=interesting_info)
        df_job_updates = df_job_updates.append(info_as_series,ignore_index=True)
      """ 
    
    # print df_job_updates
    info_as_series = pd.Series(data=interesting_info)
    df_job_updates = df_job_updates.append(info_as_series,ignore_index=True)
  
  previous_workstatus_filename = workstatus + ".previous"
  # print "Renaming current status file to %s"%previous_workstatus_filename
  try:
    os.remove(previous_workstatus_filename)
  except Exception as ex:
    pass
   
  
  
  # If the slurmInfoColumns are not alreadt in df_all_jobs_status, then add those columns 
  # will add them with NaN initializers.  If these columns are already there, do nothing
  # (I've tried concat (only works pre 0.20 pandas) and other things - and have 
  #  given up on elegant solutions)
  for column in slurmInfoColumns:
    if column not in df_all_jobs_status:
      df_all_jobs_status[column] = None
  
  # Old idea failed on different versions: df_all_jobs_status = pd.concat([df_all_jobs_status,pd.DataFrame(columns=slurmInfoColumns)])
  # Can't work because done at top: df_all_jobs_status.set_index(['uniquekey','jobid'],inplace=True)
 
  # This one here should work because we've not set the index before on it
  df_job_updates = df_job_updates.set_index('uniquekey')
  
  # Take all the non-NA values from df_job_updates, and modify df_all_jobs_status to include them
  df_all_jobs_original_status = df_all_jobs_status.copy()
  # import pdb; pdb.set_trace()
  df_all_jobs_status.update(df_job_updates)
  df_incomplete = df_all_jobs_status[df_all_jobs_status['ExitCode'] != '0']

  # Are there any changes to what we know about the status of the jobs?
  if df_all_jobs_original_status.to_msgpack() == df_all_jobs_status.to_msgpack():
    if len(df_incomplete) == 0:
      # This happens out to right of earlier printing...
      print("   All %d jobs completed successfully"%len(df_all_jobs_status))
      return
    if "casewide" in workstatus:
      print("     No updates to status of casewide jobs")
    else:
      print("     No updates to status of jobs for this mutation")
  else:
    os.rename(workstatus,previous_workstatus_filename)
  
    print("\nRecording all updates to %s"%workstatus)
    try:
      df_all_jobs_status.to_csv(workstatus,sep='\t',index=True)
    except Exception as ex:
      # If we cannot save the new csv file, then we must (attempt to!) restore the old one!
      os.remove(workstatus)
      os.rename(workstatus,previous_workstatus_filename)
      LOGGER.exception("Serious failure saving new updated file %s.  Prior file restored"%workstatus)
      sys.exit(1)
  
  if len(df_incomplete) == 0:
    print("All %d jobs completed successfully"%len(df_all_jobs_status))
  else:
    print("%d of %d jobs still incomplete:"%(len(df_incomplete),len(df_all_jobs_status)))
    print("%15s:%-20s  %s"%('Jobid','Flavor',"Info"))
    for index,row in df_incomplete.iterrows():
      infostring = ""
      if 'jobinfo' in row and len(row['jobinfo']) > 1:
        infostring = row['jobinfo']
      elif 'JobState' in row and len(row['JobState']) > 1:
        infostring = row['JobState']
  
      print("%15s:%-20s    %s"%(row['jobid'],index,infostring))
  
# Main logic here.  A period in the argument means the user wants to launch one mutation only,
# directly from a single mutation output file of psb_plan.py
if oneMutationOnly:
  monitor_one_mutation(args.projectORworkstatus)  #  The argument is a complete workstatus filename
else:
  udn_csv_filename = os.path.join(collaboration_dir,"%s_missense.csv"%args.projectORworkstatus) # The argument is an entire project UDN124356
  print("Retrieving project mutations from %s"%udn_csv_filename)
  df_all_mutations = pd.read_csv(udn_csv_filename,sep=',',index_col = None,keep_default_na=False,encoding='utf8',comment='#',skipinitialspace=True)
  print("Monitoring all jobs for %d mutations"%len(df_all_mutations))
  df_all_mutations.fillna('NA',inplace=True)
  for index,row in df_all_mutations.iterrows():
    # print without a newline - monitor_one_mutation will add one
    print("%d of %d: %-10s %-10s %-6s"%(index+1,len(df_all_mutations),row['gene'],row['refseq'],row['mutation']), end=' ')
    if 'RefSeqNotFound_UsingGeneOnly' in row['refseq']:
      row['refseq'] = 'NA'
    mutation_dir = os.path.join(collaboration_dir,"%s_%s_%s"%(row['gene'],row['refseq'],row['mutation']))
    if not os.path.exists(mutation_dir):  # python 3 has exist_ok parameter... 
      logging.critical("The specific mutation directory %s should have been created by psb_status.py.  Fatal problem."%mutation_dir)
      sys.exit(1)
       
    workstatus_filename = "%s/%s_%s_%s_workstatus.csv"%(mutation_dir,row['gene'],row['refseq'],row['mutation'])
    # print workstatus_filename
    monitor_one_mutation(workstatus_filename)  #  The argument is a complete workstatus filename
    print("-" * 80)
  # Finally monitor the job(s) under the casewide banner (Gene Dictionary creation)
  workstatus_filename = "casewide/casewide_workstatus.csv"
  if os.path.exists(workstatus_filename):
    # print workstatus_filename
    print("Casewide jobs.....                      ", end=' ')
    monitor_one_mutation(workstatus_filename)  #  The argument is a complete workstatus filename
  else:
    print("No casewide work (no %s) to monitor"%workstatus_filename)
  print("-" * 80)
