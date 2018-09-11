#!/usr/bin/env python2.7
#
# Project        : PSB Pipeline
# Filename       : psb_launch.py
# Authors        : Chris Moth and R. Michael Sivley
#project Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : chris.moth@vanderbilt.edu
# Date           : 2018-01-22
# Description    : Launch jobs specified in the input .csv file
#                : (output from psb_plan.py)
#=============================================================================#

"""\
Launch jobs on a SLURM cluster given 
   -c global.config file
   -u user.config overrides
   --verbose Log all info to screen in addition to the logfile
   --debug As --verbose, but additionally logs python 'debug' level information
   workplan.csv previously output from psb_plan.py

   --relaunch  Launch all jobs for mutation(s) regardless of their workstatus after a prior launch
               Without this flag (normal case), then jobs that have never been launched before,
               or which failed in prior launch will be launched
"""

print "%s: Pipeline launcher.  Run after psb_plan.py.   -h for detailed help"%__file__

import logging,os,pwd,sys,grp,stat
import time, datetime
import argparse,ConfigParser
import pprint
from logging.handlers import RotatingFileHandler
from logging import handlers

capra_group = grp.getgrnam('capra_lab').gr_gid

def set_capra_group_sticky(dirname):
  try:
    os.chown(dirname, -1, capra_group)
  except:
    pass

  # Setting the sticky bit on directories also fantastic
  try:
    os.chmod(dirname, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IXGRP | stat.S_ISGID);
  except:
    pass

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
# from jinja2 import Environment, FileSystemLoader

from slurm import slurm_submit

default_global_config=os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))),"config","global.config")

cmdline_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)

cmdline_parser.add_argument("-c","--config",
help="PDBMap configuration profile for database access", required=False,metavar="FILE",default=default_global_config)
cmdline_parser.add_argument("-u","--userconfig",
help="User specific settings and configuration profile overrides", required=True,metavar="FILE")
cmdline_parser.add_argument("projectORworkplan",type=str,help="Project ID (ex. 123456) or single output file from psb_plan.py  Example: ....../UDN/UDN123456/GeneName_NM_12345.1_S123A_workplan.csv")
cmdline_parser.add_argument("-v","--verbose",
help="Include routine info log entries on stderr", action = "store_true")
cmdline_parser.add_argument("-r","--relaunch",
help="Relaunch jobs, ignoring any results in a previously generated workstatus.csv file", action = "store_true")
cmdline_parser.add_argument("-d","--debug",
help="Include routine info AND 'debug' log entries on stderr", action = "store_true")

args,remaining_argv = cmdline_parser.parse_known_args()

infoLogging = False

if args.debug:
  infoLogging = True
  sh.setLevel(logging.DEBUG)
elif args.verbose:
  infoLogging = True
  sh.setLevel(logging.INFO)
# If neither option, then logging.WARNING (set at top of this code) will prevail for stdout

# A period in the argument means the user wants to launch one mutation only,
# directly from a single mutation output file of psb_plan.py
oneMutationOnly = '.' in args.projectORworkplan


config = ConfigParser.SafeConfigParser(allow_no_value=True)
config.read([args.config])
config_dict = dict(config.items("Genome_PDB_Mapper")) # item() returns a list of (name, value) pairs
try:
  globalSlurmParametersAll = dict(config.items("SlurmParametersAll"))
except Exception as ex:
  globalSlurmParametersAll = {}

try:
  globalSlurmParametersPathProx = dict(config.items("SlurmParametersPathProx"))
except Exception as ex:
  globalSlurmParametersPathProx = {}

try:
  globalSlurmParametersUDNStructure = dict(config.items("SlurmParametersUDNStructure"))
except Exception as ex:
  globalSlurmParametersUDNStructure = {}

try:
  globalSlurmParametersUDNSequence = dict(config.items("SlurmParametersUDNSequence"))
except Exception as ex:
  globalSlurmParametersUDNSequence = {}

try:
  globalSlurmParametersMakeGeneDictionaries = dict(config.items("SlurmParametersMakeGeneDictionaries"))
except Exception as ex:
  globalSlurmParametersMakeGeneDictionaries = {}

# Init all the user config dictionaries to empty

if (args.userconfig):
  userconfig = ConfigParser.SafeConfigParser() 
  userconfig.read([args.userconfig])
  config_dict.update(dict(userconfig.items("UserSpecific"))) # item() returns a list of (name, value) pairs

  # import pdb; pdb.set_trace()
  try:
    userSlurmParametersAll = dict(userconfig.items("SlurmParametersAll"))
  except:
    userSlurmParametersAll = {}

  try:
    userSlurmParametersPathProx = dict(userconfig.items("SlurmParametersPathProx"))
  except:
    userSlurmParametersPathProx = {}

  try:
    userSlurmParametersUDNSequence = dict(userconfig.items("SlurmParametersUDNSequence"))
  except:
    userSlurmParametersUDNSequence = {}

  try:
    userSlurmParametersUDNStructure = dict(userconfig.items("SlurmParametersUDNStructure"))
  except:
    userSlurmParametersUDNStructure = {}

  try:
    userSlurmParametersMakeGeneDictionaries = dict(userconfig.items("SlurmParametersMakeGeneDictionaries"))
  except:
    userSlurmParametersMakeGeneDictionaries = {}

# slurm settings for pathprox2.py with --clinvar
SlurmParametersPathProx = {}
SlurmParametersPathProx.update(globalSlurmParametersAll)
SlurmParametersPathProx.update(userSlurmParametersAll)
SlurmParametersPathProx.update(globalSlurmParametersPathProx)
SlurmParametersPathProx.update(userSlurmParametersPathProx)

# for udn_sequence.py in sequence analysis mode
SlurmParametersUDNSequence = {}
SlurmParametersUDNSequence.update(globalSlurmParametersAll)
SlurmParametersUDNSequence.update(userSlurmParametersAll)
SlurmParametersUDNSequence.update(globalSlurmParametersUDNSequence)
SlurmParametersUDNSequence.update(userSlurmParametersUDNSequence)

# for udn_sequence.py in ddG mode
SlurmParametersUDNStructure = {}
SlurmParametersUDNStructure.update(globalSlurmParametersAll)
SlurmParametersUDNStructure.update(userSlurmParametersAll)
SlurmParametersUDNStructure.update(globalSlurmParametersUDNStructure)  
SlurmParametersUDNStructure.update(userSlurmParametersUDNStructure)  

# for psb_genedicts.py
SlurmParametersMakeGeneDictionaries = {}
SlurmParametersMakeGeneDictionaries.update(globalSlurmParametersAll)
SlurmParametersMakeGeneDictionaries.update(userSlurmParametersAll)
SlurmParametersMakeGeneDictionaries.update(globalSlurmParametersMakeGeneDictionaries)  
SlurmParametersMakeGeneDictionaries.update(userSlurmParametersMakeGeneDictionaries)  

args,remaining_argv = cmdline_parser.parse_known_args()

logger.info("Command: %s"%' '.join(sys.argv))

slurm_required_settings = ['mem', 'account', 'ntasks', 'time']
for slurmSettingsDesc,slurmSettings in [('SlurmParametersPathProx',SlurmParametersPathProx), ('SlurmParametersUDNStructure',SlurmParametersUDNStructure),('SlurmParametersUDNSequence',SlurmParametersUDNSequence),('SlurmParametersMakeGeneDictionaries',SlurmParametersMakeGeneDictionaries)]:
  for req in slurm_required_settings:
    if not req in slurmSettings:
      logger.error("Can't launch jobs because you have not provided the required\nslurm setting [%s] for section %s (or SlurmParametersAll)\n in either file %s or %s\n"%
         (req,slurmSettingsDesc,args.config,args.userconfig))
      sys.exit(1)

# print "Command Line Arguments"
logger.info("Command Line Arguments:\n%s"%pprint.pformat(vars(args)))
logger.info("SlurmParametersPathProx:\n%s"%pprint.pformat(SlurmParametersPathProx))
logger.info("SlurmParametersUDNSequence:\n%s"%pprint.pformat(SlurmParametersUDNSequence))
logger.info("SlurmParametersUDNStructure:\n%s"%pprint.pformat(SlurmParametersUDNStructure))
logger.info("SlurmParametersMakeGeneDictionaries:\n%s"%pprint.pformat(SlurmParametersMakeGeneDictionaries))

# The collaboration_dir is the master directory for the case, i.e. for one patient
# Example: /dors/capra_lab/projects/psb_collab/UDN/UDN532183
udn_root_directory = os.path.join(config_dict['output_rootdir'],config_dict['collaboration'])
if oneMutationOnly:
  collaboration_dir = os.path.dirname(os.path.dirname(args.projectORworkplan))
else:
  collaboration_dir = os.path.join(udn_root_directory,args.projectORworkplan)

if not os.path.exists(collaboration_dir):  # python 3 has exist_ok parameter...
  os.makedirs(collaboration_dir)

set_capra_group_sticky(collaboration_dir)
collaboration_log_dir = os.path.join(collaboration_dir,"log")
if not os.path.exists(collaboration_log_dir):  # python 3 has exist_ok parameter...
  logger.error("%s not found.  It should have been created by psb_plan.py"%collaboration_log_dir)
  sys.exit(1)

def launch_one_mutation(workplan):
  # Load the schedule of jobs that was created by psb_plan.py
  df_all_jobs = pd.read_csv(workplan,sep='\t',keep_default_na = False,na_filter=False)
  logger.info("%d rows read from work plan file %s"%(len(df_all_jobs),workplan))

  if len(df_all_jobs) < 1:
    logger.warning("No rows in work plan file %s.  No jobs will be launched."%workplan)
    return pd.DataFrame(),workplan

  # Let's take care to not relaunch jobs that are already running, or complete, unless --relaunch flag is active
  workstatus_filename = workplan.replace('workplan','workstatus')

  df_prior_success = pd.DataFrame()
  if not args.relaunch:
    try:
      # import pdb; pdb.set_trace()
      df_workstatus = pd.read_csv(workstatus_filename,sep='\t',keep_default_na=False,na_filter=False,dtype=str)
      df_prior_success = df_workstatus[df_workstatus['ExitCode'] == '0'].copy()
      df_prior_success.set_index('uniquekey',inplace = True)

    except:
      # It's A-OK to not have a workstatus
      pass

  row0 = df_all_jobs.iloc[0]
  gene = row0['gene']
  refseq = row0['refseq']
  mutation = row0['mutation']

  geneRefseqMutation_OR_casewideString = 'casewide'
  if gene != geneRefseqMutation_OR_casewideString: # If gene is not casewide run, then fill it in specifically
    geneRefseqMutation_OR_casewideString = "%s_%s_%s"%(gene,refseq,mutation) # This is usual case for the gene by gene mutation run sets

  mutation_dir = os.path.join(collaboration_dir,geneRefseqMutation_OR_casewideString)
  if not os.path.exists(mutation_dir):  # python 3 has exist_ok parameter...
    os.makedirs(mutation_dir)
  set_capra_group_sticky(mutation_dir)
  mutation_log_dir = mutation_dir # os.path.join(mutation_dir,"log")
  if not os.path.exists(mutation_log_dir):  # python 3 has exist_ok parameter...
    os.makedirs(mutation_log_dir)
  set_capra_group_sticky(mutation_log_dir)

  # pw_name = pwd.getpwuid( os.getuid() ).pw_name # example jsheehaj or mothcw
  log_filename = os.path.join(mutation_log_dir,"psb_launch.log")


  if oneMutationOnly:
    sys.stderr.write("psb_launch log file is %s\n"%log_filename)
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
  
 
  # We are going to launch the Clinvar and COSMIC pathprox in a single slurm array
  # udn structure (ddG reports) and udn_sequence analysis each have unique requirements
  def build_slurm_dir(flavor):
    return os.path.join(mutation_dir,"slurm",("PathProx" if "PathProx" in flavor else flavor))
  def build_slurm_dir_from_row(row):
    return build_slurm_dir(row['flavor'])

  df_all_jobs['slurm_directory'] = df_all_jobs.apply(build_slurm_dir_from_row, axis=1)
  
  def build_gene_refseq_mutation_flavor(row):
    if row['gene'] == 'casewide':
      return row['flavor']
    else:
      return "%s_%s_%s_%s"%(row['gene'],row['refseq'],row['mutation'],row['flavor'])

  df_all_jobs['gene_refseq_mutation_flavor'] = df_all_jobs.apply(build_gene_refseq_mutation_flavor, axis=1)
  # import pdb; pdb.set_trace()
  
  def build_slurm_file(row):
    return "%s/%s.slurm"%(row['slurm_directory'],row['gene_refseq_mutation_flavor'])
  df_all_jobs['slurm_file'] = df_all_jobs.apply(build_slurm_file,axis=1)

  launch_strings = {}
  launch_strings['PathProx'] = []
  launch_strings['ddG'] = []
  launch_strings['SequenceAnnotation'] = []
  launch_strings['MakeGeneDictionaries'] = []

  statusdir_list = []
 
  # Build the slurm file if we are launching anew
  # import pdb; pdb.set_trace()
  allrows_cwd = None
  for index,row in df_all_jobs.iterrows():
    if (not df_prior_success.empty) and (row['uniquekey'] in df_prior_success.index):
      logger.info("Job %s was successful previously.  Slurm file will not be recreated"%row['uniquekey'])
    else:
      assert( allrows_cwd == None or allrows_cwd == row['cwd'])
      allrows_cwd = row['cwd']
      # Create the meat of the slurm script for this job
      launch_string = "%(command)s %(options)s --outdir %(outdir)s/%(flavor)s --uniquekey %(uniquekey)s"%row
 
      uniquekey_launchstring = (row['uniquekey'],launch_string) 
      if "PathProx" in row['flavor']: # Launch PathProx Clinvar and COSMIC both with same parameters
        launch_strings['PathProx'].append(uniquekey_launchstring)
      elif row['flavor'] in launch_strings: # Great for ddG, SequenceAnnocation, MakeGeneDicts
        launch_strings[row['flavor']].append(uniquekey_launchstring)
      else:
        print "I don't know this job flavor: ",row['flavor']
        sys.exit(1)

      # Before launching the slurm lob, make path to its 
      # status directory and erase all files from that
      # directory if any there from previous run
      statusdir = os.path.join(row['outdir'],row['flavor'],"status")
      if not os.path.exists(statusdir):
        try:
          os.makedirs(statusdir)
        except:
          pass
      set_capra_group_sticky(statusdir)
    
      for the_file in os.listdir(statusdir):
        file_path = os.path.join(statusdir, the_file)
        logger.info('Deleting old file %s'%file_path)
        try:
          if os.path.isfile(file_path):
            os.unlink(file_path)
        except Exception as e:
          logger.exception("Unable to delete file in status directory")
          sys.exit(1)
      # Well, we have not _actually_ submitted it _quite_ yet - but just do it
      statusdir_list.append(statusdir)
    

  # Dictionares to match each launched job unique key to it';s job id and array id
  jobids = {}
  arrayids = {}
  jobinfo = {}
  jobprogress = {}
  exitcodes = {}

  for index,row in df_all_jobs.iterrows():
    if (not df_prior_success.empty) and (row['uniquekey'] in df_prior_success.index):
      logger.info("Job %s was successful previously.  Prior results will be copied to new workstatus file"%row['uniquekey'])
      prior_success = df_prior_success.loc[row['uniquekey']]
      jobinfo[row['uniquekey']] = prior_success['jobinfo']
      jobprogress[row['uniquekey']] = prior_success['jobprogress']
      jobids[row['uniquekey']] = prior_success['jobid']
      arrayids[row['uniquekey']] = prior_success['arrayid'] if ('arrayid' in prior_success) else 0
      exitcodes[row['uniquekey']] = prior_success['ExitCode']

    # import pdb; pdb.set_trace()

  for subdir in ['PathProx', 'ddG', 'SequenceAnnotation','MakeGeneDictionaries']:
    if len(launch_strings[subdir]) == 0:
      # It's noteworth if we have a normal gene entry (not casewide) and a gene-related job is not running
      if (gene == 'casewide' and subdir == 'MakeGeneDictionaries') or (gene != 'casewide' and subdir != 'MakeGeneDictionaries'):
        logger.info("No %s jobs will be run for %s %s %s"%(subdir,gene,refseq,mutation))
    else:
      slurm_directory = build_slurm_dir(subdir)
      if not os.path.exists(slurm_directory):  # python 3 has exist_ok parameter...
        os.makedirs(slurm_directory)
      set_capra_group_sticky(slurm_directory)

      slurm_stdout_directory = os.path.join(slurm_directory,"stdout")
      if not os.path.exists(slurm_stdout_directory):  # python 3 has exist_ok parameter...
        os.makedirs(slurm_stdout_directory)
      set_capra_group_sticky(slurm_stdout_directory)

      slurm_file = os.path.join(slurm_directory,"%s_%s.slurm"%(geneRefseqMutation_OR_casewideString,subdir))
      logger.info("Creating: %s"%slurm_file)
  

      with open(slurm_file,'w') as slurmf:
        slurmf.write("""\
#!/bin/sh
#
# Project        : PSB Pipeline
# Filename       : %s
# Generated By   : %s
# For work plan  : %s
# Organization   : Vanderbilt Genetics Institute,
#                : Program in Personalized Structural Biology,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University
# Generated on   : %s
# Description    : Runs a python script on a cluster to accomplish: %s
#===============================================================================
# Slurm Parameters
"""%(slurm_file,__file__,workplan,str(datetime.datetime.now()),subdir))
      
        slurmf.write("\n")
        # It is important to make shallow copies of the default dictionaries, else mutations #7 can pick up dictionary entries set by mutation #4
        if "PathProx" in subdir:
          slurmDict = dict(SlurmParametersPathProx)
        elif subdir == "ddG":
          slurmDict = dict(SlurmParametersUDNStructure)
        elif subdir == "SequenceAnnotation":
          slurmDict = dict(SlurmParametersUDNSequence)
        elif subdir == "MakeGeneDictionaries":
          slurmDict = dict(SlurmParametersMakeGeneDictionaries)
        else:
          print "I don't know what flavor(subdir) is: ",subdir
          sys.exit(1)
         
        slurmDict['output'] = "%s/%s.out"%(slurm_stdout_directory,"%s_%%A_%%a"%geneRefseqMutation_OR_casewideString)
        slurmDict['job-name'] = "%s_%s"%(geneRefseqMutation_OR_casewideString,subdir)
        jobCount = len(launch_strings[subdir])
        if jobCount > 1:
          slurmDict['array'] = "0-%d"%(jobCount-1)
        else: # We've had trobule with modifications to the default dictionaries getting passed forward.  Make sure this did not happen
          assert 'array' not in slurmDict
             
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

        slurmf.write("""
source /dors/capra_lab/users/psbadmin/psb_prep.bash
cd %s
if [ $? != 0 ]; then
echo Failure at script launch: Unable to change to directory %s
exit 1
fi
"""%(allrows_cwd,allrows_cwd))
       

        if jobCount == 1:
          # No need to fiddle with the slurm case statement
          slurmf.write(launch_strings[subdir][0][1])
        else:
          slurmf.write("case $SLURM_ARRAY_TASK_ID in\n")
          slurm_array_id = 0
          for uniquekey_launchstring in launch_strings[subdir]:
            # Write out the "case" tag 
            slurmf.write("%d)\n"%slurm_array_id)
            # Save this to put back later
            arrayids[uniquekey_launchstring[0]] = slurm_array_id
            slurm_array_id += 1
            slurmf.write(uniquekey_launchstring[1])
            # end the case tag
            slurmf.write("\n;;\n\n")
          slurmf.write("esac\n")
      logger.info("Slurm file created: %s",slurm_file)
      jobid = slurm_submit(['sbatch',slurm_file])
      # jobid = 999 # slurm_submit(['sbatch',slurm_file])
      for uniquekey_launchstring in launch_strings[subdir]:
        jobids[uniquekey_launchstring[0]] = jobid
        jobinfo[uniquekey_launchstring[0]] = 'Submitted'
        jobprogress[uniquekey_launchstring[0]] = 'In queue'
        exitcodes[uniquekey_launchstring[0]] = None

      # Now, finally, record "Submitted" in all the status info files
      for stausdir in statusdir_list:
        with open('%s/info'%statusdir,'w') as f:
          f.write('Submitted')

  # Now create a new workstatus file
  print "Recording all jobids to %s"%workstatus_filename
 
  # import pdb; pdb.set_trace() 
  df_running_jobs = df_all_jobs.merge(pd.DataFrame(list(jobids.iteritems()),columns=['uniquekey','jobid']),on='uniquekey',how='left')
  df_running_jobs = df_running_jobs.merge(pd.DataFrame(list(arrayids.iteritems()),columns=['uniquekey','arrayid']),on='uniquekey',how='left')
  df_running_jobs = df_running_jobs.merge(pd.DataFrame(list(jobinfo.iteritems()),columns=['uniquekey','jobinfo']),on='uniquekey',how='left')
  df_running_jobs = df_running_jobs.merge(pd.DataFrame(list(jobprogress.iteritems()),columns=['uniquekey','jobprogress']),on='uniquekey',how='left')
  df_running_jobs = df_running_jobs.merge(pd.DataFrame(list(exitcodes.iteritems()),columns=['uniquekey','exitcodes']),on='uniquekey',how='left')
  df_running_jobs.set_index(['uniquekey','jobid'],inplace=True)
  df_running_jobs.to_csv(workstatus_filename,sep='\t')
  return df_running_jobs,workstatus_filename

# Main logic here.  A period in the argument means the user wants to launch one mutation only,
# directly from a single mutation output file of psb_plan.py
if oneMutationOnly:
  launch_one_mutation(args.projectORworkplan)  #  The argument is a complete workplan filename
else:
  udn_csv_filename = os.path.join(collaboration_dir,"%s_missense.csv"%args.projectORworkplan) # The argument is an entire project UDN124356
  print "Retrieving project mutations from %s"%udn_csv_filename
  df_all_mutations = pd.DataFrame.from_csv(udn_csv_filename,sep=',')
  print "Launching all jobs for %d mutations"%len(df_all_mutations)
  for index,row in df_all_mutations.iterrows():
    print "Launching %-10s %-10s %-6s"%(row['gene'],row['refseq'],row['mutation'])
    mutation_dir = os.path.join(collaboration_dir,"%s_%s_%s"%(row['gene'],row['refseq'],row['mutation']))
    if not os.path.exists(mutation_dir):  # python 3 has exist_ok parameter... 
      logger.critical("The specific mutation directory %s should have been created by psb_plan.py.  Fatal problem.")
      sys.exit(1)
       
    workplan_filename = "%s/%s_%s_%s_workplan.csv"%(mutation_dir,row['gene'],row['refseq'],row['mutation'])
    # print workplan_filename
    launch_one_mutation(workplan_filename)  #  The argument is a complete workplan filename
  casewide_workplan_filename = "casewide/casewide_workplan.csv"
  if os.path.exists(casewide_workplan_filename):
    print "Launching casewide job(s) listed in %s"%casewide_workplan_filename
    launch_one_mutation(casewide_workplan_filename)  #  The argument is a complete workplan filename
  else:
    print "No casewide file (%s) was found.  No casewide jobs will be started"%casewide_workplan_filename
  
sys.exit(0)
