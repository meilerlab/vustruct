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

psb_monitor.py interrogates the .../status and .../info directories of each individual
job, when the workstatus file says the job is incomplete

After gathering information, the workstatus.csv file for each mutation is updated.

A JSON file is created, which can in turn be read by the flask monitor - and 
displayed to the user as job progress information, through the web system

Configuration is supplied by the usual:
     -c global.config file
     -u user.config overrides

It could also be updated to  inquire status of slurm's "scontrol show job" feature
HOWEVER - this is very slow - and ACCRE staff have requested we not do this.
In 2024, I commented out this code.  It is largely replaced anyway by the squeue 
monitoring of the web interface flask code.
"""

# print("%s: Pipeline monitor for launched jobs.  -h for detailed help."%__file__)

import logging
import os,pwd,sys
import time, datetime
import argparse,configparser
import pprint
import json
from logging.handlers import RotatingFileHandler
from logging import handlers

from vustruct import VUStruct

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
import math
import hashlib
from jinja2 import Environment, FileSystemLoader

from psb_shared.ddg_repo import DDG_repo
from psb_shared.ddg_monomer import DDG_monomer
from psb_shared.ddg_cartesian import DDG_cartesian

# from slurm import SlurmJob,slurm_scontrol_show_job,slurm_jobstate_isfinished

cmdline_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)

from psb_shared import psb_config

cmdline_parser = psb_config.create_default_argument_parser(__doc__,os.path.dirname(os.path.dirname(__file__)))
cmdline_parser.add_argument("projectORworkstatus",type=str,
                                 help="Project ID (ex. UDN123456) or single workstatus.csv output file from psb_launch.py (or previous psb_monitor.py)  Example: ......$UDN/UDN123456/GeneName_NM_12345.1_S123A_workstatus.csv",
                                 default = os.path.basename(os.getcwd()),nargs='?')

args,remaining_argv = cmdline_parser.parse_known_args()
vustruct = VUStruct('monitor', args.projectORworkstatus, __file__)
vustruct.stamp_start_time()
vustruct.initialize_file_and_stderr_logging(args.debug)

LOGGER = logging.getLogger()

# Prior to getting going, we save the vustruct filename
# and then _if_ we die with an error, at least there is a record
# and psb_rep.py should be able to create a web page to that effect
vustruct.exit_code = 1
vustruct.write_file()



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

_latest_squeue_df = pd.DataFrame(columns=['job_key']).set_index('job_key')
# Try to get slurm squeue run
if 'cluster_type' in config_dict and config_dict['cluster_type'].lower() == 'slurm':
    import slurm
    parse_return = slurm.slurm_squeue(config_dict['cluster_user_id'])
    _latest_squeue_df = slurm.flatten_squeue_stdout_to_df(parse_return.stdout).set_index('job_key')


# collaboration_log_dir = os.path.join(collaboration_dir,"log")
# if not os.path.exists(collaboration_log_dir):  # python 3 has exist_ok parameter...
#  logging.error("%s not found.  It should have been created by psb_plan.py"%collaboration_log_dir)
#  sys.exit(1)

def monitor_one_mutation(workstatus_filename: str, variant_header_str: str) -> dict:
    # slurmInfoColumns = (['scontrolTimestamp','JobState','ExitCode','RunTime','TimeLimit',
    #                                'SubmitTime','EligibleTime','StartTime','EndTime',
    #                                'NodeList','BatchHost','NumNodes','NumCPUs','NumTasks','StdErr','StdOut','StdIn'])
    # set the type of all the slurmInfoColumns as string
    # type_dict = dict((c,str) for c in slurmInfoColumns)
    # Load the status of jobs that was created by psb_launch.py (or prior run of psb_monitor.py)
    # Load the schedule of jobs that was created by psb_launch.py (or previous run of psb_monitor.py)
    # keep_default_na causes empty strings to come is as such, and non pesky nan floats
    #
    #
    # We create an output file that the vustruct_flask can pick up and efficiently massage and transfer
    # out to dynamic web page 
    json_for_vustruct_flask = {}

    try:
        workstatus_file_df = pd.read_csv(workstatus_filename,delimiter='\t',dtype=str,keep_default_na=False).set_index('uniquekey')
    except FileNotFoundError:
        workstatus_file_df = pd.DataFrame(columns=['uniquekey']).set_index('uniquekey')

    if len(workstatus_file_df) > 0:
        # If we have total success already, then no need to keep going with analysis below!
        df_incomplete = workstatus_file_df[workstatus_file_df['ExitCode'] != '0']
        if len(df_incomplete) == 0:
            LOGGER.info("=== %s: All %d jobs completed successfully", variant_header_str, len(workstatus_file_df))
        else:
            # Otherwise, we have incomplete jobs - so there is more monitoring to be done - and that's A-OK
            LOGGER.info(">>> %s: %d jobs planned in file %s", variant_header_str, len(workstatus_file_df),os.path.basename(workstatus_filename))
    else:
        LOGGER.info("=0= %s: No jobs planned in file %s", variant_header_str, os.path.basename(workstatus_filename))
        return json_for_vustruct_flask # pd.DataFrame()
 
    # workstatus_file_df.set_index('uniquekey',inplace=True)
    # if len(workstatus_file_df) < 1:
    #     workstatus_filename_dirname = os.path.basename(os.path.dirname(os.path.normpath(workstatus_filename)))
    #     workstatus_file_basename = os.path.basename(workstatus_filename)
    #     message_text = "No rows(jobs) to monitor found in file %s."%os.path.join(
    #         workstatus_filename_dirname,workstatus_file_basename)
    #     LOGGER.info(message_text)

    df_all_jobs_original_status = workstatus_file_df.copy()


    # Life is easiest if we do not set the index until we are ready to update the workstatus_file_df at
    # the end of the monitor check-in.
    # df_job updates is populated by not only slurm's scontrol command, but also the status directory
    # that is filled in with files 'complete', 'progress', and 'info'
    # When we run 'scontrol show job' we get a _LOT_ of things back.
    # For now grab lots of these things in a new dataframe
    # slurm_scontrol_df = pd.DataFrame(columns=['jobid','uniquekey'] + slurmInfoColumns)
    #
    # slurm_scontrol_df[slurmInfoColumns] = slurm_scontrol_df[slurmInfoColumns].astype(str)
    df_job_updates = pd.DataFrame(columns=['jobid','uniquekey'])
   
    # Force all the columns to string to avoid conversions 
    df_job_updates[df_job_updates.columns] = df_job_updates[df_job_updates.columns].astype(str)
    
    for uniquekey,workstatus_row in workstatus_file_df.iterrows():
        # If we already have the ExitCode of the slurm job from prior run of this module
        # DO NOT attempt to further update.  There can be no new news after exit code arrives
        jobid = workstatus_row['jobid']

        # Not all of these items arrive from workstatus_row
        # For ddG calculations in the repo, these items are gathered
        # below and the dictionary updated
        json_for_vustruct_flask[uniquekey] = {
            'jobid': jobid,
            'arrayid': workstatus_row['arrayid'],
            'jobinfo': workstatus_row['jobinfo'],
            'jobprogress': workstatus_row['jobprogress'],
            'ExitCode': workstatus_row['ExitCode']
            } 
 
        interesting_info = {}
        interesting_info['uniquekey'] = uniquekey
        interesting_info['jobid'] = jobid

        # We don't bother to get additional information from the happily completed jobs
        if ('ExitCode' in workstatus_row) and workstatus_row['ExitCode'] and (len(workstatus_row['ExitCode']) >= 1) and (int(workstatus_row['ExitCode']) == 0):
            LOGGER.debug("%15s:%-20s Exit Code %s recorded previously"%(jobid,uniquekey,str(workstatus_row['ExitCode'])))
            continue

        # We now attempt to find out more "interesting information" about in process jobs
        # and we place this both on the stdout for comamnd line users and
        # add to the json_for_vustruct_flask for web users

        if 'Ddg' in workstatus_row['flavor'] and 'repo' in workstatus_row['outdir']:
            repo_calculation_flavor = 'ddG_cartesian' if 'artesian' in workstatus_row['flavor'] else 'ddG_monomer'

            # We need to set the log level to WARN to avoid endless info messages from the repo
           
            _save_log_level = LOGGER.level
            LOGGER.setLevel(logging.WARN)
            ddg_repo = DDG_repo(config_dict['ddg_config'],
                    calculation_flavor=repo_calculation_flavor)

            if workstatus_row['method'] == 'swiss':
                ddg_repo.set_swiss(workstatus_row['pdbid'],workstatus_row['chain'])
            elif workstatus_row['method'] == 'modbase':
                ddg_repo.set_modbase(workstatus_row['pdbid'],workstatus_row['chain'])
            elif workstatus_row['method'] == 'alphafold':
                ddg_repo.set_alphafold(workstatus_row['pdbid'],workstatus_row['chain'])
            elif workstatus_row['method'] == 'usermodel':
                ddg_repo.set_usermodel(workstatus_row['pdbid'],workstatus_row['chain'])
            else:
                ddg_repo.set_pdb(workstatus_row['pdbid'].lower(),workstatus_row['chain'])

            ddg_repo.set_variant(workstatus_row['pdbmut'])

            ddg_monomer_or_cartesian = \
                DDG_monomer(ddg_repo,workstatus_row['pdbmut']) if repo_calculation_flavor == 'ddG_monomer' else \
                DDG_cartesian(ddg_repo,workstatus_row['pdbmut'])
                
            LOGGER.setLevel(_save_log_level)

            # Ask ddg_monomer to return information about the job
            # The ddG jobs are under the watch of the repository, not us so much
            # It will return dictionary with ExitCode, jobprogress, jobinfo
            interesting_info.update(ddg_monomer_or_cartesian.jobstatus())
            LOGGER.info("Status of DDG job:\n%s"%str(interesting_info))
        else: # non DDG jobs
            # We dig down to the status directory outputs
            # The most reliable source of Exit=0 is the arrival of a completed file in the status directory
            progress_file_found = False
            if "PP_" in workstatus_row['flavor']:
                statusdir = os.path.join(workstatus_row['outdir'],workstatus_row['flavor'],"status")
            else: # For MusiteDeep and ScanNet we don't have an additional flavor directory
                statusdir = os.path.join(workstatus_row['outdir'],"status")
            try:
                with open("%s/progress"%statusdir) as progressfile:
                    status_text = progressfile.read().replace('\n','')
                    if status_text and len(status_text):
                        interesting_info['jobprogress'] = status_text
                        progress_file_found = True
            except FileNotFoundError as ex:
                pass
 
            info_file_found = False
            try:
                with open("%s/info"%statusdir) as infofile:
                    status_text = infofile.read().replace('\n','')
                    if status_text and len(status_text):
                        interesting_info['jobinfo'] = status_text
                        info_file_found = True
            except FileNotFoundError as ex:
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
                    LOGGER.info("    %s",interesting_info['jobprogress'])
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
                    slurm_scontrol_df = slurm_scontrol_df.append(info_as_series,ignore_index=True)
                """ 
        ################################################################################
        # Whether interesting_info populated by the ddg* routines
        # Or by our monitoring of "our" jobs, update the dataframe with all news on
        # progress of the jobs 
        # deprecated -> slurm_scontrol_df = slurm_scontrol_df.append(info_as_series,ignore_index=True)

        info_as_series = pd.Series(data=interesting_info)
        df_job_updates = pd.concat([df_job_updates,pd.DataFrame(info_as_series).transpose()],ignore_index=True)
        json_for_vustruct_flask[uniquekey].update(interesting_info)


    previous_workstatus_filename = workstatus_filename + ".previous"
    # print "Renaming current status file to %s"%previous_workstatus_filename
    try:
        os.remove(previous_workstatus_filename)
    except Exception as ex:
        pass
     
    
    
    # If the slurmInfoColumns are not alreadt in workstatus_file_df, then add those columns 
    # will add them with NaN initializers.  If these columns are already there, do nothing
    # (I've tried concat (only works pre 0.20 pandas) and other things - and have 
    #  given up on elegant solutions)
    # for column in slurmInfoColumns:
    #     if column not in workstatus_file_df:
    #         workstatus_file_df[column] = None
    
    # Old idea failed on different versions: workstatus_file_df = pd.concat([workstatus_file_df,pd.DataFrame(columns=slurmInfoColumns)])
    # Can't work because done at top: workstatus_file_df.set_index(['uniquekey','jobid'],inplace=True)
 
    # This one here should work because we've not set the index before on it
    df_job_updates = df_job_updates.set_index('uniquekey')
    
    # Take all the non-NA values from slurm_scontrol_df, and modify workstatus_file_df to include them
    df_all_jobs_original_status = workstatus_file_df.copy()


    workstatus_file_df.update(df_job_updates)
    df_incomplete = workstatus_file_df[workstatus_file_df['ExitCode'] != '0']

    # Are there any changes to what we know about the status of the jobs?
    if df_all_jobs_original_status.to_string(header=False) == workstatus_file_df.to_string(header=False):
        if len(df_incomplete) == 0:
            # This happens out to right of earlier printing...
            # LOGGER.info("All %d jobs completed successfully %s", len(workstatus_file_df), os.path.basename(workstatus_filename))
            return json_for_vustruct_flask
        if "casewide" in workstatus_filename:
            LOGGER.info("    No updates to status of casewide jobs")
        else:
            LOGGER.info("    No updates to status of jobs for this mutation: %s", os.path.basename(workstatus_filename))
    else:
        os.rename(workstatus_filename,previous_workstatus_filename)
    
        LOGGER.info("    Recording workstatus updates to %s", os.path.basename(workstatus_filename))
        try:
            workstatus_file_df.to_csv(workstatus_filename,sep='\t',index=True)
        except Exception as ex:
            # If we cannot save the new csv file, then we must (attempt to!) restore the old one!
            os.remove(workstatus_filename)
            os.rename(workstatus_filename,previous_workstatus_filename)
            LOGGER.exception("!!!! Serious failure saving new updated file %s.\n%s\n Prior file restored", workstatus_filename,str(ex))
            sys.exit(1)
    
    if len(df_incomplete) == 0: # Should never happen, because of short-circuit above
        LOGGER.info("    All %d jobs completed successfully", len(workstatus_file_df))
        pass
    else:
        incomplete_job_info = []
        incomplete_job_info.append("    %2d of %2d jobs still incomplete for %s:"%(len(df_incomplete),len(workstatus_file_df), variant_header_str))
        incomplete_job_info.append("%15s:%-7.7s %-25.25s  %s"%('Jobid','qState','Flavor',"Info"))
        for index,row in df_incomplete.iterrows():
            infostring = ""
            if 'jobinfo' in row and len(row['jobinfo']) > 1:
                infostring = row['jobinfo']
            elif 'JobState' in row and len(row['JobState']) > 1:
                infostring = row['JobState']

            # row['arrayid'] should either be a 0-length string (meaning not an array job)
            # or a string in floating point format
            array_id = row['arrayid']
            if len(array_id) == 0:
                array_id = None
                # The job_key will only be the master non-array job id for this single job per slurm file
                job_key = str(row['jobid'])
            else:
                array_id = math.trunc(float(row['arrayid']))
                job_key = "%s_%d" % (str(row['jobid']), array_id)

            try:
                squeue_row = _latest_squeue_df.loc[job_key]
                squeue_state = squeue_row['job_state']
            except KeyError:
                squeue_state = 'Unknown'

            # It _could_ be that the job is PENDING, in which case
            # Squeue returns a false - so lets do a search if there was an array at launch
            if squeue_state == 'Unknown' and array_id is not None:
                try:
                    squeue_row = _latest_squeue_df.loc[job_key.split('_')[0]]
                    squeue_state = squeue_row['job_state']
                except KeyError:
                    squeue_state = 'Unknown'
                    pass
   
            incomplete_job_info.append("%15s:%-7.7s %-25s  %-20s"%(job_key,squeue_state,row['pdbid'] + '_' + row['flavor'],infostring))
        incomplete_job_info.append("")
        LOGGER.info("%s" % "\n".join(incomplete_job_info))

    return json_for_vustruct_flask
    
# Main logic here.  A period in the argument means the user wants to launch one mutation only,
# directly from a single mutation output file of psb_plan.py
all_mutations_json_for_vustruct_flask = {}
if oneMutationOnly:
    monitor_one_mutation(args.projectORworkstatus,"")  #  The argument is a complete workstatus filename
else:
    udn_csv_filename = os.path.join(collaboration_dir,"%s_vustruct.csv"%args.projectORworkstatus) # The argument is an entire project UDN124356
    print("Retrieving project mutations from %s"%udn_csv_filename)
    df_case_vustruct_all_variants = pd.read_csv(udn_csv_filename,sep=',',index_col = None,keep_default_na=False,encoding='utf8',comment='#',skipinitialspace=True)
    print("Monitoring all jobs for %d mutations"%len(df_case_vustruct_all_variants))
    df_case_vustruct_all_variants.fillna('NA',inplace=True)
    # If this is new format, trim out non-missense variant rows from the vustruct file
    if 'effect' in df_case_vustruct_all_variants:
        df_vustruct_missense_variants = df_case_vustruct_all_variants[ df_case_vustruct_all_variants['effect'].str.contains('missense') ]
    else: # Old format is all missense.  Simple retain all of it
        df_vustruct_missense_variants = df_case_vustruct_all_variants

    for index,row in df_vustruct_missense_variants.iterrows():
        if 'RefSeqNotFound_UsingGeneOnly' in row['refseq']:
            row['refseq'] = 'NA'
        gene_refseq_mutation = "%s_%s_%s" % (row['gene'],row['refseq'],row['mutation'])
        # print without a newline - monitor_one_mutation will add one
        variant_header_str = "%d of %d: %-10.10s %-14.14s %-6.6s" % (index+1,len(df_case_vustruct_all_variants),row['gene'],row['refseq'],row['mutation'])
        mutation_dir = os.path.join(collaboration_dir,"%s"% gene_refseq_mutation)
        if not os.path.exists(mutation_dir):  # python 3 has exist_ok parameter... 
            logging.critical("The specific mutation directory %s should have been created by psb_launch.py.  Fatal problem."%mutation_dir)
            sys.exit(1)

        workstatus_filename = "%s/%s_workstatus.csv"%(mutation_dir,gene_refseq_mutation)
        # print workstatus_filename
        all_mutations_json_for_vustruct_flask[gene_refseq_mutation] = monitor_one_mutation(workstatus_filename, variant_header_str)  #  The argument is a complete workstatus filename
        # print("-" * 80)
    # Finally monitor the job(s) under the casewide banner (Gene Dictionary creation)
    workstatus_filename = "casewide/casewide_workstatus.csv"
    if os.path.exists(workstatus_filename):
        # print workstatus_filename
        variant_header_str = "Casewide jobs"
        all_mutations_json_for_vustruct_flask['casewide'] = monitor_one_mutation(workstatus_filename,variant_header_str)  #  The argument is a complete workstatus filename
    else:
        LOGGER.warning("No casewide work (no %s) to monitor", workstatus_filename)
    # print("-" * 80)

    json_filename = os.path.join(collaboration_dir,"%s_psb_monitor.json"%args.projectORworkstatus) # The argument is an entire project UDN124356
    with open(json_filename,'w') as json_f:
        json.dump(all_mutations_json_for_vustruct_flask, json_f, indent=4)

# To get here, we have a happy ending - update the .json record
vustruct.exit_code = 0
vustruct.stamp_end_time()
vustruct.write_file()
