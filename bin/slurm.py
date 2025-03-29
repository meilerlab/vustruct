#!/usr/bin/env python3
#
# Project        : PSB Pipeline
# Filename       : slurm.py
# Authors        : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : mike.sivley@vanderbilt.edu  Chris.Moth@vanderbilt.edu
# Date           : 2017-01-22                  2018 and following
# Description    : Launches and tracks SLURM submissions.
# =============================================================================#
"""
Several slurm_* functions to ease interactions between the VUStruct pipeline
and a SLURM cluster.

This code is shared by the psb_launch.py and psb_monitor.py
respective job launchers and job monitors

It is also used by the ddG mass-launch routines

"""

import subprocess
import sys
import os
import json
from time import sleep
# from copy import deepcopy

import logging
import numpy as np
import pandas as pd 

LOGGER=logging.getLogger(__name__)

# logging.basicConfig(format='%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s',
#                    datefmt='%d-%m-%Y:%H:%M:%S', )


def slurm_jobstate_isfinished(JobState):
    return JobState in ["TIMEOUT", "SPECIAL_EXIT", "PREEMPTED", "NODE_FAIL",
                        "FAILED", "COMPLETED", "CANCELLED", "BOOT_FAIL"]


# Given a job id, return the output of the slurm scontrol command as a dictionary
def slurm_scontrol_show_job(jobid):
    tries = 0
    info = {}
    while (tries < 10):  # ACCRE can fail to respond - but we'll give up after 30 minutes
        scontrol_show_job = ["scontrol", "show", "job", jobid]
        command_line = ' '.join(scontrol_show_job)
        LOGGER.info("Executing: %s", command_line)
        try:
            p = subprocess.Popen(scontrol_show_job, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except OSError as e:
            msg = "Failed to run '%s'\n%s\n--->>> You do not seem to be logged in to a slurm cluster.\n" % (
            scontrol_show_job, str(e))
            LOGGER.critical(msg)
            sys.exit(1)
        stdout, stderr = p.communicate()
        if not stderr:
            if len(stdout.strip().split()) < 2:
                LOGGER.warning(
                    'Retrying because of no stdout from slurm command %s' % command_line)
            else:
                info = dict(tuple(info.split('=')) for info in stdout.split() if len(info.split('=')) == 2)
                if ("JobState" in info) and ("JobName" in info):  # Great this is super!
                    return info
        elif b"Socket timed out" in stderr:
            # SLURM is lagging. Wait 30s then resubmit
            LOGGER.warning(
                'Retrying because of cluster "Socket timed out" no stdout from slurm command %s' % command_line)
            sleep(30)
        else:
            # Our job could easily have finished and slurm does not know about it anymore
            if (b"Invalid job id specified" in stderr):
                info["JobName"] = "UNKNOWN"
                info["JobId"] = str(jobid)
                info["JobState"] = "UNKNOWN"
                info["RunTime"] = "UNKNOWN"
                info["NodeList"] = "UNKNOWN"
                return info

            LOGGER.warning('Retrying because "%s" gave stderr: "%s"', command_line, stderr)
            sleep(10)
            # Unknown slurm submission error.

        # Otherwise, it is more likely that we had a socket problem and should retry
        tries += 1
        LOGGER.warning(
            'Try %d: "scontrol show job %s" failed.  Retrying after 10 sec sleep.' % (tries, str(jobid)))
        sleep(10)
    else:  # yuck - ran out of tries
        info["JobId"] = str(jobid)
        info["JobState"] = "UNKNOWN"
    return info


def slurm_submit(job_submit_command_line):
    sbmt_fail = True
    run_result = None
    while sbmt_fail:
        LOGGER.info(job_submit_command_line)

        try:
            # Do NOT use any of the newer parameters after python 3.6
            # because this function is copied by psb_launch.py into a script must work
            # outside of anaconda, in a generic older python environment.
            run_result = subprocess.run(
                args=job_submit_command_line,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                timeout=120, check=True) # If there is an exception, then e.stderr has stderr

        except AttributeError as e:
            msg = "Exception: %s\n\
This module requires the subprocess.run() function, which is only available in Python 3.5 and higher" % str(e)
            LOGGER.critical(msg)
            sys.exit(1)

        except subprocess.TimeoutExpired:
            LOGGER.warning('sbatch timeout on %s', job_submit_command_line)
            sleep(30)
            continue
        except (subprocess.CalledProcessError, OSError) as e:
            msg = "Failed to run '%s'\nException: %s\nstderr: %s--->>> You do not seem to be logged in to a slurm cluster.\n" % (
            job_submit_command_line, str(e), e.stderr)
            LOGGER.critical(msg)
            sys.exit(1)


        # Extract the job ID
        if run_result.returncode == 0 and run_result.stdout.startswith(b"Submitted batch job "):  # Case of clear success
            LOGGER.info('sbatch successfully launched: %s', run_result.stdout.decode('latin'))
            sbmt_fail = False
            break
        else:  # Non-zero return code
            LOGGER.warning(
                'Dubious return of proc.communicate() for job "%s" with returncode %s\nstdout: %s\nstderr: %s' % (
                job_submit_command_line, run_result.returncode, run_result.stdout.decode('latin'), run_result.stderr.decode('latin')))
        if not run_result.stderr:
            if len(run_result.stdout.strip().split()) < 2:
                LOGGER.warning('Retrying because of no stdout from slurm command %s', job_submit_command_line)
            else:
                # No stderr -> job submitted successfully
                sbmt_fail = False
        elif b"Socket timed out" in run_result.stderr:
            # SLURM is lagging. Wait 30s then resubmit
            LOGGER.warning(
                'Retrying because of cluster "Socket timed out" no stdout from slurm command %s' % job_submit_command_line)
            sleep(30)
        else:
            LOGGER.warning('Retrying because "%s" gave stderr: "%s"', job_submit_command_line, run_result.stderr.decode('latin'))
            # Unknown slurm submission error.
            # Don't give up - loop again !!!
            # Was raise Exception("%s\n"%stderr)

    cluster_jobno = None
    if not sbmt_fail:
        # Extract and return the jobid (e.g. "Submitted batch job 12345678")
        cluster_jobno = run_result.stdout.strip().split()[-1]
        # Convert the bytes to str
        cluster_jobno = cluster_jobno.decode('latin')
    return cluster_jobno


class SlurmJob:
    def __init__(self, script, args=None, array=None, account=None, debug=False):
        if args is None:
            args = []
        self.script = script  # File name or content of slurm script
        if not os.path.exists(script):
            self.script = self._create_temp(self.script)
        self.args = args  # Arguments to the slurm script
        self.array = array  # Array string (e.g. "0-99")
        self.account = account  # SLURM account name (if not user default)
        self.debug = debug  # If true, submit job to debug partition
        # Submit the job
        self.jobid = self._submit()
        self.jobname = self.get_info()['JobName']
        self.failed = False  # Fail flag for job submission failures
        self.finished = False
        self.info = None

    def _create_temp(self, script):
        """ Writes the SLURM contents to a temporary file for submission """
        fname = "../temp/psb_pathprox_%d.slurm" % np.random.uniform(100000, 999999)
        if not os.path.exists("../temp/"):
            os.makedirs("../temp")
        with open(fname, 'wb') as fout:
            fout.write(script)
        return fname

    def _submit(self):
        """ Submits the  job to SLURM"""
        self.info = None
        job_submit = ["sbatch"]
        if self.account:
            job_submit.extend(["--account", self.account])
        if self.array:
            if self.debug:
                job_submit.extend(["--array", "0-1"])
            else:
                job_submit.extend(["--array", self.array])
        if self.debug:
            job_submit.extend(["--time", "30", "--partition", "debug"])
        job_submit.append(self.script)
        job_submit.extend(self.args)
        cluster_jobno = slurm_submit(job_submit)
        LOGGER.warning('Job number %s assigned to %s', cluster_jobno, job_submit)
        return cluster_jobno

    def get_info(self, requery_slurm=True):
        if requery_slurm:
            self.info = None
        if not self.info:
            self.info = slurm_scontrol_show_job(self.jobid)
        return self.info

    def get_state(self):
        if self.failed:
            return "FAILED"
        return self.get_info(True)["JobState"]

    def get_exit_status(self):
        return self.get_info(True)["ExitCode"]

    def get_stdout(self):
        with open(self.get_info()["StdOut"], 'rb') as fin:
            return fin.readlines()

    def get_stderr(self):
        with open(self.get_info()["StdErr"], 'rb') as fin:
            return fin.readlines()

    def get_submit_time(self):
        return self.get_info()["SubmitTime"]


    def get_start_time(self):
        return self.get_info()["StartTime"]

    def get_end_time(self):
        return self.get_info()["EndTime"]

    def is_pending(self):
        return self.get_state() == "PENDING"

    def is_running(self):
        return self.get_state() == "RUNNING"

    def is_complete(self):
        return self.get_state() == "COMPLETED"

    def is_failed(self):
        return self.get_state() == "FAILED"

    def is_cancelled(self):
        return self.get_state() == "CANCELLED"

    def is_timeout(self):
        return self.get_state() == "TIMEOUT"

    def is_finished(self):
        if not self.finished:
            jobstate = self.get_state()
            self.finished = slurm_jobstate_isfinished(jobstate) or jobstate == "UNKNOWN"
            if self.finished:
                # Record job exit info and stop polling SLURM
                self.info = self.get_info()
        return self.finished


def slurm_squeue(user: str, timeout_seconds=120) -> tuple[int, str, str]:
    """
    Run the squeue command for the given cluster user account.
    Returns  the process return code (non zero is error) and stdout and stderr
    In the event of a timeout, then return 1,'',error message

    squeue cmdline options are explained at https://slurm.schedmd.com/squeue.html

    """

    # Note again that --array currently "does nothing" when --json is used
    # See the flatten_to_df code where I handle this problem
    squeue_cmd = "squeue --array --json --user %s" % user
    LOGGER.info("Running: %s" % squeue_cmd)

    try:
        parse_return = \
            subprocess.run(squeue_cmd, shell=True, encoding='UTF-8', timeout=timeout_seconds,
                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if parse_return.returncode == 0:
            LOGGER.debug(f"{squeue_cmd} finshed successfully")
        else:
            LOGGER.error(f"{squeue_cmd} finshed with exit code {parse_return.returncode}: {parse_return.stderr}")

        return (parse_return.returncode, parse_return.stdout, parse_return.stderr)
    except subprocess.TimeoutExpired as ex:
        LOGGER.error(f"{squeue_cmd} timeout with exception {std(ex)}")
        return (
            1,  # Return 1 (nonzero) to show the process failed
            "", # No stdout in this case
            str(ex) # Return the timeout exception message
            )

    unexpected_error_string = f"{squeue_cmd} failure.  Should not arrive here in slurm_squeue()"
    return(1,unexpected_error_string,unexpected_error_string)


def flatten_squeue_stdout_to_df(squeue_stdout_json: str,expand_pending_array_tasks=True) -> pd.DataFrame:
    """Process jaon returned from squeue --json so that we have a flat dataframe of information
       indexed/indexable on job_key

       "job_key" is formatted to contain job_id always, concatenated with underscore and 
        array id IF there is an array Id (if the job part of a slurm array

       The idea is that the caller, armed with a "name" can match up a running application
       to these rows of data.

       A weird aspect of the squeue --json command is that, for pending ARRAY jobs, it returns only
       a single job entry for all of them, and relies on array_task_string to enumerate the still-pending
       array members.  Setting "expand_pending_array_tasks" causes _this_ code to work around that slurm
       but which I reported here: https://support.schedmd.com/show_bug.cgi?id=21686
       """

    squeue_df_schema = {
        'job_key': str,   # The catenation of array_job_id (orjob_id) UNDERSCORE array_task_id if available
        'job_id': int,    # This is the _specific_ JOB ID to cancel - and must NOT be used in job_key
        'name': str,
        'array_job_id': pd.Int64Dtype(),  # This is the BASE job NNNNNNN# in case of an array launch
        'array_task_id': pd.Int64Dtype(), # This is returned as np NAN for non-arrays else a float 0.0, 1.0 etc for array jobs
        'array_task_string': str, # Returned in format 0-N in cases of an array job, and is used by this code for "PENDING" jobs
        'job_state': str,
        'start_time': float,
        'end_time': float,
        'time_limit': float}

    df = pd.DataFrame(columns=squeue_df_schema.keys()).astype(squeue_df_schema)

    try:
        squeue_dict = json.loads(squeue_stdout_json)
    except json.decoder.JSONDecodeError as ex:
        LOGGER.critical("Unable to decode squeue --json stdout as json: %s" % squeue_stdout_json[:40])
        return df
    except TypeError as ex:
        LOGGER.critical("Unable to decode squeue --json stdout as json: %s" % squeue_stdout_json[:40])
        return df

    # if 'jobs' not in squeue_dict or len(squeue_dict['jobs']) == 0:
    # It is entirely reasonable for squeue to return jobs as [], and we fall through for consistency in that case running
    if 'jobs' not in squeue_dict:
        LOGGER.critical("No 'jobs' key seen in squeue --json return: %s" % squeue_stdout_json[:40])
        return df

    job_keys = []
    job_ids = []
    names = []
    array_job_ids = []
    array_task_ids = []
    array_task_strings = []
    job_states = []
    start_times = []
    end_times = []
    time_limits = []

    for job_dict in squeue_dict['jobs']:

        # array_job_id and array_task_id are found in these triple-dictionaries
        # and you have to make sure the value is 'set' and is of expected form
        # HOWEVER, array_job_id is ALWAYS set.  And it can be set to 0 and _that_ means
        # that we do not yet have an array_job_id

        _array_job_id = pd.NA
        if ('array_job_id' in job_dict and 
            'set' in job_dict['array_job_id'] and  job_dict['array_job_id']['set'] and
            'number' in job_dict['array_job_id'] and
             0 != int(job_dict['array_job_id']['number'])):
            _array_job_id = int(job_dict['array_job_id']['number'])

        _job_state = ''
        if 'job_state' in job_dict and len(job_dict['job_state']) > 0:
            _job_state = job_dict['job_state'][0]

        _array_task_string = job_dict.get('array_task_string','')
        # If this one is NOT set (set element is false in source json) then it was NOT
        # part of a job array.  We won't have an array_task_id UNLESS we have an _array_job_id
        _copies = 1
        if (pd.isna(_array_job_id)):
            job_keys.append(str(job_dict['job_id'])) 
            array_task_ids.append(pd.NA)
        else:
            # We have an array - but if pending and we want to break it open manually... do that
            if (expand_pending_array_tasks and 
                _job_state == "PENDING" and 
                len(_array_task_string) > 0):
                # Deal with format of the array_task_string typically 0_10 but also 4_9 or so
                _split_first_last = _array_task_string.split('-')
                if len(_split_first_last) == 2:
                    _first_pending_array_task = int(_split_first_last[0])
                    _last_pending_array_task = int(_split_first_last[1])
                else: # Not sure what else to do ?? Maybe only one pending left to launch????

                    _first_pending_array_task = int(_split_first_last[0])
                    _last_pending_array_task = int(_split_first_last[0])
                _copies = _last_pending_array_task + 1 - _first_pending_array_task
                for array_task_id in range(_first_pending_array_task,_last_pending_array_task+1):
                    array_task_ids.append(array_task_id)
                    job_keys.append(str(_array_job_id) + '_' + str(array_task_id))
            # Else we seem to have a singular array entry that is running - so just add the one
            elif ('array_task_id' in job_dict and 
                'set' in job_dict['array_task_id'] and  job_dict['array_task_id']['set'] and
                'number' in job_dict['array_task_id']):
                 array_task_ids.append(job_dict['array_task_id']['number'])
                 job_keys.append(str(job_dict['array_job_id']['number']) + '_' + str(job_dict['array_task_id']['number']))
            else:
                array_task_ids.append(pd.NA)
                # Just save the job id without underscore - as not array references in slurm
                job_keys.append(str(job_dict['job_id'])) 

        job_ids.extend(_copies * [job_dict.get('job_id',0)])
        names.extend(_copies * [job_dict.get('name','')])

        array_job_ids.extend(_copies * [_array_job_id])

        array_task_strings.extend(_copies * [job_dict['array_task_string']])

        if 'start_time' in job_dict and 'number' in job_dict['start_time']:
            start_times.extend(_copies * [job_dict['start_time']['number']])
        else:
            start_times.extend(_copies * [''])

        if 'end_time' in job_dict and 'number' in job_dict['end_time']:
            end_times.extend(_copies * [job_dict['end_time']['number']])
        else:
            end_times.extend(_copies * [''])

        # Careful - time_limit given in minutes
        if 'time_limit' in job_dict and 'number' in job_dict['time_limit']:
            time_limits.extend(_copies * [job_dict['time_limit']['number']])
        else:
            time_limits.extend(_copies * [''])

        if 'job_state' in job_dict and len(job_dict['job_state']) > 0:
            job_states.extend(_copies * [job_dict['job_state'][0]])
        else:
            job_states.extend(_copies * [''])

    # Now populate teh columns of the final dataframe with all this good gathered stuff above, for every job
    df = pd.DataFrame({
        'job_key':job_keys,
        'job_id':job_ids,
        'name': names,
        'array_job_id': array_job_ids,
        'array_task_id': array_task_ids,
        'array_task_string': array_task_strings,
        'job_state': pd.Series(job_states,dtype='str').astype(str),
        'start_time': pd.to_datetime(start_times,unit='s'),
        'end_time': pd.to_datetime(end_times,unit='s'),
        'time_limit': pd.to_timedelta(time_limits,unit='m')},
        columns=squeue_df_schema.keys()) # , dtype=squeue_schema )


    return df;
