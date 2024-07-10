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

This code 

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


def slurm_squeue(user: str) -> subprocess.CompletedProcess:
    """
    Run the squeue command for the given cluster user account.
    Returns raw stdout in json format, and one job per every array id

    squeue cmdline options are explained at https://slurm.schedmd.com/squeue.html
    """

    squeue_cmd = "squeue --array --json --user %s" % user
    LOGGER.info("Running: %s" % squeue_cmd)

    parse_return = \
        subprocess.run(squeue_cmd, shell=True, encoding='UTF-8',
                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    LOGGER.debug("%s finshed with exit code %d", squeue_cmd, parse_return.returncode)
    return parse_return

def flatten_squeue_stdout_to_df(squeue_stdout_json: str) -> pd.DataFrame:
    squeue_dict = json.loads(squeue_stdout_json)
    # I will format "job_key" to be job_id, and then if there is an array id, add underscore and the array id

    squeue_schema = {
        'job_key': str,
        'job_id': int,
        'name': str,
        'array_job_id': float,
        'array_task_id': float,
        'job_state': str,
        'start_time': float,
        'end_time': float,
        'time_limit': float}

    df = pd.DataFrame(columns=squeue_schema.keys()).astype(squeue_schema)

    if 'jobs' not in squeue_dict or len(squeue_dict['jobs']) == 0:
        return df

    job_ids = []
    names = []
    array_job_ids = []
    array_task_ids = []
    job_states = []
    start_times = []
    end_times = []
    time_limits = []

    for job_dict in squeue_dict['jobs']:
        job_ids.append(job_dict.get('job_id',0))
        names.append(job_dict.get('name',''))

        if 'array_job_id' in job_dict and 'number' in job_dict['array_job_id']:
            array_job_ids.append(job_dict['array_job_id']['number'])
        else:
            array_job_ids.append(np.nan)

        if 'array_task_id' in job_dict and 'number' in job_dict['array_task_id']:
            array_task_ids.append(job_dict['array_task_id']['number'])
        else:
            array_task_ids.append(np.nan)

        if 'start_time' in job_dict and 'number' in job_dict['start_time']:
            start_times.append(job_dict['start_time']['number'])
        else:
            start_times.append('')

        if 'end_time' in job_dict and 'number' in job_dict['end_time']:
            end_times.append(job_dict['end_time']['number'])
        else:
            end_times.append('')

        # Careful - time_limit given in minutes
        if 'time_limit' in job_dict and 'number' in job_dict['time_limit']:
            time_limits.append(job_dict['time_limit']['number'])
        else:
            time_limits.append('')

        if 'job_state' in job_dict and len(job_dict['job_state']) > 0:
            job_states.append(job_dict['job_state'][0])
        else:
            job_states.append('')

    df = pd.DataFrame({
        'job_id':job_ids,
        'name': names,
        'array_job_id': array_job_ids,
        'array_task_id': array_task_ids,
        'job_state': pd.Series(job_states).astype(str),
        'start_time': pd.to_datetime(start_times,unit='s'),
        'end_time': pd.to_datetime(end_times,unit='s'),
        'time_limit': pd.to_timedelta(time_limits,unit='m')},
        columns=squeue_schema.keys()) # , dtype=squeue_schema )

    df['job_key'] = df['job_id'].astype(str) + '_' + df['array_task_id'].astype(str)

    return df;
