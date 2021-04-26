#!/usr/bin/env python3
#
# Project        : PSB Pipeline
# Filename       : slurm.py
# Authors        : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2017-01-22
# Description    : Launches and tracks SLURM submissions.
# =============================================================================#

import subprocess as sp
import sys
import os
from time import sleep
# from copy import deepcopy

import logging
# from logging.handlers import RotatingFileHandler
# from logging import handlers
import numpy as np

logging.basicConfig(format='%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s',
                    datefmt='%d-%m-%Y:%H:%M:%S', )


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
        logging.getLogger(__name__).info("Executing: %s", command_line)
        try:
            p = sp.Popen(scontrol_show_job, stdout=sp.PIPE, stderr=sp.PIPE)
        except OSError as e:
            msg = "Failed to run '%s'\n%s\n--->>> You do not seem to be logged in to a slurm cluster.\n" % (
            scontrol_show_job, str(e))
            logging.getLogger(__name__).critical(msg)
            sys.exit(1)
        stdout, stderr = p.communicate()
        if not stderr:
            if len(stdout.strip().split()) < 2:
                logging.getLogger(__name__).warning(
                    'Retrying because of no stdout from slurm command %s' % command_line)
            else:
                info = dict(tuple(info.split('=')) for info in stdout.split() if len(info.split('=')) == 2)
                if ("JobState" in info) and ("JobName" in info):  # Great this is super!
                    return info
        elif b"Socket timed out" in stderr:
            # SLURM is lagging. Wait 30s then resubmit
            logging.getLogger(__name__).warning(
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

            logging.getLogger(__name__).warning('Retrying because "%s" gave stderr: "%s"', command_line, stderr)
            sleep(10)
            # Unknown slurm submission error.

        # Otherwise, it is more likely that we had a socket problem and should retry
        tries += 1
        logging.getLogger(__name__).warning(
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
        logging.getLogger(__name__).info(job_submit_command_line)

        try:
            # Do NOT use any of the newer parameters after python 3.6
            # because this function is copied by psb_launch.py into a script must work
            # outside of anaconda, in a generic older python environment.
            run_result = sp.run(
                args=job_submit_command_line,
                stdout=sp.PIPE, stderr=sp.PIPE,
                timeout=120, check=True)

        except AttributeError as e:
            msg = "Exception: %s\n\
This module requires the subprocess.run() function, which is only available in Python 3.5 and higher" % str(e)
            logging.getLogger(__name__).critical(msg)
            sys.exit(1)

        except sp.TimeoutExpired:
            logging.getLogger(__name__).warning('sbatch timeout on %s', job_submit_command_line)
            sleep(30)
            continue
        except (sp.CalledProcessError, OSError) as e:
            msg = "Failed to run '%s'\nException: %s\n--->>> You do not seem to be logged in to a slurm cluster.\n" % (
            job_submit_command_line, str(e))
            logging.getLogger(__name__).critical(msg)
            sys.exit(1)


        # Extract the job ID
        if run_result.returncode == 0 and run_result.stdout.startswith(b"Submitted batch job "):  # Case of clear success
            logging.getLogger(__name__).info('sbatch successfully launched: %s', run_result.stdout)
            sbmt_fail = False
            break
        else:  # Non-zero return code
            logging.getLogger(__name__).warning(
                'Dubious return of proc.communicate() for job "%s" with returncode %s\nstdout: %s\nstderr: %s' % (
                job_submit_command_line, run_result.returncode, run_result.stdout, run_result.stderr))
        if not run_result.stderr:
            if len(run_result.stdout.strip().split()) < 2:
                logging.getLogger(__name__).warning('Retrying because of no stdout from slurm command %s', job_submit_command_line)
            else:
                # No stderr -> job submitted successfully
                sbmt_fail = False
        elif b"Socket timed out" in run_result.stderr:
            # SLURM is lagging. Wait 30s then resubmit
            logging.getLogger(__name__).warning(
                'Retrying because of cluster "Socket timed out" no stdout from slurm command %s' % job_submit_command_line)
            sleep(30)
        else:
            logging.getLogger(__name__).warning('Retrying because "%s" gave stderr: "%s"', job_submit_command_line, run_result.stderr)
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
        logging.getLogger(__name__).warning('Job number %s assigned to %s', cluster_jobno, job_submit)
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
