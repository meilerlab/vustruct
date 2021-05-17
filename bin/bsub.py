#!/usr/bin/env python3
#
# Project        : PSB Pipeline
# Filename       : bsub.py
# Author         : Chris Moth
# Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : chris.moth@vanderbilt.edu
# Date           : 2021-04-12
# Description    : Launches LSF .bsub files and returns job number.
# =============================================================================#

import subprocess as sp
import sys
import re


import logging


def bsub_submit(launch_filename):
    """
    This function is also copied into the external launch module, by psb_launch.py
    Careful when modifying this file as it must stay compatible with old Python 3.6
    """

    submit_fail = True
    # Use 'byte' strings throughout, for Python 3.6 compatability
    job_submitted_pattern = re.compile(b'Job <(\d+)> is Submitted', re.IGNORECASE)

    job_submitted_match = None
    tries = 0
    while submit_fail and tries < 10:
        logging.getLogger(__name__).info("Calling bsub < %s", launch_filename)
        cluster_jobno = ''
        bsub_command_look = "bsub < " + launch_filename
        with open(launch_filename, 'rb') as launch_f:
            try:
                # Do NOT use newer Python 3.7 etc parameters to facilitate text checking
                run_result = sp.run(
                    args=["bsub"],
                    stdin=launch_f, stdout=sp.PIPE, stderr= sp.PIPE,
                    timeout=120, check=True)
                job_submitted_match = job_submitted_pattern.match(run_result.stdout)
                if run_result.returncode == 0 and job_submitted_match:  # Case of clear success
                    # Extract and return the jobid (e.g. "Submitted batch job 12345678")
                    logging.getLogger(__name__).info('%s successfully launched: %s',
                                                     bsub_command_look, run_result.stdout.decode('latin'))
                    cluster_jobno = job_submitted_match.group(1)
                    submit_fail = False
                else:
                    if run_result.stderr:
                        logging.getLogger(__name__).warning('Retrying because "%s" gave stderr: "%s"',
                                                            bsub_command_look, run_result.stderr.decode('latin'))
                    else:
                        logging.getLogger(__name__).warning(
                            'Dubious return of sp.run for job "%s" with returncode %s\nstdout: %s\nstderr: %s',
                            bsub_command_look, run_result.returncode, run_result.stdout.decode('latin'), run_result.stderr.decode('latin'))
                # Unknown BSUB submission error.
                # Don't give up - loop again !!!
                # Was raise Exception("%s\n"%stderr)

            except sp.TimeoutExpired:
                logging.getLogger(__name__).warning('bsub timeout on %s', launch_filename)
            except sp.CalledProcessError as error:
                msg = "Failed to run '%s'\n%s\n--->>> You do not seem to be logged in to an LSF cluster.\n" % (
                    bsub_command_look, str(error))
                logging.getLogger(__name__).critical(msg)
                sys.exit(msg)

    if submit_fail:
        sys.exit("Failed to execute bsub after %d tries" % tries)

    # Convert the bytes to str
    return cluster_jobno


if __name__ == "__main__":
    bsub_submit('test.bsub')
