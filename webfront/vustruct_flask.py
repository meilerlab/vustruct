#!/usr/bin/env python3
#
# vustruct_flask
#
# Via flask, implements RESTAPI calls which
# 
# 1) return unique job identifiers get_uuid()
# 2) launch the needed psb_*.py command line programs on the compute cluster
# 3) record and report emerging website availability
# 4) inform vustruct_webupdate of running job IDs which are needing a web refresh.
# 
#     Because vustruct_flask will launch prep tasks and cluster tasks, 
#     vustrucT_flask is executed outside the container.
#
# curl samples can be found with the various calls, just search for "curl"
# in the source code below

import logging

from flask import Flask
from flask import request
from flask import jsonify
import uuid
import threading
from datetime import datetime
import time
import json
import subprocess
import os
import sys

# requests is used to read the excel or missense filess uploaded to wordpress.
import requests
from requests_file import FileAdapter
from typing import Dict

# import time
from datetime import datetime
from logging.handlers import RotatingFileHandler
from typing import Dict, Any, Tuple, List
from typing import TextIO
# Inspect to recover bsub_submit() and slurm_submit() source code for integration
# into generated final file
import inspect


import pandas as pd
# from jinja2 import Environment, FileSystemLoader

from slurm import slurm_submit
from bsub import bsub_submit

from psb_shared import psb_config
from psb_shared import psb_perms

# UDN="/dors/capra_lab/users/mothcw/UDNtests/"
UDN = os.getenv("UDN")
if not UDN:
    errormsg = "UDN environment variable must be defined before launch of %s" % UDN
    sys.exit(errormsg)


cmdline_parser = psb_config.create_default_argument_parser(__doc__, os.path.dirname(os.path.dirname(__file__)))

args, remaining_argv = cmdline_parser.parse_known_args()

# sh = logging.StreamHandler()
LOGGER = logging.getLogger()
# LOGGER.addHandler(sh)

# Now that we've added streamHandler, basicConfig will not add another handler (important!)
LOG_FORMAT_STRING = '%(asctime)s %(levelname)-4s [%(filename)16s:%(lineno)d] %(message)s'
DATE_FORMAT_STRING = '%H:%M:%S'
logging.basicConfig(format='%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s',
datefmt='%d-%m:%H:%M:%S', )

# print("Level is %s " % os.getenv("FLASK_ENV"))
if args.debug or (os.getenv("FLASK_ENV") == 'development'):
    print("Setting Debug level")
    LOGGER.setLevel(logging.DEBUG)
elif args.verbose:
    LOGGER.setLevel(logging.INFO)
else:
    LOGGER.setLevel(logging.WARNING)
# LOGGER.setFormatter(log_formatter)

LOGGER.info("%s running web-input cases in %s", __file__, UDN)
LOGGER.info("Initializing Flask(%s)", __name__)


app = Flask(__name__)

jobs_dict = {}
job_uuids_needing_website_refresh = set()


class ActiveVUstructJob:
    def __init__(self, case_id: str, job_uuid: str) -> None:
        self.case_id = case_id
        self.job_uuid = job_uuid
        # As we launch processes, we populate this with
        # the vustruct "proprocess", "plan", "launch", "report"
        self.last_module_launched = ""
        self.last_module_returncode = 0
        self.excel_file_URI = ""

        # The running of psb_rep is particularly important.
        self.last_psb_rep_start = ""
        self.last_psb_rep_end = ""

        self.working_directory = os.path.join(UDN, "external_user_%s_%s" % (case_id, job_uuid))
        self.data_format = ""

        self.init_time = datetime.now().isoformat()

    def launch_parse_udn_report(self) -> int:
        os.makedirs(self.working_directory, exist_ok=True)
        save_cwd = os.getcwd()
        LOGGER.info("chdir(%s)", self.working_directory)
        os.chdir(self.working_directory)
        LOGGER.info("Attempting to load uploaded spreadsheet from %s", self.excel_file_URI)

        # Usually we "get" the spreadsheet from wordpress using http
        # It is convenient to be able to test all this with local files
        _session = requests.Session()
        _session.mount('file://', FileAdapter())
        response = _session.get(self.excel_file_URI)

        # Write out the UDN spreadsheet in the final directory where it will be run
        with open("external_user_%s_%s.xlsx" % (self.case_id, self.job_uuid), "wb") as spreadsheet_f:
            spreadsheet_f.write(response.content)

        launch_parse_command = "parse_udn_report.py"
        LOGGER.info("Running: %s" % launch_parse_command)

        launch_parse_return = \
            subprocess.run(launch_parse_command, shell=False,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        LOGGER.debug("%s finshed with exit code %s", launch_parse_command, launch_parse_return.returncode)
        os.chdir(save_cwd)

        self.last_module_launched = "parse"
        self.last_module_returncode = launch_parse_return.returncode

        return self.last_module_returncode

    def launch_psb_rep(self) -> int:
        save_cwd = os.getcwd()
        print("chdir(%s)" % self.working_directory)
        os.chdir(self.working_directory)
        print("Attempting to run psb_rep.py in %s" % self.working_directory)

        launch_psb_rep_command = "psb_rep.py"
        print("Running: %s" % launch_psb_rep_command)

        self.last_psb_rep_start = datetime.now().isoformat()
        launch_psb_rep_return = \
            subprocess.run(launch_psb_rep_command, shell=False,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        print("%s finished" % launch_psb_rep_command)
        self.last_psb_rep_end = datetime.now().isoformat()
        os.chdir(save_cwd)
        job_uuids_needing_website_refresh.add(self.job_uuid)

        return launch_psb_rep_return.returncode


def launch_vustruct_case_thread(vustruct_job: ActiveVUstructJob) -> subprocess.CompletedProcess:
    """
    This function is launched via the thread API.  "int he background" it will run the various components
    of the entire pipeline, and return in case of terminal failure, or completion.
    Along the way, the website files will slated for updated
    """
    global job_uuids_needing_website_refresh

    jobs_dict[vustruct_job.job_uuid] = vustruct_job

    # If we need to preprocess (parse) an input to create a missense file
    # Jump on that first
    vustruct_job.last_module_returncode = 0
    if vustruct_job.data_format == 'Vanderbilt UDN Case Spreadsheet':
        vustruct_job.last_module_launched = 'preprocess'
        vustruct_job.last_module_returncode = vustruct_job.launch_parse_udn_report()
        vustruct_job.last_module_error_message = ""
        if vustruct_job.last_module_returncode != 0:
            vustruct_job.last_module_error_message = "LATER Add error message for parse_udn_report please"
        # After we complete the parse step, we definitely need to run the psb_rep program
        # so that the OTHER flask program will install the growing website if it comes to that.
        # HOWEVER, the psb_rep program must note if the parse process has terminated with bad exit code
        # and take that into consideration.  IT may not as of yet.
        launch_psb_rep_retcd = vustruct_job.launch_psb_rep()

        if launch_psb_rep_retcd != 0:
            vustruct_job.last_report_error = "psb_rep.py for %s failed with %s" % \
                                             (vustruct_job.job_uuid, launch_psb_rep_retcd)


        if vustruct_job.last_module_returncode != 0:
            return

    launch_psb_rep_retcd = vustruct_job.launch_psb_rep()
  
    if launch_psb_rep_retcd == 0:
        job_uuids_needing_website_refresh.add(vustruct_job.job_uuid)
    else:
        print("UUGH - psb_rep for %s failed with %s" % (vustruct_job.job_uuid, launch_psb_rep_retcd))

    return launch_psb_rep_retcd

# def launch_vustruct_case(uuid_str: str):
#     psb_plan_command = \
#         ("export UDN=/dors/capra_lab/users/mothcw/UDNtests/; " +
#          "mkdir -p cd $UDN/fakecase_%s; " % uuid_str) + \
#         "singularity exec --bind /dors/capra_lab/ $UDN/development.simg psb_plan.py"
#     print("UUID = %s" % uuid_str )
#     print("running: %s" % psb_plan_command)
#     psb_plan_run_return = \
#         subprocess.run(psb_plan_command, shell=True,
#                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
 

@app.route('/get_uuid', methods=['POST'])
def get_uuid():
    """
    This first REST API call from the Wordpress functions.php
    SIMPLY verifies that this FLASK module is listening.

    We do not record the uuid here in FLASK.  We return a new UUID
    and the Wordpress code will _later_ call launch_vustruct() when
    all checks out on the UI side and SUBMIT is pressed there.

    You can test this RESTAPI call easily with:
    curl -X POST http://localhost:5000/get_uuid
    """

    # LOGGER.info ("get_uuid(): The request I got is: %s" % (request.json))
    # Launch psb_plan.py

    uuid_str = str(uuid.uuid4())

    LOGGER.info("POST to /getuuid returning %s", uuid_str)

    return jsonify({"job_uuid": uuid_str})  # 'Hello ' #  + rq.get('name', 'No name'))

@app.route('/launch_vustruct', methods = ['POST', 'GET'])
def launch_vustruct():
    """
    There are a variety of ways to launch a new VUstruct case - 
    Some start frmo the missense.csv.  Others from a spreadsheet format.
    
    An example curl command would be:
    curl -X POST http://localhost:5000/launch_vustruct \ 
        -H 'Content-Type: application/json' \
        -d'{"data_format": "Vanderbilt UDN Case Spreadsheet",
            "job_uuid": "12435",
            "case_id": "test",
            "excel_file_URI": "file:////dors/capra_lab/users/mothcw/UDNtests/fakecase0/fakecase0.xlsx"}'
    """


    LOGGER.debug ("/launch_vustruct: request.json=%s", str(request.json))

    for required_key in ['case_id', 'job_uuid', 'data_format']:
        if not required_key in request.json:
            return "Request JSON is missing required key: " + required_key

    vustruct_job = ActiveVUstructJob(
        case_id=request.json['case_id'],
        job_uuid=request.json['job_uuid'])
    vustruct_job.data_format = request.json['data_format']

    if 'excel_file_URI' in request.json:
        vustruct_job.excel_file_URI = request.json['excel_file_URI']

    # This so far has been very quick.  Now run the background thread to run the entire pipeline
    # Return to caller immediately of course.

    psb_plan_thread = threading.Thread(target=launch_vustruct_case_thread, args=(vustruct_job,))
    psb_plan_thread.start()

    # Send back the entire json that was send to us AnD add a launch_status to the gemisch
    # json_to_return_after_launch = request.json | {"launch_status": "VUstruct now running on accre cluster"}
    vustruct_job_json = jsonify(vars(vustruct_job))
    LOGGER.debug("Sending back this json %s", vustruct_job_json )

    return vustruct_job_json

@app.route('/jobs_needing_refresh', methods = ['POST', 'GET'])
def jobs_needing_refresh():
    global job_uuids_needing_website_refresh
    # Reset the listof jobs needing website refresh so we don't keep going and going
    jobs_needing_website_refresh = [vars(jobs_dict[uuid]) for uuid in job_uuids_needing_website_refresh]
    # Clear out the set of jobs needing website refresh
    job_uuids_needing_website_refresh = set()
    return jsonify(jobs_needing_website_refresh)

@app.route('/peek_jobs_needing_refresh', methods = ['POST', 'GET'])
def peek_jobs_needing_refresh():
    global job_uuids_needing_website_refresh
    jobs_needing_website_refresh = [vars(jobs_dict[uuid]) for uuid in job_uuids_needing_website_refresh]
    # With this peek_ call for testing, I do not clear out the set
    LOGGER.info("peek_jobs_needing_refresh count=%d", len(jobs_needing_website_refresh))
    return jsonify(jobs_needing_website_refresh)


@app.route('/get_all_jobs', methods = ['POST'])
def get_all_jobs():
    return jsonify([vars(active_vustruct) for active_vustruct in jobs_dict.values()])

@app.route('/add_job', methods = ['POST'])
def add_job():
    new_job_info = request.json
    LOGGER.info("add_job request json = %s" % new_job_info)
    active_vustruct = ActiveVUstructJob(new_job_info['case_id'], new_job_info['job_uuid'])
    jobs_dict[new_job_info['job_uuid']] = active_vustruct


    return jsonify(vars(active_vustruct))
