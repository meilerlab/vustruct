#!/usr/bin/env python3
#
# vustruct_flask
#
# Via flask, implements RESTAPI calls which
# 
# 1) return unique job identifiers and
# 2) launch the needed psb_*.py command line programs on the compute cluster
# 3) record and report emerging website availability
# 4) inform vustruct_webupdate of running job IDs which are needing a web refresh.
# 
#     Because vustruct_flask will launch prep tasks and cluster tasks, 
#     vustrucT_flask is executed outside the container.

import logging
from flask import Flask
from flask import request
from flask import jsonify
import uuid
import threading
import time
import json
import subprocess
import os
import requests
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
UDN=os.getenv("UDN")
if not UDN:
    errormsg = "UDN environment variable must be defined before launch of %s" % UDN
    sys.exit(errormsg)


cmdline_parser = psb_config.create_default_argument_parser(__doc__, os.path.dirname(os.path.dirname(__file__)))

args, remaining_argv = cmdline_parser.parse_known_args()

if args.debug:
    sh.setLevel(logging.DEBUG)
elif args.verbose:
    sh.setLevel(logging.INFO)


sh = logging.StreamHandler()
LOGGER = logging.getLogger()
LOGGER.addHandler(sh)

# Now that we've added streamHandler, basicConfig will not add another handler (important!)
LOG_FORMAT_STRING = '%(asctime)s %(levelname)-4s [%(filename)16s:%(lineno)d] %(message)s'
DATE_FORMAT_STRING = '%H:%M:%S'
log_formatter = logging.Formatter(LOG_FORMAT_STRING, DATE_FORMAT_STRING)

LOGGER.setLevel(logging.DEBUG)
sh.setLevel(logging.WARNING)
sh.setFormatter(log_formatter)

LOGGER.info("%s running web-input cases in %s" , __file__, UDN)
LOGGER.info("Initializing Flask(%s)", __name__)
app = Flask(__name__)

jobs_dict = {}
jobs_needing_website_refresh = []

class ActiveVUstructJob:
    def __init__(self, case_id: str, job_uuid: str) -> None:
        self.case_id = case_id
        self.job_uuid = job_uuid
        self.working_directory = os.path.join(UDN,"external_user_%s_%s" % (case_id, job_uuid))
  
    def launch_parse_udn_report(self, excel_file_url: str):
        os.makedirs(self.working_directory,exist_ok=True)
        save_cwd = os.getcwd()
        print("chdir(%s)" % self.working_directory)
        os.chdir(self.working_directory)
        print("Attempting to load uploaded spreadsheet from %s" % excel_file_url)
        response = requests.get(excel_file_url);
        with open("external_user_%s_%s.xlsx" % (self.case_id, self.job_uuid), "wb") as spreadsheet_f:
            spreadsheet_f.write(response.content)

        launch_parse_command = "parse_udn_report.py"
        print("Running: %s" % launch_parse_command)
        launch_parse_return = \
            subprocess.run(launch_parse_command, shell=False,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        print("%s finished" % launch_parse_command)
        os.chdir(save_cwd)

        return launch_parse_return

    def launch_psb_rep(self):
        save_cwd = os.getcwd()
        print("chdir(%s)" % self.working_directory)
        os.chdir(self.working_directory)
        print("Attempting to run psb_rep.py in %s" % self.working_directory)

        launch_psb_rep_command = "psb_rep.py"
        print("Running: %s" % launch_psb_rep_command)
        launch_psb_rep_return = \
            subprocess.run(launch_psb_rep_command, shell=False,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        print("%s finished" % launch_psb_rep_command)
        os.chdir(save_cwd)

        return launch_psb_rep_return





def launch_vustruct_case(response_json: Dict[str,str]) -> subprocess.CompletedProcess :
    global jobs_needing_website_refresh

    active_vustruct_job = ActiveVUstructJob(
        response_json['case_id'],
        response_json['job_uuid']
    )
   
    launch_parse_retcd = 0 
    if 'data_format' in response_json and response_json['data_format'] == 'Vanderbilt UDN Case Spreadsheet':
        launch_parse_retcd = active_vustruct_job.launch_parse_udn_report(response_json['excel_file_URL'])
        if launch_parse_retcd != 0:
            print("UUGH - parse_udn_Report for %s failed with %s" % (response_json['job_uuid'], launch_parse_retcd))

    launch_psb_rep_retcd = active_vustruct_job.launch_psb_rep()
  
    if launch_psb_rep_retcd == 0:
        jobs_needing_website_refresh.append(response_json['job_uuid'])
    else:
        print("UUGH - psb_rep for %s failed with %s" % (response_json['job_uuid'], launch_psb_rep_retcd))

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
 
@app.route('/get_uuid', methods = ['POST', 'GET'])
def get_uuid():
   print ("get_uuid(): The request I got is: %s" % (request.json))
   # Launch psb_plan.py

   uuid_str = str(uuid.uuid4())

   return jsonify({"job_uuid": uuid_str}) # 'Hello ' #  + rq.get('name', 'No name')

@app.route('/launch_vustruct', methods = ['POST', 'GET'])
def launch_vustruct():
   print ("launch_vustruct(): The request I got is: %s" % (request.json))
   psb_plan_thread = threading.Thread(target=launch_vustruct_case, args=(request.json,))
   psb_plan_thread.start()

   # Send back the entire json that was send to us AnD add a launch_status to the gemisch
   json_to_return_after_launch = request.json | {"launch_status": "VUstruct now running on accre cluster"}
   print("Sending back this json %s" % json_to_return_after_launch)

   return jsonify(json_to_return_after_launch)

@app.route('/jobs_needing_refresh', methods = ['POST', 'GET'])
def jobs_needing_refresh():
    global jobs_needing_website_refresh
    json_jobs_needing_website_refresh = jsonify(jobs_needing_website_refresh)
    # Reset the listof jobs needing website refresh so we don't keep going and going
    jobs_needing_website_refresh = []
    return json_jobs_needing_website_refresh
