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

from functools import wraps
from flask import Flask
from flask import abort
from flask import request
from flask import jsonify
import uuid
import threading
from datetime import datetime
from datetime import timedelta
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

app = Flask(__name__)
app.logger.setLevel(logging.DEBUG)

# UDN="/dors/capra_lab/users/mothcw/VUStruct/"
UDN = os.getenv("UDN")
if not UDN:
    errormsg = "UDN environment variable must be defined before launch of %s" % UDN
    sys.exit(errormsg)

config = {
    'report_interval_initial_minutes': 1.0,
    'report_interval_stretch_factor': 1.5,
    'report_interval_max_minutes': 30,
    'report_timeout_days': 4.0
}
    



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
# pw_name = pwd.getpwuid( os.getuid() ).pw_name # example jsheehaj or mothcw
os.makedirs('log', exist_ok=True)
log_filename = os.path.join('log', "vustruct_flask.log")

LOGGER.info("Additionally logging to file %s" % log_filename)
needRoll = os.path.isfile(log_filename)

local_fh = RotatingFileHandler(log_filename, backupCount=7)
formatter = logging.Formatter('%(asctime)s %(levelname)-4s [%(filename)20s:%(lineno)d] %(message)s',
                              datefmt="%H:%M:%S")
local_fh.setFormatter(formatter)
local_fh.setLevel(logging.DEBUG)
LOGGER.addHandler(local_fh)

if needRoll:
    local_fh.doRollover()

LOGGER.info("info test")

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



cases_dict = {}
case_uuids_needing_website_refresh = set()


class VUstructCaseManager:
    def __init__(self, case_id: str, case_uuid: str) -> None:
        self.case_id = case_id
        self.case_uuid = case_uuid
        # As we launch processes, we populate this with
        # the vustruct "proprocess", "plan", "launch", "report"
        self.last_module_launched = ""
        self.last_module_returncode = 0
        self.upload_file_URI = ""

        # The running of psb_rep is particularly important.
        self.last_psb_rep_start = ""
        self.last_psb_rep_end = ""

        # This prefix results in a sibling directory to our interactive runs, with 
        # a very distinctive  and long unique name however
        # This prefix is used through to the website creation.
        # However, it is rather removed from the actual text of the webpages themselves
        # during final website installation
        self.external_user_prefix = "external_user_%s_%s" % (case_id, case_uuid)
        self.working_directory = os.path.join(UDN, self.external_user_prefix)
        self.data_format = ""

        self.init_time = datetime.now().isoformat()

        # Convenient to push the current directory, do work in the case, and pop back
        self.save_cwd = None

    def fetch_input_file(self,file_extension: str) -> None:
        os.makedirs(self.working_directory, exist_ok=True)
        self.save_cwd = os.getcwd()
        LOGGER.info("chdir(%s)", self.working_directory)
        os.chdir(self.working_directory)
        LOGGER.info("Attempting to load uploaded file from %s", self.upload_file_URI)

        # Usually we "get" the spreadsheet from wordpress using http
        # It is convenient to be able to test all this with local files
        _session = requests.Session()
        _session.mount('file://', FileAdapter())
        response = _session.get(self.upload_file_URI)

        # Write out the uploaded spreadseet or vcf in the final directory where it will be run
        with open("%s.%s" % (self.external_user_prefix, file_extension), "wb") as spreadsheet_f:
            spreadsheet_f.write(response.content)

        return

    def run_coordinates_parser_command_line(self, parser_python_filename: str) -> subprocess.CompletedProcess:
        LOGGER.info("Running: %s" % parser_python_filename)

        parse_return = \
            subprocess.run(parser_python_filename, shell=True, encoding='UTF-8',
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        LOGGER.debug("%s finshed with exit code %d", parser_python_filename, parse_return.returncode)
        os.chdir(self.save_cwd)
        self.save_cwd = None
        return parse_return




    def parse_udn_report(self) -> int:
        self.fetch_input_file('xlsx')
      
        parse_return = self.run_coordinates_parser_command_line("parse_udn_report.py")
        #parse_command_line="""
        #export UDN=/dors/capra_lab/users/mothcw/VUStruct; cd %s; singularity exec ../development.simg parse_udn_report.py"""%\
        #self.working_directory

        # self.last_module_launched = "parse"
        # self.last_module_returncode = parse_return.returncode

        return parse_return.returncode

    def vcf2missense(self) -> int:
        self.fetch_input_file('vcf')
      
        parse_return = self.run_coordinates_parser_command_line("vcf2missense.py")
        #parse_command_line="""
        #export UDN=/dors/capra_lab/users/mothcw/VUStruct; cd %s; singularity exec ../development.simg vcf2missense.py"""%\
        #self.working_directory

        # self.last_module_launched = "parse"
        # self.last_module_returncode = parse_return.returncode

        return parse_return.returncode

    def vcf2missense_with_liftover(self) -> int:
        self.fetch_input_file('vcf')
      
        parse_return = self.run_coordinates_parser_command_line("vcf2missense.py --liftover")
        #parse_command_line="""
        #export UDN=/dors/capra_lab/users/mothcw/VUStruct; cd %s; singularity exec ../development.simg vcf2missense.py"""%\
        #self.working_directory

        # self.last_module_launched = "parse"
        # self.last_module_returncode = parse_return.returncode

        return parse_return.returncode

    def psb_plan(self) -> int:
        os.makedirs(self.working_directory, exist_ok=True)
        save_cwd = os.getcwd()
        LOGGER.info("chdir(%s)", self.working_directory)
        os.chdir(self.working_directory)

        # If the user uploaded a missense csv file, then we will
        # fetch that here.  BUT, if they started with another format
        # then the misesense file should already be in place..
        missense_filename = "%s_missense.csv" % self.external_user_prefix
        if self.data_format == 'MissenseCSV':
            LOGGER.info("Saving missense data from wordpress UI to local missense filename %s", missense_filename)

            with open(missense_filename,'w') as f:
                f.write(self.missense_csv);
        else:
            if not os.path.exists(missense_filename):
                LOGGER.error('%s file must be present for psb_plan - but it is not there', os.path.join(self.working_directory,missense_filename))

        launch_parse_command = "psb_plan.py"
        LOGGER.info("Running: %s" % launch_parse_command)

        launch_parse_return = \
            subprocess.run(launch_parse_command, shell=False, encoding='UTF-8',
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        LOGGER.debug("%s finshed with exit code %s", launch_parse_command, launch_parse_return.returncode)
        os.chdir(save_cwd)


        return self.last_module_returncode

    def psb_launch(self) -> int:
        # shell out to 
        # $ singularity exec (psb_launch.py --nolaunch)
        #
        # Then run the generated ./launch_...case.py script
        save_cwd = os.getcwd()
        LOGGER.info("chdir(%s)", self.working_directory)
        os.chdir(self.working_directory)

        psb_launch_command = "psb_launch.py --nolaunch"

        # psb_launch_command = "export UDN=/dors/capra_lab/users/mothcw/VUStruct/; singularity exec --bind /dors/capra_lab $UDN/development.simg psb_launch.py --nolaunch"

        LOGGER.info("Running: %s" % psb_launch_command)

        psb_launch_return = \
            subprocess.run(psb_launch_command, shell=True, encoding='UTF-8',
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        LOGGER.debug("%s finshed with exit code %s", psb_launch_command, psb_launch_return.returncode)

        def line_count(s: str):
            if not s:
                return 0
            return s.count('\n') + 1
        LOGGER.debug("#stdout lines = %d   #stderr lines = %d", line_count(psb_launch_return.stdout), line_count(psb_launch_return.stderr))

        self.last_module_launched = "launch"
        self.last_module_returncode = psb_launch_return.returncode

        generated_launch_command = "./launch_%s.py" % self.external_user_prefix
        if psb_launch_return.returncode != 0: 
            LOGGER.warn("Exiting psb_launch() without running: %s", generated_launch_command)
        else:
            # Crank up slurm! with the generated
            # I.E. Run the launch_...py that was generated by psb_launch.py from 
            # inside the container
            LOGGER.info("Running: %s" % psb_launch_command)
            all_jobs_slurm_launch_return = \
                subprocess.run(generated_launch_command, shell=True, encoding='UTF-8',
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            LOGGER.debug("%s finshed with exit code %s", \
                 generated_launch_command, \
                 all_jobs_slurm_launch_return.returncode)
            LOGGER.debug("#stdout lines = %d   #stderr lines = %d",  \
                 line_count(all_jobs_slurm_launch_return.stdout), \
                 line_count(all_jobs_slurm_launch_return.stderr))

            self.last_module_returncode = all_jobs_slurm_launch_return.returncode

        os.chdir(save_cwd)


        return self.last_module_returncode

    def psb_monitor(self) -> int:
        save_cwd = os.getcwd()
        LOGGER.info("psb_monitor: chdir(%s)", self.working_directory)
        os.chdir(self.working_directory)

        psb_monitor_command_line = "psb_monitor.py"
        LOGGER.info("psb_monitor: %s", psb_monitor_command_line)

        self.last_psb_monitor_start = datetime.now().isoformat()
        launch_psb_monitor_return = \
            subprocess.run(psb_monitor_command_line, shell=False,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        LOGGER.info("psb_monitor finished with %d", launch_psb_monitor_return.returncode)
        self.last_psb_monitor_end = datetime.now().isoformat()
        os.chdir(save_cwd)

        return launch_psb_monitor_return.returncode

    def psb_rep(self, last_flag:bool=False, seconds_remaining:int=4*24*3600) -> int:
        global case_uuids_needing_website_refresh
        save_cwd = os.getcwd()
        LOGGER.info("psb_rep: chdir(%s)", self.working_directory)
        os.chdir(self.working_directory)

        # When we run psb_rep, we want to --transform out the uuid from the filenames
        # And we want the external_user filename prefixes to be removed as well.
        psb_rep_command_line = ['psb_rep.py','--tar_only', '--strip_uuid', self.case_uuid, '--web_case_id', self.case_id]
        if not last_flag:
            psb_rep_command_line.append('--embed_refresh')
            if seconds_remaining > 0:
                psb_rep_command_line.extend(['--seconds_remaining','%d' % seconds_remaining])

        LOGGER.info("psb_rep: Starting %s", ' '.join(psb_rep_command_line))

        self.last_psb_rep_start = datetime.now().isoformat()
        launch_psb_rep_return = \
            subprocess.run(psb_rep_command_line, shell=False,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        LOGGER.info("psb_rep: Returned %d from %s", launch_psb_rep_return.returncode, ' '.join(psb_rep_command_line))
        if launch_psb_rep_return.returncode != 0:
            LOGGER.warn("psb_rep failure %d log is %s", launch_psb_rep_return.returncode, os.path.join(os.getcwd(),"psb_rep.log"))
        self.last_psb_rep_end = datetime.now().isoformat()
        os.chdir(save_cwd)
        case_uuids_needing_website_refresh.add(self.case_uuid)
        LOGGER.info("case_uids needing refresh = %s", '\n'.join(case_uuids_needing_website_refresh))

        return launch_psb_rep_return.returncode

def monitor_report_loop(vustruct_case: VUstructCaseManager):
    # For 4 days loop, creating a new report every ?30? minutes
    # as the cluster churns away.  This is a separate functio to support
    # entry into _just_ this function (restart report creation)
    vustruct_case.timeout_time = datetime.now() + timedelta(days=config['report_timeout_days'])

    last_psb_rep = False

    last_psb_retcd  = 0
    report_sleep_seconds = 60 * config['report_interval_initial_minutes']

    while not last_psb_rep:
        time_remaining = (vustruct_case.timeout_time - datetime.now()).total_seconds()
        if time_remaining <= 0:
            # Time's up - this is the last report run
            time_remaining = 0
            last_psb_rep = True

        psb_monitor_retcd = vustruct_case.psb_monitor()

        if psb_monitor_retcd != 0 :
            LOGGER.error("psb_monitor.py for %s failed with %s" % (vustruct_case.case_uuid, launch_psb_rep_retcd))

        # We will build a new report if psb_monitor was successful OR
        # if this is the very last psb_rep run
        if (psb_monitor_retcd == 0) or last_psb_rep :
            psb_rep_retcd = vustruct_case.psb_rep(last_flag=last_psb_rep,seconds_remaining=time_remaining)
  
            if psb_rep_retcd == 0:
                case_uuids_needing_website_refresh.add(vustruct_case.case_uuid)
            else:
                LOGGER.error("psb_rep for %s failed with %s" % (vustruct_case.case_uuid, launch_psb_rep_retcd))
        time.sleep(report_sleep_seconds) # Wait a few minutes before trying again 
        report_sleep_seconds = min(
            report_sleep_seconds * config['report_interval_stretch_factor'],
            config['report_interval_max_minutes'] * 60.0)

    return last_psb_retcd



def launch_vustruct_case_thread(vustruct_case: VUstructCaseManager) -> subprocess.CompletedProcess:
    """
    This function is launched via the thread API.  "int he background" it will run the various components
    of the entire pipeline, and return in case of terminal failure, or completion.
    Along the way, the website files will slated for updated
    """

    cases_dict[vustruct_case.case_uuid] = vustruct_case

    # If we need to preprocess (parse) an input to create a missense file
    # Jump on that first and set last_module_launched to 'preprocess'
    # If no preprocessing is required by the user, then we flow through 
    # to psb_plan and process the missense.csv files...
    vustruct_case.last_module_returncode = 0
    if vustruct_case.data_format == 'Vanderbilt UDN Case Spreadsheet':
        vustruct_case.last_module_launched = 'preprocess'
        vustruct_case.last_module_returncode = vustruct_case.parse_udn_report()
    elif vustruct_case.data_format == 'VCF GRCh38':
        vustruct_case.last_module_launched = 'preprocess'
        vustruct_case.last_module_returncode = vustruct_case.vcf2missense()
    elif vustruct_case.data_format == 'VCF GRCh37 (runs liftover)':
        vustruct_case.last_module_launched = 'preprocess'
        vustruct_case.last_module_returncode = vustruct_case.vcf2missense_with_liftover()


    # Whatever we might have done with genomic starting coordinates... proceed!
    if vustruct_case.last_module_launched == 'preprocess':
        vustruct_case.last_module_error_message = ""
        if vustruct_case.last_module_returncode != 0:
            vustruct_case.last_module_error_message = "LATER Add error message for parse_udn_report please"
        # After we complete the parse step, we definitely need to run the psb_rep program
        # so that the OTHER flask program will install the growing website if it comes to that.
        # HOWEVER, the psb_rep program must note if the parse process has terminated with bad exit code
        # and take that into consideration.  IT may not as of yet.
        launch_psb_rep_retcd = vustruct_case.psb_rep()

        if launch_psb_rep_retcd != 0:
            vustruct_case.last_report_error = "psb_rep.py for %s failed with %s" % \
                                             (vustruct_case.case_uuid, launch_psb_rep_retcd)


        if vustruct_case.last_module_returncode != 0:
            return None
    

    # Whether we parsed genomic coordinates or not, the NEXT thing
    # to do is run psb_plan.
    vustruct_case.last_module_launched = 'plan'
    vustruct_case.last_module_returncode = vustruct_case.psb_plan()
    vustruct_case.last_module_error_message = ""
    if vustruct_case.last_module_returncode != 0:
        vustruct_case.last_module_error_message = "LATER Add error message for psb_plan please"

    # Error or not with psb_plan...
    # Now that we did psb_plan, we need to report
    launch_psb_rep_retcd = vustruct_case.psb_rep()

    if launch_psb_rep_retcd != 0:
        vustruct_case.last_report_error = "psb_rep.py for %s failed with %s" % \
                                         (vustruct_case.case_uuid, launch_psb_rep_retcd)


    if vustruct_case.last_module_returncode != 0:
        return None

    # Now we need to do launch (create the .slurm files as first step)
    # and the same code runs the generated ./launch...py outside the container

    vustruct_case.last_module_launched = 'launch'
    vustruct_case.last_module_returncode = vustruct_case.psb_launch()
    vustruct_case.last_module_error_message = ""
    if vustruct_case.last_module_returncode != 0:
        vustruct_case.last_module_error_message = "LATER Add error message for psb_launch please"

    # Error or not with psb_launch ...
    # Now that we did psb_plan, we need to report
    launch_psb_rep_retcd = vustruct_case.psb_rep()

    if launch_psb_rep_retcd != 0:
        vustruct_case.last_report_error = "psb_rep.py for %s failed with %s" % \
                                         (vustruct_case.case_uuid, launch_psb_rep_retcd)
    if vustruct_case.last_module_returncode != 0:
        return None

    monitor_report_loop(vustruct_case)

    return 0 # For now - not sure what to do after 4 days psb_rep_retcd

# def launch_vustruct_case(uuid_str: str):
#     psb_plan_command = \
#         ("export UDN=/dors/capra_lab/users/mothcw/VUStruct/; " +
#          "mkdir -p cd $UDN/fakecase_%s; " % uuid_str) + \
#         "singularity exec --bind /dors/capra_lab/ $UDN/development.simg psb_plan.py"
#     print("UUID = %s" % uuid_str )
#     print("running: %s" % psb_plan_command)
#     psb_plan_run_return = \
#         subprocess.run(psb_plan_command, shell=True,
#                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

AUTH_BASIC_CREDENTIAL="vustruct_password"

def require_authn(func):

    @wraps(func)
    def protected_endpoint(*args, **kwargs):
        # print("%s" % request.headers)
        # We only accept posts from localhost OR via the special Trafik gateway setup by Eric Appelt
        if request.headers.get('X-Forwarded-Host', '') == 'api.vgi01.accre.vanderbilt.edu':
            # Appears that this request coming through trafik - but let's be sure!
            if request.headers.get('X-Forwarded-Server','') != 'vm-infr-traefik.vampire':
                abort(401)
            elif request.headers.get('X-Forwarded-For','') != '129.59.141.42':
                abort(401)
        #Only other option is that we are in development mode on local host
        elif request.headers.get('Host','') != 'localhost:3080':
            abort(405)

        # print("%s" % request.headers.get('Authorization', '') )
        if request.headers.get('Authorization', '') != \
    AUTH_BASIC_CREDENTIAL:
            abort(401)
        return func(*args, **kwargs)

    return protected_endpoint

@app.route('/get_uuid', methods=['POST'])
@require_authn
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

    LOGGER.info ("In /get_uuid: Received data: %s json: %s" % (request.data, request.json))
    # Launch psb_plan.py

    uuid_str = str(uuid.uuid4())

    LOGGER.info("POST to /getuuid returning %s", uuid_str)

    return jsonify({"case_uuid": uuid_str})  # 'Hello ' #  + rq.get('name', 'No name'))

@app.route('/launch_vustruct', methods = ['POST', 'GET'])
@require_authn
def launch_vustruct():
    """
    There are a variety of ways to launch a new VUstruct case - 
    Some start frmo the missense.csv.  Others from a spreadsheet format.
    
    An example curl command would be:
    curl -X POST http://localhost:5000/launch_vustruct \ 
        -H 'Content-Type: application/json' \
        -d'{"data_format": "Vanderbilt UDN Case Spreadsheet",
            "case_uuid": "12435",
            "case_id": "test",
            "upload_file_URI": "file:////dors/capra_lab/users/mothcw/VUStruct/fakecase0/fakecase0.xlsx"}'
    """

    LOGGER.debug ("/launch_vustruct: request.json=\n%s", json.dumps(request.json,indent=3))

    for required_key in ['case_id', 'case_uuid', 'data_format']:
        if not required_key in request.json:
            return "Request JSON is missing required key: " + required_key

    vustruct_case = VUstructCaseManager(
        case_id=request.json['case_id'],
        case_uuid=request.json['case_uuid'])
    vustruct_case.data_format = request.json['data_format']
    vustruct_case.missense_csv = request.json['missense_csv']

    if 'upload_file_URI' in request.json:
        vustruct_case.upload_file_URI = request.json['upload_file_URI']

    # This so far has been very quick.  Now run the background thread to run the entire pipeline
    # Return to caller immediately of course.

    psb_plan_thread = threading.Thread(target=launch_vustruct_case_thread, args=(vustruct_case,))
    psb_plan_thread.start()

    # Send back the entire json that was send to us AnD add a launch_status to the gemisch
    # json_to_return_after_launch = request.json | {"launch_status": "VUstruct now running on accre cluster"}
    vustruct_case_json = jsonify(vars(vustruct_case))
    LOGGER.debug("Sending back this json %s", vustruct_case_json )

    return vustruct_case_json

@app.route('/dev_monitor_report_loop', methods = ['POST'])
@require_authn
def dev_monitor_report_loop():
    """
    For development purposes, it can be very helpful to loop on psb_monitor.py and psb_rep.py after 
    Therefore, there is this target to permit a "restart" of course
    
    An example curl command would be:
    curl -X POST http://localhost:5000/report_vustruct \ 
        -H 'Content-Type: application/json' \
        -d'{"case_uuid": "12435",
            "case_id": "test"
    """

    LOGGER.debug ("/dev_report_vustruct: request.json=%s", str(request.json))

    for required_key in ['case_id', 'case_uuid']:
        if not required_key in request.json:
            return "Request JSON is missing required key: " + required_key

    vustruct_case = VUstructCaseManager(
        case_id=request.json['case_id'],
        case_uuid=request.json['case_uuid'])

    cases_dict[vustruct_case.case_uuid] = vustruct_case

    vustruct_case.last_module_launched = 'launch'
    vustruct_case.last_module_returncode = 0
    vustruct_case.last_module_error_message = ""

    # This so far has been very quick.  Now run the background thread to run the entire pipeline
    # Return to caller immediately of course.

    dev_report_vustruct_thread = threading.Thread(target=monitor_report_loop, args=(vustruct_case,))
    dev_report_vustruct_thread.start()

    # Send back the entire json that was send to us AnD add a launch_status to the gemisch
    # json_to_return_after_launch = request.json | {"launch_status": "VUstruct now running on accre cluster"}
    vustruct_case_json = jsonify(vars(vustruct_case))
    LOGGER.debug("Sending back this json %s", vustruct_case_json )

    return vustruct_case_json


@app.route('/cases_needing_refresh', methods = ['POST', 'GET'])
@require_authn
def cases_needing_refresh():
    global case_uuids_needing_website_refresh
    # Reset the listof jobs needing website refresh so we don't keep going and going
    cases_needing_website_refresh = [vars(cases_dict[uuid]) for uuid in case_uuids_needing_website_refresh]
    # Clear out the set of jobs needing website refresh
    case_uuids_needing_website_refresh = set()
    return jsonify(cases_needing_website_refresh)

@app.route('/peek_cases_needing_refresh', methods = ['POST', 'GET'])
@require_authn
def peek_cases_needing_refresh():
    global case_uuids_needing_website_refresh
    cases_needing_website_refresh = [vars(cases_dict[uuid]) for uuid in case_uuids_needing_website_refresh]
    # With this peek_ call for testing, I do not clear out the set
    LOGGER.info("peek_cases_needing_refresh count=%d", len(cases_needing_website_refresh))
    return jsonify(cases_needing_website_refresh)


@app.route('/get_all_cases', methods = ['POST'])
@require_authn
def get_all_jobs():
    return jsonify([vars(active_vustruct) for active_vustruct in cases_dict.values()])

@app.route('/add_case', methods = ['POST'])
@require_authn
def add_job():
    new_case_info = request.json
    LOGGER.info("add_case request json = %s" % new_case_info)
    active_vustruct = VUstructCaseManager(new_case_info['case_id'], new_case_info['case_uuid'])
    cases_dict[new_case_info['case_uuid']] = active_vustruct


    return jsonify(vars(active_vustruct))
