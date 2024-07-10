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
#     vustruct_flask is executed outside the container.
#
# curl samples can be found with the various calls, just search for "curl"
# in the source code below
#
# The class CaseManager keeps the details of each
# progressing case organized - and is periodically saved in case 
# for reloading in case of server restart
#
# FILESYSTEM INTERACTION
# To allow for stop and restart of vustruct_flask.py....
# vustruct_flask.py records bookkepping to the filesystem
#
# As cases are added to the global set of monitored cases,
# as well as when casea compelte...  
# the global "./vustruct_flask_active_cases.json" file is
# updated
#
# The progress each of running case is recorded in the
# working_directory of each running flask case under
# management

from __future__ import annotations

import uuid
import threading
import logging
import os
import sys
from io import StringIO
from datetime import datetime
from datetime import timedelta
import time
import json
import subprocess
import tempfile
import csv
import string
import pprint

import pandas as pd
import numpy as np

from functools import wraps
from flask import Flask
from flask import abort
from flask import request
from flask import jsonify
from flask import Blueprint
from flask import current_app
from flask import make_response
from werkzeug.exceptions import HTTPException
import werkzeug
werkzeug.serving._log_add_style = False
# requests is used to read the excel or missense filess uploaded to wordpress.
import requests
from requests_file import FileAdapter
from typing import Dict

# import time
from logging.handlers import RotatingFileHandler
from typing import Dict, Any, Tuple, List
from typing import TextIO
# Inspect to recover bsub_submit() and slurm_submit() source code for integration
# into generated final file
import inspect



###############################################
# Logging must be initiated in Flask
# before we include any other modules that
# engage the usual logging system

from logging.config import dictConfig

os.makedirs('log', exist_ok=True)

# Setup the logging.  By calling dictConfig, the flask app.logger
# will automatically pick up these settings
# This is the "way to go" with a flask application
LOG_FILENAME='log/vustruct_flask.log'
logging.config.dictConfig({
    'version': 1,
    'formatters': {
        'default': {
            'format': '%(asctime)s %(levelname)-4.4s [%(filename)20s:%(lineno)d] %(message)s',
            'datefmt': '%H:%M:%S'
        }
     },
    'handlers': {
        'wsgi': { # NOT USED AT THE MOMENT
            'class': 'logging.StreamHandler',
            'stream': 'ext://flask.logging.wsgi_errors_stream',
            'formatter': 'default'
        },
        'stderr': {
            'class': 'logging.StreamHandler',
            'stream': sys.stderr,
            'formatter': 'default'
        },
        'rotating_to_file': {
            'level': 'DEBUG',
            'class': "logging.handlers.RotatingFileHandler",
            'formatter': 'default',
            "filename": LOG_FILENAME,
            "maxBytes": 100000,
            "backupCount": 10,

        }
     },
    'root': {
        'level': 'INFO',
        'handlers': ['stderr', 'rotating_to_file']
    }
})
app = Flask(__name__)

app.logger.info("Flask app created with Flask(\'%s\')" % __name__)
app.logger.info("Logging level %d is additionally output to %s" % (app.logger.getEffectiveLevel(), LOG_FILENAME))
# if os.getenv("FLASK_ENV") == 'development':
assert app.logger.getEffectiveLevel() == logging.DEBUG

# from jinja2 import Environment, FileSystemLoader

import slurm
from bsub import bsub_submit

from psb_shared import psb_config
from psb_shared import psb_perms

# Flask programs cannot easily digest command line arguments - so we 
# the calling bash script sets environment variables in stead
vustruct_config_filename = os.getenv("VUSTRUCT_CONFIG")
if not vustruct_config_filename:
    errormsg = "VUSTRUCT_CONFIG environment variable must be defined before launch of %s " % __name__
    sys.exit(errormsg)
vustruct_userconfig_filename = os.getenv("VUSTRUCT_USERCONFIG")
if not vustruct_userconfig_filename:
    errormsg = "VUSTRUCT_USERCONFIG environment variable must be defined before launch of %s " % __name__
    sys.exit(errormsg)

# In order to use the exising config file parsing machinery, it was helpful
# to create this "fakeargs" little class
# print "Command Line Arguments"
class FakeArgs(object):
    pass
fakeargs = FakeArgs()
fakeargs.config = vustruct_config_filename
fakeargs.userconfig = vustruct_userconfig_filename
fakeargs.caseconfig = ''
# UDN="/dors/capra_lab/users/mothcw/VUStruct/"
# UDN = os.getenv("UDN")
# if not UDN:
#     errormsg = "UDN environment variable must be defined before launch of %s" % UDN
#     sys.exit(errormsg)

# These items moved to the global.config file
# config = {
#     'report_interval_initial_minutes': 1.0,
#     'report_interval_stretch_factor': 1.5,
#     'report_interval_max_minutes': 30,
#     'case_timeout_days': 4.0,
#     'flask_case_filename': 'flask_case_progress.json',
#     'monitored_caselist_file': '/dors/capra_lab/users/mothcw/VUStruct/vustruct_flask_monitored_caselist.txt'
# }

class GlobalActiveCases:

    # static variables shared among all GlobalActiveCases members
    save_filename = os.path.join(os.getcwd(), "GlobalActiveCases_vustruct_flask.json")
    global_lock = threading.Lock()
    with global_lock:
        # Given a uuid, all the requests can find the full VUStructFlask instance
        uuid_VUStructFlask_dict = {}

        # When vustruct_webupdate calls for updates
        # this set is returned to it
        uuids_needing_website_refresh = set()

        # When a case has run psb_rep.py for the last time...
        # Then it needs to be added to this set so that the
        # user can be informed that no more will be coming
        uuids_final_website_update = set()

        save_needed = True

    @staticmethod
    def add_flask_case(flask_case):
        with GlobalActiveCases.global_lock:
            GlobalActiveCases.uuid_VUStructFlask_dict[flask_case.case_uuid] = flask_case
        GlobalActiveCases.save_needed = True

    @staticmethod
    def get_flask_case_or_none(case_uuid: str):
        found_flask_case = None
        with GlobalActiveCases.global_lock:
            found_flask_case = GlobalActiveCases.uuid_VUStructFlask_dict.get(case_uuid, None)
        return found_flask_case

    @staticmethod
    def remove_flask_case_uuid(flask_case_uuid: str):
        with GlobalActiveCases.global_lock:
            if flask_case_uuid in GlobalActiveCases.uuid_VUStructFlask_dict:
                del GlobalActiveCases.uuid_VUStructFlask_dict[flask_case_uuid]
        GlobalActiveCases.save_needed = True

    @staticmethod
    def save_active_cases() -> None:
        if not GlobalActiveCases.save_needed:
            return
        
        # Save list of all the active cases for reloading in case of 
        # vustruct_flask.py restart
        _active_cases = []
        with GlobalActiveCases.global_lock:
            for case_uuid in GlobalActiveCases.uuid_VUStructFlask_dict:
                _active_cases.append({
                    'case_id': GlobalActiveCases.uuid_VUStructFlask_dict[case_uuid].case_id,
                    'case_uuid': case_uuid
                    }) 

        _active_cases_json = json.dumps(
            _active_cases,
            indent=4)

        _tmp_filename = None 
        with tempfile.NamedTemporaryFile(
            mode='w+t',
            prefix=GlobalActiveCases.save_filename,
            suffix=".tmp",
            delete=False,
            ) as f:
                f.write(_active_cases_json)
                _tmp_filename = f.name
        os.rename(_tmp_filename, GlobalActiveCases.save_filename)
        GlobalActiveCases.save_needed = False
    
    @staticmethod
    def reload_from_save_file():
        app.logger.info("GlobalActiveCases: Attempting to restart monitoring of cases in %s" % \
            GlobalActiveCases.save_filename)

        try:
            with open(GlobalActiveCases.save_filename,'r+t') as f:
                _active_cases = json.load(f)
       
        except FileNotFoundError as ex:
            app.logger.warning("Unable to open %s", GlobalActiveCases.save_filename)
            _active_cases = []

        app.logger.info("reload_from_save_file.  %d active cases in %s" % (
            len(_active_cases), GlobalActiveCases.save_filename))

        for active_case in _active_cases:
            flask_case = CaseManager(
                active_case['case_id'],
                active_case['case_uuid'])

            app.logger.info("Attempting reload of case %s_%s" % (
                flask_case.case_id,
                flask_case.case_uuid))

            try:
                flask_case.load_vustruct_csv()
            except Exception as ex:
                app.logger.error("Could not load vustruct.csv file %s %s:\n%s",      
                    active_case['case_id'],
                    active_case['case_uuid'],
                    str(ex))
                continue

            try:
                flask_case.load_from_json()
            except Exception as ex:
                app.logger.error("Could not load json file for %s %s:\n%s",      
                    active_case['case_id'],
                    active_case['case_uuid'],
                    str(ex))
                continue

            flask_case.load_psb_monitor_json()
            if flask_case.psb_monitor_results_dict:
                flask_case.psb_load_jobids_from_psb_monitor_results()
            else:
                flask_case.psb_load_jobids_from_workstatus_file()


            for (unique_id, job_id_int, array_id_int) in flask_case.launched_jobs:
                if array_id_int is not None:
                    job_key = "%s_%s" % (job_id_int, array_id_int)
                else:
                    job_key = str(job_id_int)
                    # import pdb; pdb.set_trace()
                    temp = global_squeue_df 


            app.logger.info("Case %s_%s to be reactivated for monitoring and reporting" %
                (active_case['case_id'],
                active_case['case_uuid']))
            # Add back this case to the global_cases_dict, in case it has been lost

            GlobalActiveCases.add_flask_case(flask_case)
        
            flask_case.last_module_launched = 'launch'
            flask_case.last_module_returncode = 0
            flask_case.last_module_error_message = ""
        
            # This so far has been very quick.  Now run the background thread to run the entire pipeline
            # Return to caller immediately of course.
        
            dev_report_vustruct_thread = threading.Thread(target=monitor_report_loop, args=(flask_case,False))
            dev_report_vustruct_thread.start()
        
            # Great we loaded the case... now we have to think about what comes next - restarting a thread likely required

    @staticmethod
    def add_website_refresh_needed(flask_case_uuid: str):
        with GlobalActiveCases.global_lock:
            GlobalActiveCases.uuids_needing_website_refresh.add(flask_case_uuid)
        app.logger.info("case_uids needing refresh = %s", '\n'.join(GlobalActiveCases.uuids_needing_website_refresh))

    @staticmethod
    def mark_final_website_update(flask_case_uuid: str) -> None:
        with GlobalActiveCases.global_lock:
            GlobalActiveCases.uuids_final_website_update.add(flask_case_uuid)
        app.logger.info("case_uids final websites = %s", '\n'.join(GlobalActiveCases.uuids_final_website_update))

    @staticmethod
    def cases_needing_website_refresh(clear = True) -> list[dict]:
        """
        Called by a request handler - so must be very fast
        """

        cases_needing_website_refresh = []

        with GlobalActiveCases.global_lock:
            # With the lock on, it would be disaster to write to file system or such
            for case_uuid in GlobalActiveCases.uuids_needing_website_refresh:
                # If somehow we've lost the case from our FlaskDict then we can't use that one
                flask_case = GlobalActiveCases.uuid_VUStructFlask_dict.get(case_uuid, None)
                if flask_case is not None:
                    # We don't return everything to vustruct_webdate.py - just what it needs to
                    # find the report tar files and unpack them
                    cases_needing_website_refresh.append({
                       'case_id': flask_case.case_id,
                       'case_uuid': case_uuid,
                       'working_directory': flask_case.working_directory
                       })

                    # If this webupdate is the LAST one (after what could have been days)
                    # Then we are performing a final update and this is important
                    if clear and case_uuid in GlobalActiveCases.uuids_final_website_update:
                        del GlobalActiveCases.uuid_VUStructFlask_dict[case_uuid]
                        # We dare not write to file system now - but we could stand to do it later
                        GlobalActiveCases.save_needed = True
                        GlobalActiveCases.uuids_final_website_update.remove(case_uuid)

            # Clear out the set of jobs needing website refresh
            if clear:
                GlobalActiveCases.uuids_needing_website_refresh = set()

        return cases_needing_website_refresh 

app.logger.info("Config file locations:\n%s" % pprint.pformat(vars(fakeargs)))

if not os.path.exists(fakeargs.config):
    app.logger.critical("Global config file not found: " + fakeargs.config)
    sys.exit(1)

required_config_items = ['output_rootdir', 'collaboration', 'cluster_user_id']

config, config_dict = psb_config.read_config_files(fakeargs, required_config_items)

VUStructFlask_config = dict(config.items("VUStructFlask"))

VUStructFlask_floats = [
    'report_interval_initial_minutes',
    'report_interval_stretch_factor',
    'report_interval_max_minutes',
    'case_timeout_days']

#Before we get cranked up, check that the convig items are there and convertable to float
for numerical_item in VUStructFlask_floats:
    VUStructFlask_config[numerical_item] = float(
        VUStructFlask_config[numerical_item])


# print("Level is %s " % os.getenv("FLASK_ENV"))
# else:
#     app.logger.setLevel(logging.INFO)
# app.logger.setFormatter(log_formatter)

UDN = os.path.join(config_dict['output_rootdir'], config_dict['collaboration'])
app.logger.info("%s running web-input cases in %s", __file__, UDN)

def kill_vustruct_flask(caller_info: str = "No caller info provided"):

    print("Shutting down")
    signal_number = 15
    our_pid = os.getpid()
    app.logger.info(f"{caller_info}\nShutting down with kill(pid={our_pid},signal={signal_number})")
    os.kill(our_pid,signal_number)

# We cann squeue and load the slurm status of all running jobs
# The status is indexed on job_key, the slurm concatenation
# of a job id, and then if array id, _array_id
# So that for any _particular_ case, we can quickly
# get the squeue 
global_squeue_df = pd.DataFrame()
def global_squeue_df_refresh() -> None:
    global global_squeue_df
    _slurm_user_id = config_dict['cluster_user_id']

    # Shell out to squeue, to gather status of all our running jobs
    # across all cases
    parse_return = slurm.slurm_squeue(_slurm_user_id)

    if parse_return.returncode != 0:
        app.logger.error("Unable to run squeue return_code=%d:\n%s", parse_return.returncode, parse_return.stderr)
        return

    # Set the index to job_key
    _latest_squeue_df = slurm.flatten_squeue_stdout_to_df(parse_return.stdout).set_index('job_key')
    app.logger.info("%d active slurm jobs for user %s", len(_latest_squeue_df), config_dict['cluster_user_id'])
    global_squeue_df = _latest_squeue_df

def global_squeue_df_refresh_thread() -> None:
    while True:
        _sleep_seconds = int(VUStructFlask_config['squeue_interval_minutes']) * 60
        app.logger.info("Next squeue refresh in %d seconds" % _sleep_seconds)
        time.sleep(_sleep_seconds)
        global_squeue_df_refresh()
    
global_squeue_df_refresh()    

squeue_refresh_thread = threading.Thread(target=global_squeue_df_refresh_thread)
squeue_refresh_thread.start()


# An unhandled exception demands a total shutdown of this process
# after log attempt.  Presumably the service manager will restart
# us
# See https://stackoverflow.com/questions/64976156/flask-generic-exception-handler-unable-to-catch-405-error
# for details about not choking on HTTP errors from evil posters
@app.errorhandler(Exception)
def handle_exception(e):
    # Return bad HTTP errors or goofed  json requests
    # If we get a bad POST Json, keep going in that case too!
    if isinstance(e, HTTPException):
        response = e.get_response()
        payload = json.dumps({
            "code": e.code,
            "name": e.name,
            "description": e.description,
        })
        response.data = f"{payload}\n"
        response.content_type = "application/json"
        return make_response(response, e.code)

    if isinstance(e, json.decoder.JSONDecodeError):
        return make_response({
            "JSONDecodeError": str(e),
            "msg": e.msg,
            "doc": e.doc,
            "pos": e.pos,
            "lineno": e.lineno,
            "colno": e.colno
            }, 999)
    # now handle non-HTTP exceptions only
    # response_body = {"code": 500, "message": "unhandled exception occurred", "details": str(e)}
    # return make_response(response_body, 500)

    app.logger.exception("Unhandled exception %s", str(e))
    kill_vustruct_flask("global exception")
    return "Unhandled global exception %s" % str(e)

class CaseManager:
    # We need to be able to save case management information to disk when this process is terminated
    # and restore it on restart.
    # We'll leverage python's getattr and setattr machinery to speed transition to and from an
    # external dictionary.  Time and other non-json elements will get special treatment

    _static_vustruct_member_strings_and_ints = [
        'case_id',
        'case_uuid',
        'last_module_launched',   # One of preprocess, plan, launch
        'last_module_returncode', # If non-zero, then that means we died at that point and did not coninue to launch
        'upload_file_URI',
        'last_psb_rep_start',
        'last_psb_rep_end',
        'last_psb_monitor_end',
        'data_format',
        'init_time',
        'report_sleep_seconds',
        'save_cwd',
        'launched_jobs']

    def __init__(self, case_id: str, case_uuid: str) -> None:
        for member_name in CaseManager._static_vustruct_member_strings_and_ints:
            setattr(self, member_name, "")

        self.case_id = case_id
        self.case_uuid = case_uuid
        # As we launch processes, we populate this with
        # the vustruct "proprocess", "plan", "launch", "report"
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
        # external_user_prefix replaced with property functions to eliminate
        # refundancy.  It is NOT needed to be part of JSON returned by restAPIs
        # self.external_user_prefix = "external_user_%s_%s" % (case_id, case_uuid)

        self.working_directory = os.path.join(UDN, self.external_user_prefix)
        self.data_format = ""

        # After launch, we (again) generate a report
        # Then we extend the refresh interval - because less and less changes 
        # Out to a maximum of the config file entry.
        # This here is our starting entry.
        self.report_sleep_seconds = 60 * VUStructFlask_config['report_interval_initial_minutes']


        self.init_time = datetime.now().isoformat()

        # Init to empty.  This gets loaded  from the case's final input file
        # ...._vustruct.csv, after successful preprocess or direct uplodate
        self.vustruct_csv_df = pd.DataFrame()

        # After the launch, launched jobs are triples that contain the job id and array id and the uniquekey of the jobs
        self.launched_jobs = [] # pd.DataFrame(columns=['unique_key', 'jobid', 'arrayid'])

        # After each call to psb_monitor, we need to harvest the status of all case-wide jobs
        self.psb_monitor_results_dict = {}

    @property
    def external_user_prefix(self) -> str:
        return "external_user_%s_%s" % (self.case_id, self.case_uuid)

    def to_dict(self) -> Dict:
        """ Serialize a CaseManager instance to JSON ready dictionary """
        as_dict = {}
        for member_name in CaseManager._static_vustruct_member_strings_and_ints:
            as_dict[member_name] = getattr(self, member_name)
        return as_dict

    @classmethod
    def from_dict(cls, dict) -> CaseManager:
        """ Build a CaseManager from a dictionary read from json """
        new_CaseManager = cls(dict['case_id'],dict['case_uuid'])
        for member_name in CaseManager._static_vustruct_member_strings_and_ints:
            if member_name in dict:
                setattr(new_CaseManager, member_name, dict[member_name])
        return new_CaseManager

    @property
    def fullpath_case_progress_filename(self):
        return  os.path.join(
            self.working_directory, 
            "%s_flask_progress.json" % self.external_user_prefix)

    def save_as_json(self):
        """
        Save the progress of a vustruct case so that if vustruct_flask is restarted, it can resume
        where it left off
        """
        flask_case_json = json.dumps(
            self.to_dict(),
            indent=4)
    
        _tmp_filename = None 
        with tempfile.NamedTemporaryFile(
            mode='w+t',
            prefix=os.path.join(
                self.working_directory, 
                "%s_progress" % (self.case_id)),
            suffix=".tmp",
            delete=False,
            ) as f:
                f.write(flask_case_json)
                _tmp_filename = f.name
    
    
        os.rename(_tmp_filename, self.fullpath_case_progress_filename)
        app.logger.info("VUStruct case progress saved to %s" % self.fullpath_case_progress_filename )

    def load_from_json(self):
        """
        Load the progress of a vustruct case, presumably after restart of vustruct_flask
        """

        app.logger.info("VUStruct case progress loading from %s" % self.fullpath_case_progress_filename )
        with open(self.fullpath_case_progress_filename,'r') as f:
            flask_case_as_dict = json.load(f)

        self.from_dict(flask_case_as_dict)
    
    def create_failure_website(self, body_html: str, website_filelist: list[str]):
        # When something goes terribly wrong in flask processing, flask can create 
        # the website update itself - but it is not pretty
        with open("index.html","w") as f:
            f.write(f"""<html>\
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
  "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html lang="en-us" xmlns="http://www.w3.org/1999/xhtml" >
<head>
<meta charset="utf-8">

    <title>VUStruct case {self.case_id}</title>
</head>
<body>
{body_html}
</body>
</html>
""")



        # The temp filename has a hideous UTC timestamp on it

        website_filelist_filename = os.path.join(UDN,self.external_user_prefix,
                                             "%s_website_files.list" % self.external_user_prefix) 

        with open(website_filelist_filename, "w") as f:
            f.write('\n'.join((os.path.relpath(
                website_file,  # os.path.realpath(website_file),
                os.path.realpath(UDN))
                for website_file in website_filelist)))
        app.logger.info("A filelist for creating a website is in %s", website_filelist_filename)


        website_tar_temp_filename = "%s_%s.tar.gz" % (
            self.external_user_prefix,
            datetime.now().replace(microsecond=0).isoformat().translate(str.maketrans('','',string.punctuation)))

        website_tar_filename = "%s.tar.gz" % self.external_user_prefix

        search_string = os.path.join(self.external_user_prefix, self.external_user_prefix)
        replace_string = f"{self.external_user_prefix}/{self.case_id}"
        tar_transformer = f"--transform 's[{search_string}[{replace_string}[g' --show-transformed-names"
        tar_maker = 'cd ..; tar cvzf %s --files-from %s %s --mode=\'a+rX,go-w,u+w\' > %s.stdout; cd -; mv %s %s; mv %s.stdout %s.stdout' % (
            os.path.join(self.external_user_prefix, website_tar_temp_filename),
            os.path.join(self.external_user_prefix, website_filelist_filename),
            tar_transformer,
            os.path.join(self.external_user_prefix, website_tar_temp_filename), # Last file before cd -
            website_tar_temp_filename,website_tar_filename, # The temp and final tar files
            website_tar_temp_filename,website_tar_filename) # The temp and final stdout files
        app.logger.info("Creating .tar website file with: %s", tar_maker)
        subprocess.call(tar_maker, shell=True)

    def fetch_input_file(self,file_extension: str) -> str:
        """
        Retury a string if this fails because of problems with wordpress caller
        """
        if not self.upload_file_URI:
            self.create_failure_website(body_html=\
"""<H1>VUStruct failure while transferring your input file</H1>
<p>Please use your browser back button.  Then, reselect your input file</p>""",
        				website_filelist = ['index.html'])

            return "Input file URI did not transfer"

        app.logger.info("Attempting to load uploaded file from %s", self.upload_file_URI)

        # Usually we "get" the spreadsheet from wordpress using http
        # It is convenient to be able to test all this with local files
        _session = requests.Session()
        _session.mount('file://', FileAdapter())
        response = _session.get(self.upload_file_URI)

        # Write out the uploaded spreadseet or vcf in the final directory where it will be run
        final_VUStruct_input_filename = os.path.join(
              self.working_directory,
             "%s.%s" % (self.external_user_prefix, file_extension))

        with open(final_VUStruct_input_filename, "wb") as f:
            f.write(response.content)

        return None

    def run_coordinates_parser_command_line(self, parser_python_filename: str) -> subprocess.CompletedProcess:
        os.makedirs(self.working_directory, exist_ok=True)
        app.logger.info("Running: %s in \n%s", parser_python_filename, self.working_directory)

        parse_return = \
            subprocess.run(parser_python_filename, shell=True, encoding='UTF-8',
                       cwd=self.working_directory,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        app.logger.debug("%s finshed with exit code %d", parser_python_filename, parse_return.returncode)
        return parse_return




    def parse_udn_report(self) -> int:
        parse_return = self.run_coordinates_parser_command_line("parse_udn_report.py")
        return parse_return.returncode

    def parse_wustl(self) -> int:
        parse_return = self.run_coordinates_parser_command_line("parse_wustl.py")
        return parse_return.returncode

    def vcf2vustruct(self) -> int:
        parse_return = self.run_coordinates_parser_command_line("vcf2vustruct.py")
        return parse_return.returncode

    def vcf2vustruct_with_liftover(self) -> int:
        parse_return = self.run_coordinates_parser_command_line("vcf2vustruct.py --liftover")
        return parse_return.returncode
 
    @property
    def vustruct_filename(self) -> str:
        return os.path.join(self.working_directory, "%s_vustruct.csv" % self.external_user_prefix)

    def load_vustruct_csv(self) -> None:
        """
        following successful psb_plan, it is very helpful to simply have the
        entire vustruct.csv loaded as a dataframe.  That is invariant
        over the entire run anyway

        """

        with open(self.vustruct_filename,'r') as f:
            vustruct_csv_data = f.read()
            self.vustruct_csv_df = pd.read_csv(StringIO(vustruct_csv_data), sep=',', index_col=None, keep_default_na=False, encoding='utf8',
                                       comment='#', skipinitialspace=True)

    def psb_plan(self) -> int:
        if not os.path.exists(self.vustruct_filename):
            app.logger.error('%s file must be present for psb_plan - but it is not there', self.vustruct_filename)
            return 1

        launch_parse_command = "psb_plan.py"
        app.logger.info("Running: %s" % launch_parse_command)

        launch_parse_return = \
            subprocess.run(launch_parse_command, shell=False, encoding='UTF-8',
                       cwd = self.working_directory,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        app.logger.debug("%s finshed with exit code %s", launch_parse_command, launch_parse_return.returncode)


        return self.last_module_returncode

    @staticmethod
    def _parse_job_id_array_id_from_row(job_status_row) -> Tuple[int, int]:
        job_id_int = None
        jobid_fail = None
        try:
            job_id_int = int(job_status_row['jobid'])
            if job_id_int <= 0:
                jobid_fail = "Zero job id in %s" % job_status_row
          
        except Exception as ex:
            jobid_fail = "After launch, could not convert job id in %s" % job_status_row
    
        if jobid_fail:
            app.logger.critical(jobid_fail)
            return (None, None)
    
        array_id_int = None
        arrayid_fail = None
        # It's fine for arrayid to be missing - but NOT to be 
        # a weird number if it is there
        # 0 is a totally legit arrayid - so we have to take care
        # to not think of that as null/missing
        if job_status_row['arrayid']:
            try:
                # Offen, array_id is not used - and will be nan
                # Because pandas dataframes cannot store None in int columns
                # I store the array ID as a float... and that means a bit of
                # type conversion to get back to the int.
                array_id_int = int(float(job_status_row['arrayid']))
                if array_id_int < 0:
                    arrayid_fail = "Bad array ID in %s" % job_status_row
            except Exception as ex:
                arrayid_fail = "Could not convert array id [%s] from data:\n%s" % (
                    str(ex), job_status_row)
    
        if arrayid_fail:
            app.logger.critical(arrayid_fail)
            return (None, None)
    
        return (job_id_int, array_id_int)

    def psb_load_jobids_from_psb_monitor_results(self) -> bool:
        # It is very convenient to load the jobs idfs from the .json
        # file emitted from psb_monitor.  If we have that....
        _launched_jobs = []
        app.logger.info("Attempting load_jobids_from_psb_monitor_results %s_%s" % (
            self.case_id,
            self.case_uuid))

        _psb_monitor_results_dict = self.psb_monitor_results_dict.copy()
        for gene_refseq_mutation in _psb_monitor_results_dict:
            # app.logger.info("NOW DOING %s" % gene_refseq_mutation)
            # for unique_key in _psb_monitor_results_dict[gene_refseq_mutation]:
            #     app.logger.info("1 %s %s %s" % (threading.get_ident(), unique_key, _psb_monitor_results_dict[gene_refseq_mutation] is None))
            # for unique_key in _psb_monitor_results_dict[gene_refseq_mutation]:
            #     app.logger.info("2 %s %s %s" % (threading.get_ident(), unique_key, _psb_monitor_results_dict[gene_refseq_mutation] is None))
            if _psb_monitor_results_dict[gene_refseq_mutation] is None:
                continue
            for unique_key in _psb_monitor_results_dict[gene_refseq_mutation]:
                # app.logger.info("3 %s %s %s" % (threading.get_ident(), unique_key, _psb_monitor_results_dict[gene_refseq_mutation] is None))
                job_row=_psb_monitor_results_dict[gene_refseq_mutation][unique_key]

                job_id_int, array_id_int = self._parse_job_id_array_id_from_row(job_row)

                if (job_id_int is None) or (array_id_int is None):
                    continue

                _launched_jobs.append(
                    [ unique_key,
                      job_id_int,
                      array_id_int ])

        self.launched_jobs = _launched_jobs

        app.logger.info("From psb_monitor.py %d job/array ids have been loaded" % len(self.launched_jobs))


    def psb_load_jobids_from_workstatus_file(self):
        # Called as part of psb_launch, this code harvests the job IDs from the workstatus files, and 
        # builds an in-memory database of ids which are monitored and reported to websites via restapi calls
        # This is rather a "legacy" function since psb_monitor gathers a nice JSON that we can use to
        # get all this after psb_monitor with one file read

        _launched_jobs = []

        app.logger.info("Attempting legacy load_jobids_from_workstatus files %s_%s" % (
            self.case_id,
            self.case_uuid))

        df_all_mutations = self.vustruct_csv_df.fillna('NA')

        # each of our gene entries in teh vustruct file will have its own jobIDs and array IDs to harvest
        for index,row in df_all_mutations.iterrows():
            if 'RefSeqNotFound_UsingGeneOnly' in row['refseq']:
                row['refseq'] = 'NA'
            mutation_dir = os.path.join(self.working_directory,"%s_%s_%s"%(row['gene'],row['refseq'],row['mutation']))
            if not os.path.exists(mutation_dir):  # python 3 has exist_ok parameter... 
                logging.critical("Skipping missing mutation directory %s which should have been created by psb_launch.py",
                   mutation_dir)
                continue
    
            workstatus_filename = "%s/%s_%s_%s_workstatus.csv"%(mutation_dir,row['gene'],row['refseq'],row['mutation'])

            workstatus_rows = []
            try:
                with open(workstatus_filename) as f:
                    reader = csv.DictReader(f,delimiter='\t')
                    for csv_row in reader:
                        workstatus_rows.append(csv_row)
            except FileNotFoundError as ex:
                workstatus_rows = [] # pd.DataFrame() # The empty dataframe will not be iterated below

            # Iterate through workstatus: Figure out slurm array indices
            for job_status_row in workstatus_rows:
                (job_id_int, array_id_int) = self._parse_job_id_array_id_from_row(job_status_row)

                if job_id_int is None or array_id_int is None:
                    continue
                 
                _launched_jobs.append(
                    [ job_status_row['uniquekey'], 
                      job_id_int,
                      array_id_int ])


        self.launched_jobs = _launched_jobs

        app.logger.info("From workstatus files %d job/array ids have been loaded" % len(self.launched_jobs))

        # app.logger.info("WOW - what did we get for the actual job ids????")
        # app.logger.info("The dataframe of job IDs and array IDs is\n%s" % str(self.launched_jobs)) # .to_string())
                 

    def psb_launch(self) -> int:
        # shell out to 
        # $ singularity exec (psb_launch.py --nolaunch)
        # Then run the generated ./launch_...case.py script

        # psb_launch_command = "export UDN=/dors/capra_lab/users/mothcw/VUStruct/; singularity exec --bind /dors/capra_lab $UDN/development.simg psb_launch.py --nolaunch"
        psb_launch_command = ['psb_launch.py', '--nolaunch']

        app.logger.info("Running: %s from directory:\n%s" % (psb_launch_command, self.working_directory))

        psb_launch_return = \
            subprocess.run(psb_launch_command, shell=False, encoding='UTF-8',
                       cwd = self.working_directory,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        app.logger.debug("%s finshed with exit code %s", psb_launch_command, psb_launch_return.returncode)
        if psb_launch_return.returncode != 0:
            app.logger.critical("stdout=\n%s\n---stderr=\n%s" % (psb_launch_return.stdout, psb_launch_return.stderr))

        def line_count(s: str):
            if not s:
                return 0
            return s.count('\n') + 1
        app.logger.debug("#stdout lines = %d   #stderr lines = %d", line_count(psb_launch_return.stdout), line_count(psb_launch_return.stderr))

        self.last_module_launched = "launch"
        self.last_module_returncode = psb_launch_return.returncode

        generated_launch_command = "./launch_%s.py" % self.external_user_prefix
        if psb_launch_return.returncode != 0: 
            app.logger.warn("Exiting psb_launch() without running: %s", generated_launch_command)
        else:
            # Crank up slurm! with the generated
            # I.E. Run the launch_...py that was generated by psb_launch.py from 
            # inside the container
            app.logger.info("Running: %s" % generated_launch_command)
            all_jobs_slurm_launch_return = \
                subprocess.run(generated_launch_command, shell=False, encoding='UTF-8',
                       cwd = self.working_directory,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            app.logger.debug("%s finshed with exit code %s", \
                 generated_launch_command, \
                 all_jobs_slurm_launch_return.returncode)
            app.logger.debug("#stdout lines = %d   #stderr lines = %d",  \
                 line_count(all_jobs_slurm_launch_return.stdout), \
                 line_count(all_jobs_slurm_launch_return.stderr))
            if all_jobs_slurm_launch_return.returncode != 0:
                 app.logger.critical("stdout=\n%s\n---stderr=\n%s" % (
                     all_jobs_slurm_launch_return.stdout, 
                     all_jobs_slurm_launch_return.stderr))

            self.last_module_returncode = all_jobs_slurm_launch_return.returncode

 

        return self.last_module_returncode

    def load_psb_monitor_json(self):
        psb_monitor_json_filename = os.path.join(self.working_directory, 
                           self.external_user_prefix + "_psb_monitor.json")
        try:
            with open(psb_monitor_json_filename, 'r+t') as json_f:
                # Break into 2 steps in case of exception
                _temp_psb_monitor_results_dict = json.load(json_f)
                self.psb_monitor_results_dict = _temp_psb_monitor_results_dict
                return self.psb_monitor_results_dict
        except Exception as ex:
            # This is bad  - but not enough to crash out entirely
            app.logger.critical("Can't read %s after psb_monitor.py" % psb_monitor_json_filename)
            return {}


    def psb_monitor(self) -> int:
        """
        Run psb_monitor.py, required before psb_rep.py is run
        This code also gathers the written json for flask as part of dynamic job
        status replies
        """

        psb_monitor_command_line = "psb_monitor.py"
        app.logger.info("psb_monitor: %s", psb_monitor_command_line)

        self.last_psb_monitor_start = datetime.now().isoformat()
        launch_psb_monitor_return = \
            subprocess.run(psb_monitor_command_line, shell=False,
                       cwd=self.working_directory,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        app.logger.info("psb_monitor finished with %d", launch_psb_monitor_return.returncode)
        self.last_psb_monitor_end = datetime.now().isoformat()

        if 0 == launch_psb_monitor_return.returncode:
            self.load_psb_monitor_json()

        return launch_psb_monitor_return.returncode

    def sleep_after_report(self):
        if self.report_sleep_seconds < 60 * VUStructFlask_config['report_interval_initial_minutes']:
            app.logger.critical('Sleep seconds was not properly initialized')
            self.report_sleep_seconds = 60 * VUStructFlask_config['report_interval_initial_minutes']
        time.sleep(self.report_sleep_seconds) # Wait a few minutes before trying again 

    def increase_report_sleep_seconds(self):
        self.report_sleep_seconds = min(
            self.report_sleep_seconds * VUStructFlask_config['report_interval_stretch_factor'],
            VUStructFlask_config['report_interval_max_minutes'] * 60.0)

    @property
    def psb_rep_refresh_interval_seconds(self):
        return 5*60 # Need to think about this - annoying when it refreshes...

    def psb_rep(self, last_flag:bool=False, seconds_remaining:int=4*24*3600, refresh_interval_seconds=30) -> int:
        # Once we ask for a regenerated website it is critically important that that process complete.
        # In the event that the called psb_rep.py program fails, then we must create a web page with that information

        # !!!! remember to (generally) call psb_monitor() before calling this psb_rep

        # When we run psb_rep, we want to --transform out the uuid from the filenames
        # And we want the external_user filename prefixes to be removed as well.
        psb_rep_command_line = ['psb_rep.py','--tar_only', '--strip_uuid', self.case_uuid, '--web_case_id', self.case_id]
        if not last_flag:
            psb_rep_command_line.append('--embed_refresh')
            psb_rep_command_line.extend(['--refresh_interval_seconds', "%d" % refresh_interval_seconds])
            if seconds_remaining > 0:
                psb_rep_command_line.extend(['--seconds_remaining','%d' % seconds_remaining])

        app.logger.info("%s psb_rep: Starting %s in directory:\n%s", 
                        threading.get_ident(),
                        ' '.join(psb_rep_command_line),
                        self.working_directory)

        self.last_psb_rep_start = datetime.now().isoformat()
        launch_psb_rep_return = \
            subprocess.run(psb_rep_command_line, shell=False,
                       cwd = self.working_directory,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        app.logger.info("psb_rep: Returned %d from %s", launch_psb_rep_return.returncode, ' '.join(psb_rep_command_line))

        # Did psb_rep fail utterly?
        if launch_psb_rep_return.returncode != 0:
            # We now create a fail html page
            # This should be incredibly rare - a massive problem if psb_rep.py cannot deliver a webpage

            self.create_failure_website(body_html=\
f"""<H1>VUStruct failure while running the report generator</H1>
<p>This error should never happen.  Please let {VUStructFlask_config['webmaster']} know about this</p>
<A href=psb_rep.log>VUStruct Detailed Report Failure Log</A>""",
        				website_filelist = ['index.html', 'psb_rep.log'])

            app.logger.warn("psb_rep failure %d log is %s", launch_psb_rep_return.returncode, os.path.join(os.getcwd(),"psb_rep.log"))
        self.last_psb_rep_end = datetime.now().isoformat()

        GlobalActiveCases.add_website_refresh_needed(self.case_uuid)
        return launch_psb_rep_return.returncode

def monitor_report_loop(flask_case: CaseManager, relaunch:bool = False):
    """
    Primarily for development and testing, this allows a total relaunch, or start from monitor only
    of a web run.  Avoids needing to endlessly start new cases.
    """

    if relaunch:
        # Now we need to do launch (create the .slurm files as first step)
        # and the same code runs the generated ./launch...py outside the container
    
        flask_case.last_module_launched = 'launch'
        flask_case.last_module_returncode = flask_case.psb_launch()
        flask_case.last_module_error_message = ""
        if flask_case.last_module_returncode != 0:
            flask_case.last_module_error_message = "LATER Add error message for psb_launch please"
        # Even as these newly launched jobs are PENDING, go ahead and get that info
        global_squeue_df_refresh()

    psb_monitor_retcd = flask_case.psb_monitor()

    if psb_monitor_retcd != 0 :
        app.logger.error("psb_monitor.py for %s failed with %s", 
            flask_case.case_uuid, psb_monitor_retcd)


    # Error or not with psb_launch ...
    # We need to update the report page
    launch_psb_rep_retcd = flask_case.psb_rep(refresh_interval_seconds=flask_case.psb_rep_refresh_interval_seconds)


    # For 4 days loop, creating a new report every ?30? minutes
    # as the cluster churns away.  This is a separate functio to support
    # entry into _just_ this function (restart report creation)
    flask_case.timeout_time = datetime.now() + timedelta(days=float(VUStructFlask_config['case_timeout_days']))

    last_psb_rep = False

    last_psb_retcd  = 0

    while not last_psb_rep:
        if flask_case.last_module_returncode == 0:
            if flask_case.psb_monitor_results_dict:
                flask_case.psb_load_jobids_from_psb_monitor_results()
            else:
                flask_case.psb_load_jobids_from_workstatus_file()

        time_remaining = (flask_case.timeout_time - datetime.now()).total_seconds()
        if time_remaining <= 0:
            # Time's up - this is the last report run
            time_remaining = 0
            last_psb_rep = True

        psb_monitor_retcd = flask_case.psb_monitor()

        if psb_monitor_retcd != 0 :
            app.logger.error("psb_monitor.py for %s failed with %s", 
                flask_case.case_uuid, psb_monitor_retcd)

        # We will build a new report if psb_monitor was successful OR
        # if this is the very last psb_rep run
        if (psb_monitor_retcd == 0) or last_psb_rep :
            psb_rep_retcd = flask_case.psb_rep(
                last_flag=last_psb_rep,
                seconds_remaining=time_remaining,
                refresh_interval_seconds=flask_case.psb_rep_refresh_interval_seconds)
  
            if psb_rep_retcd == 0:
                GlobalActiveCases.add_website_refresh_needed(flask_case.case_uuid)
            else:
                app.logger.error("psb_rep for %s failed with %s", 
                    flask_case.case_uuid, psb_rep_retcd)

            # This delay loop here is important to avoid hammering the report generator
            flask_case.sleep_after_report()
            flask_case.increase_report_sleep_seconds()

    return last_psb_retcd

def squeue_flask_thread() -> None:
    """
    Every 2 minutes, we need to refresh the squeue status of all the jobs
    """
    
def scancel_flask_case_thread(flask_case: CaseManager) -> subprocess.CompletedProcess:
    """
    For every job# and array_id assigned to the case, scancel it
    """
    app.logger.info("Performing scancel on all of these running job ids:\n%s", str(flask_case.launched_jobs))

    job_id_set = set()
    scancel_cmdline = ['scancel'] 
    # Append job_ids to build complete command line to scancel
    # We do not specify array_ids here - because we want to kill everything 
    for launched_job in flask_case.launched_jobs:
        unique_key = launched_job[0]
        job_id =     launched_job[1]
        array_id =   launched_job[2]
        if job_id not in job_id_set:
            scancel_cmdline.append(str(job_id))
            job_id_set.add(job_id)

    scancel_return = \
        subprocess.run(scancel_cmdline, shell=False, encoding='UTF-8',
                   timeout=2*60,
                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    app.logger.info("scancel finshed with exit code %d", scancel_return.returncode)
    return scancel_return
        
def launch_flask_case_thread(flask_case: CaseManager, initial_launch: bool = True) -> subprocess.CompletedProcess:
    """
    This function is launched via the thread API.  "in the background".
    It runs the various components of the entire pipeline, and returns in case of module failure
    or after all the days that the pipeline is monitored

    Along the way, the website files will be periodically regenerated and the "needs web update" flag is set
 
    Importantly, this code is called for initial launches and vustruct_flask resumptions (in case of a restart of this service)
    Hence, the initial_launch flag, which is false if it is a resumption
    """

    GlobalActiveCases.add_flask_case(flask_case) # Likely a repeat - but no worries
    GlobalActiveCases.save_active_cases()

    if initial_launch:
        try:
            os.makedirs(flask_case.working_directory, exist_ok=True)
        except Exception as e:
            app.logger.exception(f"Exception during os.makedirs({flask_case.working_directory}): {e}")
            raise e
        flask_case.save_as_json()

    _fetch_error_msg = None

    # If we are relaunching an existing run, then it is impossible that
    # preprocessing has not been already performed
    if (initial_launch) or (not flask_case.last_module_launched):
        flask_case.last_module_returncode = 1
        # With the VUStruct format, we have skipped preprocess
        # And we just need to copy in our _vustruct.csv
        if flask_case.data_format == 'VUStruct CSV':
            flask_case.last_module_launched = 'plan'
            _fetch_error_msg = flask_case.fetch_input_file('csv')
            if not _fetch_error_msg:
                flask_case.last_module_returncode = 1
        else:
            # If we need to preprocess (parse) an input to create a missense file
            # Jump on that first and set last_module_launched to 'preprocess'
            # If no preprocessing is required by the user, then we flow through 
            # to psb_plan and process the missense.csv files...
            flask_case.last_module_returncode = 0
            if flask_case.data_format == 'Vanderbilt UDN Case Spreadsheet':
                flask_case.last_module_launched = 'preprocess'
                _fetch_error_msg = flask_case.fetch_input_file('xlsx')
                if not _fetch_error_msg:
                    flask_case.last_module_returncode = flask_case.parse_udn_report()
            elif flask_case.data_format == 'VCF GRCh38':
                flask_case.last_module_launched = 'preprocess'
                _fetch_error_msg = flask_case.fetch_input_file('vcf')
                if not _fetch_error_msg:
                    flask_case.last_module_returncode = flask_case.vcf2vustruct()
            elif flask_case.data_format == 'VCF GRCh37 (runs liftover)':
                flask_case.last_module_launched = 'preprocess'
                _fetch_error_msg = flask_case.fetch_input_file('vcf')
                if not _fetch_error_msg:
                    flask_case.last_module_returncode = flask_case.vcf2vustruct_with_liftover()
            elif flask_case.data_format == 'WUSTL Case Spreadsheet':
                flask_case.last_module_launched = 'preprocess'
                _fetch_error_msg = flask_case.fetch_input_file('xlsx')
                if not _fetch_error_msg:
                    flask_case.last_module_returncode = flask_case.parse_wustl()
            flask_case.save_as_json()

    if not _fetch_error_msg and flask_case.last_module_returncode == 0:
        # Read the resulting pipeline input vustruct.csv file into our DataFrame
        flask_case.load_vustruct_csv()

    # Whatever we might have done with genomic starting coordinates... proceed!
    if flask_case.last_module_launched == 'preprocess':
        flask_case.last_module_error_message = ""
        preprocess_returncode = flask_case.last_module_returncode
        if preprocess_returncode != 0:
            flask_case.last_module_error_message = "LATER Add error message for parse_udn_report please"
        # After we complete the parse step, we definitely need to run the psb_rep program
        # so that the OTHER flask program will install the growing website if it comes to that.
        # HOWEVER, the psb_rep program must note if the parse process has terminated with bad exit code
        # and take that into consideration.  In the event of error return code, we return from this launch_thread
        # and we populate the report with last_flag set because there are not going to be any more psb_rep.py calls
        launch_psb_rep_retcd = flask_case.psb_rep(last_flag = (flask_case.last_module_returncode != 0),refresh_interval_seconds=10)

        if launch_psb_rep_retcd != 0:
            flask_case.last_report_error = "psb_rep.py for %s failed with %s" % \
                                             (flask_case.case_uuid, launch_psb_rep_retcd)

        # In the event the preprocessor failed, then the case has finished.
        # This thread must exit and the global dictionary of active cases cleared
        if preprocess_returncode != 0:
            GlobalActiveCases.mark_final_website_update(flask_case.case_uuid)
            return None
    

    # Whether we parsed genomic coordinates or not, the NEXT thing
    # to do is run psb_plan.
    flask_case.last_module_launched = 'plan'


    # Whether directly, or by a preprocess step - we now psb_plan
    flask_case.last_module_returncode = flask_case.psb_plan()
    flask_case.last_module_error_message = ""
    if flask_case.last_module_returncode != 0:
        flask_case.last_module_error_message = "psb_plan failed.  Review log files"
        return None


    flask_case.save_as_json()

    # Error or not with psb_plan...
    # Now that we did psb_plan, we need to report
    # So that the user sees psb_plan is complete, can review logs, etc.
    launch_psb_rep_retcd = flask_case.psb_rep(refresh_interval_seconds=20)

    if launch_psb_rep_retcd != 0:
        flask_case.last_report_error = "psb_rep.py for %s failed with %s" % \
                                         (flask_case.case_uuid, launch_psb_rep_retcd)


    # Plan can fail for lots of reasons.  In that case, user will see that the case is over
    # and we halt flask case management right here
    if flask_case.last_module_returncode != 0:
        return None

    # Now we need to do launch (create the .slurm files as first step)
    # and the same code runs the generated ./launch...py outside the container

    flask_case.last_module_launched = 'launch'
    flask_case.last_module_returncode = flask_case.psb_launch()
    flask_case.last_module_error_message = ""
    if flask_case.last_module_returncode != 0:
        flask_case.last_module_error_message = "LATER Add error message for psb_launch please"

    # Error or not with psb_launch ...
    # Now that we did psb_plan, we need to report
    launch_psb_rep_retcd = flask_case.psb_rep(refresh_interval_seconds=flask_case.psb_rep_refresh_interval_seconds)

    if launch_psb_rep_retcd != 0:
        flask_case.last_report_error = "psb_rep.py for %s failed with %s" % \
                                         (flask_case.case_uuid, launch_psb_rep_retcd)
    if flask_case.last_module_returncode != 0:
        return None

    monitor_report_loop(flask_case)

    return 0 # For now - not sure what to do after 4 days psb_rep_retcd

# def launch_flask_case(uuid_str: str):
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
    """
    Implement an authentication layer for endpoints per discussions with ACCRE
    """
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

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
#  Restart after a halted session - super important    #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
GlobalActiveCases.reload_from_save_file()




# Simply getting back a response can tell caller
# That we are running
@app.route('/alive', methods = ['POST'])
@require_authn
def alive():
    return "{}"

# def shutdown_server():
# Does not work since several years ago - but idea is not terrible
# Which is to try to use the shutdown of the werkzeug environment
#     func = request.environ.get('werkzeug.server.shutdown')
#     if func is None:
#         raise RuntimeError('Not running with the Werkzeug Server')
#     func()

@app.route('/halt', methods=['POST','GET'])
@require_authn
def halt():
    """
    Exit the server
    """
    kill_vustruct_flask("The /halt endpoint was posted to")

    return "Caller won't see this because this process is now dead"

@app.route('/case_status', methods=['POST','GET'])
@require_authn
def case_status():
    return {'CaseStatus': "Keep goin on vu!"}
 
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

    app.logger.info ("In /get_uuid: Received data: %s json: %s" % (request.data, request.json))
    # Launch psb_plan.py

    uuid_str = str(uuid.uuid4())

    app.logger.info("POST to /getuuid returning %s", uuid_str)

    return jsonify({"case_uuid": uuid_str})  # 'Hello ' #  + rq.get('name', 'No name'))

@app.route('/launch_vustruct', methods = ['POST', 'GET'])
@require_authn
def launch_vustruct():
    """
    There are a variety of ways to launch a new VUStruct case - 
    Some start frmo the missense.csv.  Others from a spreadsheet format.
    
    An example curl command would be:
    curl -X POST http://localhost:5000/launch_vustruct \ 
        -H 'Content-Type: application/json' \
        -d'{"data_format": "Vanderbilt UDN Case Spreadsheet",
            "case_uuid": "12435",
            "case_id": "test",
            "upload_file_URI": "file:////dors/capra_lab/users/mothcw/VUStruct/fakecase0/fakecase0.xlsx"}'
    """

    app.logger.debug ("/launch_vustruct: request.json=\n%s", json.dumps(request.json,indent=3))

    for required_key in ['case_id', 'case_uuid', 'data_format', 'upload_file_URI']:
        if not required_key in request.json:
            error_msg = "Request JSON is missing required key: " + required_key
            app.logger.critical(error_msg)
            return {'error': error_msg}

    flask_case = CaseManager(
        case_id=request.json['case_id'],
        case_uuid=request.json['case_uuid'])
    flask_case.data_format = request.json['data_format']

    known_data_formats = [
        'VUStruct CSV',
        'VCF GRCh38',
        'VCF GRCh37 (runs liftover)',
        'Vanderbilt UDN Case Spreadsheet',
        'WUSTL Case Spreadsheet']

    if  flask_case.data_format not in known_data_formats:
        error_msg = "Data format %s is not recognized.  Must be one of: %s" % (
            flask_case.data_format,','.join(known_data_formats) )
        app.logger.critical(error_msg)
        return {'error': error_msg}


    # NO - vustruct_csv is no longer data typed in manually
    # It's an uploaded file like anythign else
    # flask_case.vustruct_csv = request.json['missense_csv']
    # ^ If we do something else with wordpress later, then this changes

    if 'upload_file_URI' in request.json:
        flask_case.upload_file_URI = request.json['upload_file_URI']

    # This so far has been very quick.  Now run the background thread to run the entire pipeline
    # Return to caller immediately of course.
    # We pass in "True" to tell the thread it is an initial launch of the pipeline
    
    parse_plan_launch_thread = threading.Thread(target=launch_flask_case_thread, args=(flask_case, True))
    parse_plan_launch_thread.start()

    # Send back the entire json that was send to us AnD add a launch_status to the gemisch
    # json_to_return_after_launch = request.json | {"launch_status": "VUStruct now running on accre cluster"}
    flask_case_json = flask_case.to_dict()
    app.logger.debug("Sending back this json %s", pprint.pformat(flask_case_json) )

    return flask_case_json

@app.route('/dev_monitor_report_loop', methods = ['POST'])
@require_authn
def dev_monitor_report_loop():
    """
    For development purposes, it can be very helpful to loop on psb_monitor.py and psb_rep.py after 
    Therefore, there is this target to permit a "restart" of course

    If jobs already running, then do NOT include a relaunch element
    
    An example curl command would be:
    curl -X POST http://localhost:5000/report_vustruct \ 
        -H 'Content-Type: application/json' \
        -d'{"case_uuid": "12435",
            "case_id": "test"

    You can also add relaunch: "True" or anything
    """

    app.logger.debug ("/dev_report_vustruct: request.json=%s", str(request.json))

    for required_key in ['case_id', 'case_uuid']:
        if not required_key in request.json:
            return "Request JSON is missing required key: " + required_key

    flask_case = CaseManager(
        case_id=request.json['case_id'],
        case_uuid=request.json['case_uuid'])

    relaunch = ('relaunch' in request.json)

    # If we cannot read this file, then we are clearly having a major problem
    # and the user needs to understand that we cannot proceed
    try:
        flask_case.load_vustruct_csv()
    except Exception as ex:
        return jsonify({"Error": str(ex)})

    # Even though we are going to start assuming that launch went great
    # We need to load the last state of things
    try:
        flask_case.load_from_json()
    except Exception as ex:
        return jsonify({"Error": str(ex)})

    # Add back this case to the global_cases_dict, in case it has been lost
    GlobalActiveCases.add_flask_case(flask_case)

    flask_case.last_module_launched = 'launch'
    flask_case.last_module_returncode = 0
    flask_case.last_module_error_message = ""

    # This so far has been very quick.  Now run the background thread to run the entire pipeline
    # Return to caller immediately of course.

    dev_report_vustruct_thread = threading.Thread(target=monitor_report_loop, args=(flask_case,relaunch))
    dev_report_vustruct_thread.start()

    # Send back the entire json that was send to us AnD add a launch_status to the gemisch
    # json_to_return_after_launch = request.json | {"launch_status": "VUStruct now running on accre cluster"}
    app.logger.debug("Sending back this json %s", flask_case.to_dict() )

    return flask_case.to_dict()


@app.route('/cases_needing_website_refresh', methods = ['POST', 'GET'])
@require_authn
def cases_needing_website_refresh():
    """
    Every 10 to 30 seconds, vustruct_webupdate.py calls this endpoint to learn of the
    all the cases where psb_rep.py has generated a new report, and thus need to be 
    copied into the filesystem of the webserver

    The return value is a list of dictionaries with case_id/case_uuid/working_directory
    """

    return jsonify(GlobalActiveCases.cases_needing_website_refresh(clear=True))
    

@app.route('/peek_cases_needing_refresh', methods = ['POST', 'GET'])
@require_authn
def peek_cases_needing_refresh():
    """
    Return cases needing refresh, without clearing
    This endpoint is not used currently
    """

    return jsonify(GlobalActiveCases.cases_needing_website_refresh(clear=False))


@app.route('/add_case', methods = ['POST'])
@require_authn
def add_case():
    """
    An example curl command would be:
    curl -X POST http://localhost:5000/add_case \ 
        -H 'Content-Type: application/json' \
        -d'{"data_format": "Vanderbilt UDN Case Spreadsheet", \
            "case_uuid": "b60003f4-08bb-4821-83c5-bf6ddaa7e5fa",
            "case_id": "SAMPLEJUNE2024_3PM"}'
    """

    new_case_info = request.json
    app.logger.info("add_case request json = %s" % new_case_info)
    _flask_case = CaseManager(new_case_info['case_id'], new_case_info['case_uuid'])
    GlobalActiveCases.add_flask_case(_flask_case)
    # Relaunch the thread - with the initial flag set to False so that it knows to look for previous
    # parse results.  Assumes that the raw user suppled file was already placed where it needs to be
    parse_plan_launch_thread = threading.Thread(target=launch_flask_case_thread, args=(_flask_case,False))
    parse_plan_launch_thread.start()

    return jsonify(vars(_flask_case))

@app.route('/scancel_case', methods = ['POST'])
@require_authn
def scancel_case():
    """
    An example curl command would be:
    curl -X POST http://localhost:5000/scancel_case \ 
        -H 'Content-Type: application/json' \
        -d'{"case_uuid": "b60003f4-08bb-4821-83c5-bf6ddaa7e5fa"}'
    """

    scancel_case_info = request.json
    app.logger.info("scancel_case request json = %s", scancel_case_info)


    _flask_case = GlobalActiveCases.get_flask_case_or_none(scancel_case_info['case_uuid'])
    _scancel_return = {'flask_error': 'Code failed to initialize'}
    if _flask_case is not None:
        slurm_job_count = len(_flask_case.launched_jobs)
        if not slurm_job_count:
            _scancel_return = {
                'case_uuid': scancel_case_info['case_uuid'],
                'flask_error': "No launched slurm jobs to scancel"
            }
        else:
            scancel_thread = threading.Thread(target=scancel_flask_case_thread, args=(_flask_case,))
            scancel_thread.start()
            _scancel_return =  vars(_flask_case) | {'flask_success': 'scanceled %d slurm jobs' % slurm_job_count}
    else: # Return error!
       _scancel_return = {
            'case_uuid': scancel_case_info['case_uuid'],
            'flask_error': "Not found in GlobalActiveCases dictionary"
            }
              

    return jsonify(_scancel_return)

@app.route('/squeue_monitor', methods = ['POST'])
@require_authn
def squeue_monitor():
    """
    vustruct_flask.py is constantly harving squeue details
    and also json returned from psb_monitor.py after each run.

    By combinging these details, we create an end user UI table
    to see progress being made by all the jobs

    An example curl command would be:
    curl -X POST http://localhost:5000/squeue_monitor \ 
        -H 'Content-Type: application/json' \
        -d'{"case_uuid": "b60003f4-08bb-4821-83c5-bf6ddaa7e5fa"}'
    """

    app.logger.info("squeue_monitor request json = %s", request.json)

    _flask_case = GlobalActiveCases.get_flask_case_or_none(request.json['case_uuid'])
    _sm_return = {}

    _sm_dict = {}
    if _flask_case is None:
        _sm_return = {
            'case_uuid': request.json['case_uuid'],
            'flask_error': "Case  not found"}
        return _sm_return

    if len(_flask_case.vustruct_csv_df) == 0:        
        _sm_return = {
            'case_uuid': request.json['case_uuid'],
            'flask_error': "Case somehow lacks a vustruct.csv file.  Cannot continue"}
        return _sm_return

    # For each variant row in the vustruct file
    for index, row in _flask_case.vustruct_csv_df.iterrows():
        # Now load all the job information for the variant, as returned by
        # psb_monitor
        gene_refseq_mutant = "%s_%s_%s" % (row['gene'], row['refseq'], row['mutation'])
        gfm_jobs_info = _flask_case.psb_monitor_results_dict.get(gene_refseq_mutant, None)
        if gfm_jobs_info:
            _sm_return[gene_refseq_mutant]  = {}
            for unique_id, psb_monitor_single_job_info in gfm_jobs_info.items():
                job_id_int = int(psb_monitor_single_job_info['jobid'])
                array_id_int = psb_monitor_single_job_info['arrayid']
                if array_id_int is not None and len(array_id_int):
                    job_key = "%s_%s" % (job_id_int, int(float(array_id_int)))
                else:
                    job_key = str(job_id_int)

                try:
                    # Combine what psb_monitor.py returned about the jobs
                    # last time we ran it
                    # With the dynamic squeue data we got from $squeue command
                    matching_squeue_entry = global_squeue_df.loc[[job_key]].iloc[0].to_dict()
                    # We need to convert the timetypes to strings in the UI
                    # We chop out year and everything after minutes with [5:16]
                    matching_squeue_entry['start_time'] = \
                        str(matching_squeue_entry['start_time'])[5:16]
                    matching_squeue_entry['end_time'] = \
                        str(matching_squeue_entry['start_time'])[5:16]

                    # For the time_limit, divide by one minute of time to get integer
                    # Then convert to hours float by dividing by 60.0
                    matching_squeue_entry['time_limit'] = str(float(
                       matching_squeue_entry['time_limit']  / 
                       np.timedelta64(1,'m'))/60.0) + 'h'
                       
                    _sm_return[gene_refseq_mutant][unique_id] = \
                             psb_monitor_single_job_info | \
                             matching_squeue_entry | \
                             {'job_key': job_key}
                except KeyError as ex:
                    # So, there is no slurm information about this jobid.  That's A-OK!
                    pass
                
    return _sm_return
