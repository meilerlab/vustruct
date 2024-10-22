#!/usr/bin/env python3
#
# vustruct_webdate
#
# A second flask program which
# - runs local alongside the wordpress server and
# - implements RESTAPI calls which:
#   1. record case identifiers which are in process
#   2. rebuild local websites in file system from emerging vustruct .tar.gz files
#
# For now, set the CASE...BASE directories in the lines below
from __future__ import annotations

import uuid
import threading

CASE_URL_BASE="https://staging.meilerlab.org/vustruct"
CASE_FILESYSTEM_BASE="/var/www/staging.meilerlab.org/vustruct"
VUSTRUCT_FLASK="https://api.vgi01.accre.vanderbilt.edu"

import logging
from logging.config import dictConfig


from flask import Flask
case_dict = {}
case_uuids_needing_website_refresh = set()

import atexit
from apscheduler.schedulers.background import BackgroundScheduler

import tarfile

from flask import request
from flask import jsonify
import uuid
import threading
from datetime import datetime
import time
import json
import subprocess
import os
import shutil
import sys

from slurm import slurm_submit
from bsub import bsub_submit

# requests is used to read the excel or missense filess uploaded to wordpress.
import requests
from requests_file import FileAdapter
# import time
from datetime import datetime
from logging.handlers import RotatingFileHandler
from typing import Dict, Any, Tuple, List
from typing import TextIO
# Inspect to recover bsub_submit() and slurm_submit() source code for integration
# into generated final file
import inspect

from psb_shared import psb_config
from psb_shared import psb_perms

os.makedirs('log', exist_ok=True)

# Setup the logging.  By calling dictConfig, the flask app.logger
# will automatically pick up these settings
# This is the "way to go" with a flask application
LOG_FILENAME='log/vustruct_webupdate.log'
logging.config.dictConfig({
    'version': 1,
    'formatters': {'default': {
        'format': '%(asctime)s %(levelname)-4.4s [%(filename)16s:%(lineno)d] %(message)s',
        'datefmt': '%H:%M:%S'
    }},
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
        'handlers': ['wsgi']
    }
})

app = Flask(__name__)

app.logger.info("Flask app created with Flask(\'%s\')" % __name__)
assert app.logger.getEffectiveLevel() == logging.DEBUG
app.logger.setLevel(logging.DEBUG)
# 
app.logger.info("Created Flask app")
assert app.logger.getEffectiveLevel() == logging.DEBUG


# UDN="/dors/capra_lab/users/mothcw/VUStruct/"
UDN = os.getenv("UDN")
if not UDN:
    errormsg = "UDN environment variable must be defined before launch of %s" % UDN
    sys.exit(errormsg)



cmdline_parser = psb_config.create_default_argument_parser(__doc__, os.path.dirname(os.path.dirname(__file__)))

args, remaining_argv = cmdline_parser.parse_known_args()

# sh = logging.StreamHandler()
# LOGGER = logging.getLogger()
# app.logger.addHandler(sh)

# Now that we've added streamHandler, basicConfig will not add another handler (important!)
# LOG_FORMAT_STRING = '%(asctime)s %(levelname)-4s [%(filename)16s:%(lineno)d] %(message)s'
# DATE_FORMAT_STRING = '%H:%M:%S'
# logging.basicConfig(format='%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s',
#     datefmt='%d-%m:%H:%M:%S', )

os.makedirs('log', exist_ok=True)

#
# log_filename = os.path.join('log', 'vustruct_webupdate.log')
#
# app.logger.info("Additionally logging to file %s" % log_filename)
# needRoll = os.path.isfile(log_filename)
#
# local_fh = RotatingFileHandler(log_filename, backupCount=7)
# local_fh.setLevel(logging.DEBUG)
# app.logger.addHandler(local_fh)

# if needRoll:
#   local_fh.doRollover()

# print("Level is %s " % os.getenv("FLASK_ENV"))
if args.debug or (os.getenv("FLASK_ENV") == 'development'):
    print("Setting Debug level")
    app.logger.setLevel(logging.DEBUG)
elif args.verbose:
    app.logger.setLevel(logging.INFO)
else:
    app.logger.setLevel(logging.WARNING)
# app.logger.setFormatter(log_formatter)



app.logger.info("%s updating the main website %s", __file__, UDN)
app.logger.info("Initializing Flask(%s)", __name__)

# Create a mapping of the unique uuids for the running cases to user and job names
# so that we know where to go fish out emerging websites.`
case_dict = {}

@app.route('/health_check', methods=['POST'])
def health_check():
    health_good_dict = {'health': 'good'}
    app.logger.info("/health_check.. returning  %s", health_good_dict)
    return health_good_dict

# Return the latest slurm squeue details about a particular job
# This request is directly forwarded on to vustruct_flask
# and that response is returned to the PHP code calling this
# entry point
@app.route('/case_status', methods=['GET'])
def case_status():
    app.logger.info("IN CASE STATUS ")
    _title = request.args.get('title')
    _uuid = request.args.get('uuid')
    app.logger.info("HERE IN CASE_STATUS! %s %s",request.args.get('title'),request.args.get('uuid'))
    vustruct_flask_case_status = requests.get(
        os.path.join(VUSTRUCT_FLASK,'case_status'),
        params={'title': _title, 'uuid': _uuid},
        headers={'Authorization': 'vustruct_password'}
        )

    if vustruct_flask_case_status.status_code != 200:
        app.logger.warning("vustruct_flask.py not responding (%d).  Check environment!", \
            vustruct_flask_case_status.status_code)

    return vustruct_flask_case_status.json()

# Create a "placeholder" webpage, quickly, and return the https:// URL
# to the user so that after the [SUBMIT] the user can start to see their emerging 
# page
@app.route('/create_placeholder_webpage', methods=['POST'])
def create_placeholder_webpage():
    app.logger.info("/create_placeholder_webpage: %s\n" , json.dumps(request.json,indent=2))

    case_url = "None"

    if not ('case_uuid' in request.json and 'case_id' in request.json and 'data_format' in request.json):
        app.logger.info("/create_placeholder_webpage: missing case_uuid, case_id, or data_format")
        return {}

    case_url = os.path.join(CASE_URL_BASE,
                            request.json['case_uuid'],
                            request.json['case_id'])

    webpage_home = os.path.join(CASE_FILESYSTEM_BASE,
                                request.json['case_uuid'],
                                request.json['case_id'])
    os.makedirs(webpage_home, exist_ok=True)

    initial_task_information = "VUStruct processing is starting..."
    if ('pasted_vustruct_csv' in request.json and request.json['pasted_vustruct_csv']):
       vustruct_lines = len(request.json['pasted_vustruct_csv'].split('\n'))
       vustruct_variants = vustruct_lines - 1 # Subtract one because we don't count the header
       initial_task_information = "VUStruct is planning calculations on best available PDB and model structures for your "
       if vustruct_variants == 1:
           initial_task_information += "variant."
       else:
           initial_task_information += "%d variants." % vustruct_variants
    else:
       initial_task_information = "VUStruct is converting file: <b>%s</b> (%s) genomic coordinates to protein variants." % (
           os.path.basename(request.json['upload_file_URI']), request.json['data_format']);

    index_html = os.path.join(webpage_home, "index.html")
    with open(index_html, "w") as initial_webpage_f:
        local_refresh_variable_name = 'refresh_count_' + request.json['case_uuid']
        initial_webpage_f.write(f"""\
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
  "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html lang="en-us" xmlns="http://www.w3.org/1999/xhtml" >
<head>
        <meta charset="utf-8">

        <title>VUStruct case {request.json['case_id']}</title>

      <script type = "text/JavaScript">
            window.refreshCount = localStorage.getItem('{local_refresh_variable_name}');
            if (window.refreshCount)
                window.refreshCount = parseInt(window.refreshCount)
            else
                window.refreshCount = 0;
      </script>
"""

            + f"""
</head>
   <body onload = "JavaScript:AutoRefresh(10000);">
<img src="../../html/PSBgraphicHeader.png" alt="VUStruct Graphic"/>
<p>
<hr>
        <H1>VUStruct case {request.json['case_id']}</H1>

"""

+ "<p>" + initial_task_information+ "<p>" + 

"""<br>
      <p>This page will refresh every 10 seconds.</p>

<script type="text/JavaScript">
if (window.refreshCount > 0) {
    document.write('This page has been refreshed ');
    if (window.refreshCount == 1)
        document.write('once');
    else {
        document.write(window.refreshCount.toString())
        document.write(' times.');
      }
}
"""
+
   f""" 
function incrementRefreshCountAndReload() {{
    localStorage.setItem(
          '{local_refresh_variable_name}',
           (window.refreshCount + 1).toString());
    location.reload(true); 
    }};


function AutoRefresh( timeout_milliseconds ) {{
   setTimeout( incrementRefreshCountAndReload, timeout_milliseconds);
}}
</script>

<!--
 <button>
    <a href="javascript:{{incrementRefreshCountAndReload()}}"> Click Refresh the Page</a>
</button> 
-->

<script type="text/javascript" id="hourglass">
document.body.style.cursor = 'wait';
</script>


   </body>
</html>
""")


    # Launch psb_plan.py

    create_placeholder_webpage_return = {'webpage_url': case_url}
    app.logger.info("/create_placeholder_webpage: returning %s", create_placeholder_webpage_return)
    return create_placeholder_webpage_return

def xfer_to_web_thread(case_needing_refresh) -> subprocess.CompletedProcess:
    app.logger.info("Updating website for case: %s  case_uuid: %s", case_needing_refresh['case_id'], case_needing_refresh['case_uuid'])

    # app.logger.info("Copying to live web: %s", case_needing_refresh)
    case_uuids_needing_website_refresh.remove(case_needing_refresh['case_uuid'])

    internal_case_name = "external_user_%s_%s" % (
        case_needing_refresh['case_id'],
        case_needing_refresh['case_uuid'])
    tar_filename = os.path.join(
                   case_needing_refresh['working_directory'],
                   "%s.tar.gz" % internal_case_name)

    # open file
    app.logger.info("Extracting files in %s to %s", tar_filename, CASE_FILESYSTEM_BASE)
    with tarfile.open(tar_filename) as tar_f:
        tar_f.extractall(CASE_FILESYSTEM_BASE)

    webpage_home = os.path.join(CASE_FILESYSTEM_BASE,
                                  case_needing_refresh['case_uuid'],
                                  case_needing_refresh['case_id'])
    os.makedirs(webpage_home, exist_ok=True)

    # Now recursively move all the files in the tar x directory to the correct directory
    # LATER - make this a mv.  For now brute copy
    cp_source = os.path.join(CASE_FILESYSTEM_BASE, internal_case_name )
    cp_dest = os.path.join(CASE_FILESYSTEM_BASE, case_needing_refresh['case_uuid'], case_needing_refresh['case_id']);
    app.logger.info("shutil.copytree(src=%s,dst=%s,dirs_exist_ok=True", cp_source, cp_dest)
    shutil.copytree(src=cp_source,dst=cp_dest,dirs_exist_ok=True)
    # for root, dirs, files in os.walk(os.path.join('/var/www/html/vustruct/')

    # Clean up the index.html to just show the case name - and NOT
    # the uuid and so forth
    index_html_filename = os.path.join(cp_dest,"index.html")
    app.logger.info("Cleaning %s", index_html_filename)
    index_html_content = ""
    with open(index_html_filename, 'r') as f:
        index_html_content = f.read();
    index_html_content  = index_html_content.replace(internal_case_name, case_needing_refresh['case_id'])
    with open(index_html_filename,'w') as f:
        f.write(index_html_content);

    # Get rid of any old file that is there
    app.logger.info("Make sure that %s is totally removed now", cp_source)

    return None

def refresh_case_websites():
    # Go fetch the cases that our vustruct_flask application is managing, and which have just
    # run psb_rep.py
    # app.logger.info("In refresh_case_websites()")
    cases_needing_refresh = requests.get(
        os.path.join(VUSTRUCT_FLASK,'cases_needing_website_refresh'),
        headers={'Authorization': 'vustruct_password'}
        )
    if cases_needing_refresh.status_code != 200:
        app.logger.warning("vustruct_flask.py not responding (%d).  Check environment!", \
            cases_needing_refresh.status_code)
    else:
        cases_needing_refresh = cases_needing_refresh.json()
        # In most cases, there won't be any cases to refresh
        if len(cases_needing_refresh):
            app.logger.info("%d cases need website transfers" % len(cases_needing_refresh))

            # When there are cases to refresh, launch I/O blocking threads to do this "in the background"
            for case_needing_refresh in cases_needing_refresh:
                app.logger.info("Details of case being refreshed to website:  %s", json.dumps(case_needing_refresh, indent=3))
                case_uuids_needing_website_refresh.add(case_needing_refresh['case_uuid'])
                _thread = threading.Thread(target=xfer_to_web_thread, args=(case_needing_refresh,))
                _thread.start()
    # sys.exit(0)

# Turn off the super noisy messages from the scheduler every 15 seconds
logging.getLogger('apscheduler.executors.default').propagate = False
refresh_case_websites()
scheduler = BackgroundScheduler()
scheduler.add_job(func=refresh_case_websites, trigger="interval", seconds=15)
scheduler.start()

# Shut down the scheduler when exiting the app
atexit.register(lambda: scheduler.shutdown())
