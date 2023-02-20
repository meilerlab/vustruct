#!/usr/bin/env python3
#
# vustruct_webdate
#
# A second flask program which
# - runs local alongside the wordpress server and
# - implements RESTAPI calls which:
#   1. record job identifiers which are in process
#   2. rebuild local websites in file system from emerging vustruct .tar.gz files
#
# For now, set the CASE...BASE directories in the lines below
CASE_URL_BASE="https://staging.meilerlab.org/vustruct"
CASE_FILESYSTEM_BASE="/var/www/staging.meilerlab.org/vustruct"
VUSTRUCT_FLASK="https://api.vgi01.accre.vanderbilt.edu"

import logging

from flask import Flask
jobs_dict = {}
job_uuids_needing_website_refresh = set()

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


# UDN="/dors/capra_lab/users/mothcw/UDNtests/"
UDN = os.getenv("UDN")
if not UDN:
    errormsg = "UDN environment variable must be defined before launch of %s" % UDN
    sys.exit(errormsg)



cmdline_parser = psb_config.create_default_argument_parser(__doc__, os.path.dirname(os.path.dirname(__file__)))

args, remaining_argv = cmdline_parser.parse_known_args()

# sh = logging.StreamHandler()
logging.basicConfig(format='%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s',
datefmt='%d-%m:%H:%M:%S', )
LOGGER = logging.getLogger()
# LOGGER.addHandler(sh)

# Now that we've added streamHandler, basicConfig will not add another handler (important!)
LOG_FORMAT_STRING = '%(asctime)s %(levelname)-4s [%(filename)16s:%(lineno)d] %(message)s'
DATE_FORMAT_STRING = '%H:%M:%S'
log_formatter = logging.Formatter(LOG_FORMAT_STRING, DATE_FORMAT_STRING)

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

# Create a mapping of the unique uuids for running jobs to user and job names
# so that we know where to go fish out emerging websites.`
jobs_dict = {}

@app.route('/health_check', methods=['POST'])
def health_check():
    health_good_dict = {'health': 'good'}
    LOGGER.info("/health_check.. returning  %s", health_good_dict)
    return health_good_dict

@app.route('/add_uuid', methods=['POST'])
def add_uuid():
    LOGGER.info("/add_uuid: %s" % request.json)

    case_url = "None"

    if not ('job_uuid' in request.json and 'case_id' in request.json):
        LOGGER.info("/add_uuid: missing job_uuid or case_id")
        return {}

    case_url = os.path.join(CASE_URL_BASE,
                            request.json['job_uuid'],
                            request.json['case_id'])

    webpage_home = os.path.join(CASE_FILESYSTEM_BASE,
                                request.json['job_uuid'],
                                request.json['case_id'])
    os.makedirs(webpage_home, exist_ok=True)

    index_html = os.path.join(webpage_home, "index.html")
    with open(index_html, "w") as initial_webpage_f:
        initial_webpage_f.write("""\
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
  "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html lang="en-us" xmlns="http://www.w3.org/1999/xhtml" >
<head>
        <meta charset="utf-8">
        <title>Personal Structural Biology case %s</title>

      <script type = "text/JavaScript">
            var refreshCount = 0;
            function AutoRefresh( t ) {
               setTimeout("{refreshCount++; location.reload(true);}", t);
            }
      </script>
</head>
<body>
<img src="../../html/PSBgraphicHeader.png" alt="VUStruct Graphic"/>
</body>
<p>
<hr>
<br>
Webpage has been refreshed
<script>document.write(refreshCount)</script>
times
        <H1>Personal Structural Biology case %s</H1>
      <script type = "text/JavaScript">
var date = new Date();
        var current_time = date.getHours()+":"+date.getMinutes()+":"+ date.getSeconds();
document.write("Time: " + current_time);
</script>
<button>
    <a href="javascript:{refreshCount++; location.reload(true)}"> Click Refresh the Page</a>
</button>

 </head>
   <body onload = "JavaScript:AutoRefresh(10000);">
      <p>Your Page will refresh every 10 seconds!</p>
   </body>
</html>
""" % ( request.json['case_id'], request.json['case_id'] ) )


    # Launch psb_plan.py

    add_uuid_return = {'webpage_url': case_url}
    LOGGER.info("/add_uuid: returning %s", add_uuid_return)
    return add_uuid_return

def xfer_to_web_thread(job_needing_refresh) -> subprocess.CompletedProcess:
    LOGGER.info("Updating website for case: %s  job_uuid: %s", job_needing_refresh['case_id'], job_needing_refresh['job_uuid'])

    # LOGGER.info("Copying to live web: %s", job_needing_refresh)
    job_uuids_needing_website_refresh.remove(job_needing_refresh['job_uuid'])

    internal_case_name = "external_user_%s_%s" % (
        job_needing_refresh['case_id'],
        job_needing_refresh['job_uuid'])
    tar_filename = os.path.join(
                   job_needing_refresh['working_directory'],
                   "%s.tar.gz" % internal_case_name)

    # open file
    with tarfile.open(tar_filename) as tar_f:
        tar_f.extractall(CASE_FILESYSTEM_BASE)

    webpage_home = os.path.join(CASE_FILESYSTEM_BASE,
                                  job_needing_refresh['job_uuid'],
                                  job_needing_refresh['case_id'])
    os.makedirs(webpage_home, exist_ok=True)

    # Now recursively move all the files in the tar x directory to the correct directory
    # LATER - make this a mv.  For now brute copy
    cp_source = os.path.join(CASE_FILESYSTEM_BASE, internal_case_name )
    cp_dest = os.path.join(CASE_FILESYSTEM_BASE, job_needing_refresh['job_uuid'], job_needing_refresh['case_id']);
    shutil.copytree(src=cp_source,dst=cp_dest,dirs_exist_ok=True)
    # for root, dirs, files in os.walk(os.path.join('/var/www/html/vustruct/')

    # Clean up the index.html to just show the case name - and NOT
    # the uuid and so forth
    index_html_filename = os.path.join(cp_dest,"index.html")
    LOGGER.info("Cleaning %s", index_html_filename)
    index_html_content = ""
    with open(index_html_filename, 'r') as f:
        index_html_content = f.read();
    index_html_content  = index_html_content.replace(internal_case_name, job_needing_refresh['case_id'])
    with open(index_html_filename,'w') as f:
        f.write(index_html_content);



    # Get rid of any old file that is there



    return None

def refresh_case_websites():
    # Go fetch the jobs that our vustruct_flask application is managing, and which have just
    # run psb_rep.py
    jobs_needing_refresh = requests.get(
        os.path.join(VUSTRUCT_FLASK,'jobs_needing_refresh'),
        headers={'Authorization': 'vustruct_password'}
        )
    if jobs_needing_refresh.status_code != 200:
        LOGGER.warning("vustruct_flask.py not responding.  Check environment!")
    else:
        jobs_needing_refresh = jobs_needing_refresh.json()
        LOGGER.info("%d jobs need website transfers" % len(jobs_needing_refresh))

        for job_needing_refresh in jobs_needing_refresh:
            LOGGER.info("%s", job_needing_refresh)
            job_uuids_needing_website_refresh.add(job_needing_refresh['job_uuid'])
            psb_plan_thread = threading.Thread(target=xfer_to_web_thread, args=(job_needing_refresh,))
            psb_plan_thread.start()
    # sys.exit(0)

refresh_case_websites()
scheduler = BackgroundScheduler()
scheduler.add_job(func=refresh_case_websites, trigger="interval", seconds=10)
scheduler.start()

# Shut down the scheduler when exiting the app
atexit.register(lambda: scheduler.shutdown())
