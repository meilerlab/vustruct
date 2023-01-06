#!/usr/bin/env python3
#
# vustruct_webdate
#
# A second flask program which
# - runs local alongside the wordpress server and
# - implements RESTAPI calls which:
#   1. record job identifiers which are in process
#   2. rebuild local websites in file system from emerging vustruct .tar.gz files

import argparse
import logging
from flask import Flask
from flask import request
from flask import jsonify
import uuid
import sys
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

# UDN="/dors/capra_lab/users/mothcw/UDNtests/"
UDN = os.getenv("UDN")
if not UDN:
    errormsg = "UDN environment variable must be defined before launch of %s" % UDN
    sys.exit(errormsg)



cmdline_parser = argparse.ArgumentParser()
cmdline_parser.add_argument(
    "-v", "--verbose",
    help="Include routine info log entries on stderr", default=False, action="store_true")
cmdline_parser.add_argument(
    "-d", "--debug",
    help="Include routine info AND 'debug' log entries on stderr",
    default=False, action="store_true")

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

LOGGER.info("%s running web-input cases in %s", __file__, UDN)
LOGGER.info("Initializing Flask(%s)", __name__)

print("HELLO")

app = Flask(__name__)

print("HELLO after Flask")

# Create a mapping of the unique uuids for running jobs to user and job names
# so that we know where to go fish out emerging websites.`
jobs_dict = {}
jobs_needing_website_refresh = []

@app.route('/add_uuid', methods=['POST', 'GET'])
def add_uuid():
    print("add_uuid(): The request I got is: %s" % request.json)

    case_url = "None"

    if 'job_uuid' in request.json and 'case_id' in request.json:
        case_url = os.path.join("http://localhost/vustruct/",
                                request.json['job_uuid'],
                                request.json['case_id'])

    webpage_home = os.path.join("/var/www/html/vustruct/",
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

    return case_url

