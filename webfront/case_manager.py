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

    def __init__(self, UDN: str, case_id: str, case_uuid: str) -> None:
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

    def to_dict(self) -> dict:
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

    def fetch_input_file(self,file_extension_including_period: str) -> str:
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
             "%s%s" % (self.external_user_prefix, file_extension_including_period))

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

        launch_plan_command = "psb_plan.py"
        app.logger.info("%s Running: %s" , threading.get_ident(), launch_plan_command)

        launch_plan_return = \
            subprocess.run(launch_plan_command, shell=False, encoding='UTF-8',
                       cwd = self.working_directory,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        app.logger.debug("%s %s finshed with exit code %s", threading.get_ident(), launch_plan_command, launch_plan_return.returncode)


        return launch_plan_return.returncode

