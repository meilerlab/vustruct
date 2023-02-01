#!/usr/bin/env python3
#
# Project        : PSB Pipeline
# _filename       : vustruct.py
# Authors        : Chris Moth
# Organization   : Meiler Lab
#                : Vanderbilt University
# Email          : chris.moth@vanderbilt.edu
# Date           : 2023-01
# Description
# class VUstruct, in essence, oversees all the individual VUstruct
# pipeline command line components.
#
# Each psb_*.py command line program "reports in" to class VUstruct
# as progress is made.  VUstruct stores records of the invocation
# of each vustruct main program, the location of generated  log files,
# and other details helpful to a final website, which must include
# user-accessible log files to replace tradition .log files accessible
# from the computer cluster filesystem.
#
# These storage records are written to the file $case_id + "_vustruct" + ".json"
#
# The individual pipeline programs tracked by class VUstruct are:
#
# 0) Optional Genomic Coordinates preprocessing:
#         parse_udn_report.py, parse_wustl.py, or vcf2missense.py
# 1) psb_plan.py
# 2) psb_launch.py
# 3) psb_monitor.py
# 3) psb_report.py

# =============================================================================#

import json
import copy
import sys
from datetime import datetime
import logging
LOGGER = logging.getLogger(__name__)


class VUstruct:

    def __init__(self, vustruct_module: str, case_id: str, executable_file: str):
        possible_vustruct_modules = [
            'preprocess',
            'plan',
            'launch',
            'report']

        assert vustruct_module in possible_vustruct_modules

        self._module = vustruct_module

        self._filename = case_id + "_vustruct.json"

        # Create a dictioary of dictionaries, empty, for each of the possible modules that are run
        # from the command line
        self._empty_vustruct_dict = {
            command_line_module: {
                'executable': '',
                'start_time': '',
                'logfile': '',
                'exit_code': 1  # Default file is marked as not exiting well

            } for command_line_module in possible_vustruct_modules
        }

        self._vustruct_dict = copy.deepcopy(self._empty_vustruct_dict)

        self.read_file()

        self.logfile = ''
        self._vustruct_dict[self._module]['executable'] = executable_file

    # The final 'output' of this class is a dictionary that can be fed to
    # a jinja2 template as part of website creation, to display the various
    # logs and so forth from the run of all the modules.
    def dict_for_jinja2(self):
        return copy.deepcopy(self._vustruct_dict)

    def stamp_start_time(self):
        self._vustruct_dict[self._module]['start_time'] = datetime.now().isoformat()

    @property
    def start_time(self):
        return self._vustruct_dict[self._module]['start_time']

    @property
    def exit_code(self):
        return self._vustruct_dict[self._module]['exit_code']

    @exit_code.setter
    def exit_code(self, value):
        self._vustruct_dict[self._module]['exit_code'] = value

    @property
    def logfile(self):
        return self._vustruct_dict[self._module]['logfile']

    @logfile.setter
    def logfile(self, logfile_name):
        self._vustruct_dict[self._module]['logfile'] = logfile_name

    @property
    def preprocess(self):
        return self._vustruct_dict['preprocess']

    @property
    def executable_file(self):
        return self._vustruct_dict[self._module]['executable']

    def read_file(self):
        try:
            with open(self._filename, 'r') as f:
                self._vustruct_dict = json.load(f)
        except FileNotFoundError:
            self._vustruct_dict = copy.deepcopy(self._empty_vustruct_dict)

    def write_file(self):
        with open(self._filename, 'w') as f:
            json.dump(self._vustruct_dict, f, indent=2)


if __name__ == '__main__':
    print("Only for use as library")
    sys.exit(1)
