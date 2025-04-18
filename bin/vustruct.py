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
#
# class VUStruct, in essence, oversees all the individual VUStruct
# pipeline command line components:
# 'preprocess', 'plan', 'launch', 'report'
#
# As "side effect" of insantiating this class, an error exit for unhandled exceptions is created
#
# Each psb_*.py command line program "reports in" to class VUStruct
# as progress is made.  VUStruct stores records of the invocation
# of each vustruct main program, the location of generated  log files,
# and other details helpful to a final website, which must include
# user-accessible log files to replace tradition .log files accessible
# from the computer cluster filesystem.
#
# These storage records are written to the file $case_id + "_vustruct" + ".json"
#
# The individual pipeline programs tracked by class VUStruct are:
#
# 0) Optional Genomic Coordinates preprocessing:
#         parse_udn_report.py, parse_wustl.py, or vcf2vustruct.py
# 1) psb_plan.py      Structure analysis and job planning
# 2) psb_launch.py    Launch planned jobs on the cluster
# 3) psb_monitor.py   (checks vustruct.launch_successful)
# 4) psb_report.py    Generate user-facing webpage report
#
# In addition to various member functions which can answer the status of the case
# and early progress results...
# The final 'output' of this class is a dictionary that can be fed to
# a jinja2 template as part of website creation, to display the various
# logs and so forth from the run of all the modules.
# 
# class VUStruct is independent of class VUS_flask_case_manager (vustruct_flask.py)
# VUStruct instances evolve on disk in the same way, whether applications are 
# run at the command line or via the website (vustruct_flask.py)
# =============================================================================#

import logging
from logging.handlers import RotatingFileHandler
# from logging.handlers import FileHandler

import json
import copy
import sys
import os
from datetime import datetime

class VUStruct:

    def __init__(self, vustruct_module: str, case_id: str, executable_file: str):
        """Create logging streams compatible with web-site building in the log/
           directory.  Also setup an exception handler that routes unhandled
           exceptions to the logs, on the way to program exit"""

        self._stream_handler = None             # Outputs to stderr I think
        self._rotating_file_handler = None      # Detailed logs
        self._critical_error_handler = None     # File in errorlogs/ useful to display to users on abrupt end of psb_plan.py

        self._psb_plan_temp_filehandler = None

        possible_vustruct_modules = [
            'preprocess',
            'plan',
            'launch',
            'monitor',
            'report']

        assert vustruct_module in possible_vustruct_modules

        self._module = vustruct_module

        self._filename = case_id + "_vustruct.json"

        # Create a dictioary of dictionaries, empty, for each of the possible modules that are run
        # from the command line
        self._empty_vustruct_dict = {
            command_line_module: {
                'executable': '',  # Empty string until a module is _attempted_
                'start_time': '',
                'end_time': '',
                'log_filename': '',
                'critical_error_log_filename': '',
                'input_filename': '',
                'outputs': {},
                'exit_code': 1  # Default file is marked as not exiting well

            } for command_line_module in possible_vustruct_modules
        }

        self._vustruct_dict = copy.deepcopy(self._empty_vustruct_dict)

        self.read_file()

        self.log_filename = ''
        self.critical_error_log_filename = ''
        self._vustruct_dict[self._module]['executable'] = executable_file
        sys.excepthook = VUStruct.handle_unhandled_exception_and_exit

        old_log_record_factory = logging.getLogRecordFactory()
        def vustruct_log_record_factory(*args, **kwargs):
            record = old_log_record_factory(*args,**kwargs)
            record.vustruct_case_id = case_id
            return record

        logging.setLogRecordFactory(vustruct_log_record_factory)


    @staticmethod
    def handle_unhandled_exception_and_exit(exc_type, exc_value, exc_traceback):
        logging.getLogger().error("Uncaught exception %s %s"  % ( exc_type.__name__, exc_value),  exc_info=(exc_type, exc_value, exc_traceback))
        sys.exit(1)

    @property
    def standard_log_formatter(self):
        return logging.Formatter(
            '%(asctime)s %(levelname)-8s [%(filename)20s:%(lineno)-4d] %(message)s',
            datefmt="%H:%M:%S")

    @property
    def critical_error_log_formatter(self):
        return logging.Formatter(
            '-'*80 + '\n' + 
            'VUStruct Case %(vustruct_case_id)s\n%(levelname)s error while processing input: %(asctime)s in program %(filename)s at line %(lineno)d\n%(message)s' +
            '\n' + '-'*80,
            datefmt="%H:%M:%S")


    def initialize_file_and_stderr_logging(self, debug: bool) -> str:
        """
        For many applications, we need a rotating file handler in log/
        named for the mainline program.  Return the name of the create log_filename
        """
        assert self._rotating_file_handler == None
        assert self._stream_handler == None
        assert self._critical_error_handler == None

        program_name = os.path.splitext(os.path.basename(self.executable_file))[0]
        self.log_filename = os.path.join("./log", "%s.log" % program_name)
        self.critical_error_log_filename = os.path.join("./log", "%s.critical_error.log" % program_name)

        # If the main log needs a file roll, note that and we'll get to it...
        need_roll = os.path.isfile(self.log_filename)

        # Build out any directories down to the ./log directory
        os.makedirs(os.path.dirname(self.log_filename), exist_ok=True)
        sys.stderr.write("VUStruct log file set to: %s\n" % self.log_filename)

        root_logger = logging.getLogger()
        root_logger.setLevel(logging.DEBUG)

        # We typically log INFO and above to console
        self._stream_handler = logging.StreamHandler()
        self._stream_handler.setLevel(logging.DEBUG if debug else logging.INFO)
        log_formatter = logging.Formatter('%(asctime)s %(levelname)-8s [%(filename)20s:%(lineno)-4d] %(message)s',
                                      datefmt="%H:%M:%S")
        self._stream_handler.setFormatter(self.standard_log_formatter)
        root_logger.addHandler(self._stream_handler)

        # In the saved file, we capture DEBUG level and above
        self._rotating_file_handler = RotatingFileHandler(self.log_filename, backupCount=7)
        self._rotating_file_handler.setLevel(logging.DEBUG if debug else logging.INFO)
        root_logger.addHandler(self._rotating_file_handler)
        self._rotating_file_handler.setFormatter(self.standard_log_formatter)

        if need_roll:
            self._rotating_file_handler.doRollover()

        root_logger.info("Log file opened by %s", self.executable_file)

        log_symlink = os.path.basename(self.log_filename)
        try:
            os.remove(log_symlink)
        except FileNotFoundError:
            pass
        except OSError as ex:
            if (ex.errno != errno.ENOENT):
                raise

        try:
            os.symlink(self.log_filename, log_symlink)
        except:
            root_logger.info(f"Unable to complete os.symlink('{self.log_filename}', '{log_symlink}')")
            pass

        # We need one more log file, this one ONLY for the critical and error entries
        # We create this file anew with every run - not an append
        self._critical_error_handler = logging.FileHandler(self.critical_error_log_filename, mode='w')
        self._critical_error_handler.setLevel(logging.ERROR)
        root_logger.addHandler(self._critical_error_handler)
        self._critical_error_handler.setFormatter(self.critical_error_log_formatter)


        return self.log_filename

    def set_custom_log_formatter(self, log_format_string: str, date_format_string: str):
        _log_formatter = logging.Formatter(log_format_string,
                                      datefmt=date_format_string)

        self._rotating_file_handler.setFormatter(_log_formatter)
        self._stream_handler.setFormatter(_log_formatter)

    def add_temp_logger_to_psb_plan(self, mutation_log_dir: str) -> str:
        """
        It is very convenient for psb_plan to open an echo .log file in each 
        subdirectory as it works through the missense.csv list...
        """
        assert self._psb_plan_temp_filehandler == None
        mutation_log_filename = os.path.join(mutation_log_dir, "psb_plan.log")
        _needRoll = os.path.isfile(mutation_log_filename)

        self._psb_plan_temp_filehandler = RotatingFileHandler(mutation_log_filename, backupCount=7)
        self._psb_plan_temp_filehandler.setFormatter(self.standard_log_formatter)
        self._psb_plan_temp_filehandler.setLevel(logging.INFO)
        logging.getLogger().addHandler(self._psb_plan_temp_filehandler)
        if _needRoll:
            self._psb_plan_temp_filehandler.doRollover()
        return mutation_log_filename

    def remove_temp_logger_from_psbplan(self):
        # Close out the local log file for this mutation
        self._psb_plan_temp_filehandler.flush()
        self._psb_plan_temp_filehandler.close()
        logging.getLogger().removeHandler(self._psb_plan_temp_filehandler)
        self._psb_plan_temp_filehandler = None

    def dict_for_jinja2(self):
        # Return all information about all planning and launch phases
        # for report generation
        jinja2_dict = copy.deepcopy(self._vustruct_dict)
        return jinja2_dict

    def stamp_start_time(self):
        self._vustruct_dict[self._module]['start_time'] = datetime.now().isoformat()

    @property
    def start_time(self):
        return self._vustruct_dict[self._module]['start_time']

    def stamp_end_time(self):
        self._vustruct_dict[self._module]['end_time'] = datetime.now().isoformat()

    @property
    def end_time(self):
        return self._vustruct_dict[self._module]['end_time']


    @property
    def exit_code(self):
        return self._vustruct_dict[self._module]['exit_code']

    @exit_code.setter
    def exit_code(self, value: int):
        assert type(value) == int
        self._vustruct_dict[self._module]['exit_code'] = value

    @property
    def log_filename(self):
        return self._vustruct_dict[self._module]['log_filename']

    @log_filename.setter
    def log_filename(self, log_filename):
        self._vustruct_dict[self._module]['log_filename'] = log_filename

    @property
    def input_filename(self):
        return self._vustruct_dict[self._module]['input_filename']

    @input_filename.setter
    def input_filename(self, input_filename):
        self._vustruct_dict[self._module]['input_filename'] = input_filename

    @staticmethod
    def _read_failure_info_file(self_preprocess_or_plan) -> list[str]:
        critical_error_log_filename = os.path.join("./log", 
            "%s.critical_error.log" % os.path.basename(self_preprocess_or_plan['executable']).split('.')[0])

        failure_info = []
        try:
            with open(critical_error_log_filename) as f:
                for line in f:
                    failure_info.append(line.rstrip())
        except Exception as ex:
            failure_info =  ["Unable to open %s:" % critical_error_log_filename] + str(ex).split('\n')

        return failure_info

    @property
    def preprocess(self) -> dict:
        return self._vustruct_dict['preprocess']

    @property
    def preprocess_attempted(self) -> bool:
        return bool(self.preprocess['executable']) 

    @property 
    def preprocess_failed(self) -> bool:
        """
        Return True if one of the preprocessors
        was run, and exited non-zero, indicating failure
        """
        return (self.preprocess_attempted and 
               self.preprocess['exit_code'] != 0)

    @property 
    def preprocess_succeeded(self) -> bool:
        """
        Return True if one of the preprocessors
        was run, and exited zero, indicating success
        """
        return (self.preprocess_attempted and 
                self.preprocess['exit_code'] == 0)

    @property 
    def preprocess_failure_info(self) -> list[str]:
        return self._read_failure_info_file(self.preprocess)

    @property
    def plan_attempted(self) -> bool:
        return bool(self.plan['executable']) 

    @property 
    def plan_failed(self) -> bool:
        """
        Return True if vustruct_plan.py 
        was run, but exited non-zero, indicating failure
        """
        return (self.plan_attempted and 
               self.plan['exit_code'] != 0)


    @property 
    def plan_succeeded(self) -> bool:
        """
        Return True if plan
        was run, and exited zero, indicating success
        """
        return (self.plan_attempted and 
                self.plan['exit_code'] == 0)

    @property
    def plan_failure_info(self) -> list[str]:
        return self._read_failure_info_file(self.plan)

    @property
    def plan(self):
        return self._vustruct_dict['plan']

    @property
    def launch_attempted(self) -> bool:
        return bool(self.launch['executable']) 

    @property 
    def launch_failed(self) -> bool:
        """
        Return True if vustruct_launch.py 
        was run, but exited non-zero, indicating failure
        """
        return (self.launch_attempted and 
               self.launch['exit_code'] != 0)

    @property 
    def launch_succeeded(self) -> bool:
        """
        Return True if plan
        was run, and exited zero, indicating success
        """
        return (self.launch_attempted and 
                self.launch['exit_code'] == 0)

    @property
    def launch(self):
        return self._vustruct_dict['launch']

    @property
    def monitor(self):
        return self._vustruct_dict['monitor']

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
