#!/usr/bin/env python
"""Class PsbShared exposes static functions that psb_*.py modules employ
to handle command line argument processing, default parameters, and config file reads"""

import argparse
import configparser
import os
import stat
import getpass
import sys
import inspect

import logging

LOGGER = logging.getLogger(__name__)

statusdir = None


def create_default_argument_parser(callers_docstring: str,
                                   global_config_dirname: str = "../config",
                                   user_config_dirname: str = "..", add_help: bool = True) -> argparse.ArgumentParser:
    """All psb_pipeline command lines share common features which are initialized
    in this function.  These include auto-setting of the default case name to the
    current last component of the path, setting of a global configuration file based on
    whether the environment is development or production, and setting the default user
    configuration overrides to ../username.config"""

    script_realpath = os.path.realpath(__file__)

    # As of April 2021, the global config file will always be searched first in $UDN/config (i.e. ../config)
    if global_config_dirname != "../config":
        if os.path.exists("../config/global.config"):
            global_config_dirname = "../config/"

    # is_production = global_config_dirname.endswith("/capra_lab/users/psbadmin")

    default_global_config = None
    default_global_config = os.path.join(
        global_config_dirname,
        "global.config") # if is_production else "dev_global.config")
    cmdline_parser = argparse.ArgumentParser(
        description=callers_docstring,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=add_help)
    cmdline_parser.add_argument(
        "-c", "--config",
        help="PDBMap configuration profile for database access and global parameters",
        required=False, metavar="FILE", default=default_global_config)
    cmdline_parser.add_argument(
        "-u", "--userconfig",
        help="User specific settings and configuration profile overrides",
        required=False, metavar="FILE", default=os.path.join(user_config_dirname, getpass.getuser() + ".config"))
    cmdline_parser.add_argument(
        "-g", "--caseconfig",
        help="Case-specific settings and configuration profile overrides",
        required=False, metavar="FILE",
        default=os.path.join(".", "%s.config" % os.path.basename(os.getcwd()).replace(os.path.sep, '')))
    cmdline_parser.add_argument(
        "-v", "--verbose",
        help="Include routine info log entries on stderr", default=False, action="store_true")
    cmdline_parser.add_argument(
        "-d", "--debug",
        help="Include routine info AND 'debug' log entries on stderr",
        default=False, action="store_true")

    return cmdline_parser


def read_config_files(args, required_config_items=None):
    """All psb_pipeline applications parse global and user-specific config files.  This
    function cals ConfigParser.SafeConfigParser to get the job done"""

    config = configparser.SafeConfigParser(allow_no_value=True)

    config_file_list = [args.config, args.userconfig, args.caseconfig]
    LOGGER.info("Attempting to read config from " + str(config_file_list))
    parsed_file_list = config.read(config_file_list)
    LOGGER.info("Successfully read config from " + str(parsed_file_list))
    if len(parsed_file_list) < 1:
        LOGGER.critical(
            """Config files %s were not found.
            Parser cannot proceed without at least one valid config file.""",
            str(config_file_list))
        sys.exit(1)

    # config.items() returns a list of (name, value) pairs - convert to dict
    config_dict = dict(config.items("Genome_PDB_Mapper"))

    # It is an odd thing that we poulate from a UserSpecific section - fix later
    try:
        config_dict.update(dict(config.items("UserSpecific")))
    except (KeyError,configparser.NoSectionError):
        LOGGER.warning("No UserSpecific section was found in the read configuration files")

    if required_config_items:
        missing_keys = [name for name in required_config_items if name not in config_dict]

        if len(missing_keys):
            LOGGER.error("Can't proceed without configuration file options set for: %s", str(missing_keys))
            sys.exit(1)

    return config, config_dict
