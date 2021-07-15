#!/usr/bin/env python3
#
# Project        : PSB Pipeline
# Filename       : psb_launch.py
# Authors        : Chris Moth and R. Michael Sivley
# project Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : chris.moth@vanderbilt.edu
# Date           : 2018-01-22
# Description    : Launch jobs specified in the input .csv file
#                : (output from psb_plan.py)
# =============================================================================#

"""\
Launch pipeline jobs on a Slurm or LSF cluster given the UDN case ID
   or a specific workplan.csv previously output from psb_plan.py

   For help
        psb_plan.py --help
"""

import logging
import os
import sys
import grp
import stat
import copy

# import time
from datetime import datetime
import pprint
from logging.handlers import RotatingFileHandler
from typing import Dict, Any, Tuple, List
from typing import TextIO

# Inspect to recover bsub_submit() and slurm_submit() source code for integration
# into generated final file
import inspect

# Dataframes
import pandas as pd
# from jinja2 import Environment, FileSystemLoader

from slurm import slurm_submit
from bsub import bsub_submit

from psb_shared import psb_config
from psb_shared import psb_perms


# THIS BELOW IS CRAP
try:
    capra_group = grp.getgrnam('capra_lab').gr_gid
except KeyError:
    capra_group = os.getegid()


# THIS ABOVE IS CRAP


def set_capra_group_sticky(dirname):
    return
#    try:
#        os.chown(dirname, -1, capra_group)
#    except OSError:
#        pass
#
#    # Setting the sticky bit on directories also fantastic
#    try:
#        os.chmod(dirname,
#                 stat.S_IRUSR|stat.S_IWUSR|stat.S_IXUSR|stat.S_IRGRP|stat.S_IWGRP|stat.S_IXGRP|stat.S_ISGID)
#    except OSError:
#        pass
#

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

cmdline_parser = psb_config.create_default_argument_parser(__doc__, os.path.dirname(os.path.dirname(__file__)))

cmdline_parser.add_argument("projectORworkplan", type=str,
                            help="Project ID (ex. UDN123456) or single output file from psb_plan.py\
                            Example: ......$UDN/UDN123456/GeneName_NM_12345.1_S123A_workplan.csv",
                            default=os.path.basename(os.getcwd()), nargs='?')
cmdline_parser.add_argument("-r", "--relaunch",
                            help="Relaunch jobs, ignoring any results in a previously generated workstatus.csv file",
                            action="store_true")
# cmdline_parser.add_argument("-l", "--lsf",
#                            help="Create LSF .bsub files instead of Slurm .slurm.", default=False,
#                             action="store_true")
cmdline_parser.add_argument("-n", "--nolaunch",
                            help="Create slurm(LSF) files, but do NOT exec sbatch(bsub) to launch them.", default=False,
                            action="store_true")

args, remaining_argv = cmdline_parser.parse_known_args()

if args.debug:
    sh.setLevel(logging.DEBUG)
elif args.verbose:
    sh.setLevel(logging.INFO)
# If neither option, then logging.WARNING (set at top of this code) will prevail for stdout

# A period in the argument means the user wants to launch one mutation only,
# directly from a single mutation output file of psb_plan.py
oneMutationOnly = ('.' in args.projectORworkplan and os.path.isfile(args.projectORworkplan))

required_config_items = ['output_rootdir', 'collaboration', 'ddg_config']

config, config_dict = psb_config.read_config_files(args, required_config_items)

LOGGER.info("Command: %s", ' '.join(sys.argv))

# Most will launch on Slurm
# WUSTL is LSF
# Also, someday, hook in launching on high core boxes directly
CLUSTER_TYPE = 'LSF' if 'cluster_type' in config_dict and config_dict['cluster_type'].upper() == 'LSF' else 'Slurm'
CONTAINER_TYPE = None
SINGULARITY_IMAGE = None
if 'container_type' in config_dict:
    if config_dict['container_type'].upper() == 'DOCKER':
        CONTAINER_TYPE = 'Docker'
    elif config_dict['container_type'].upper() == 'SINGULARITY':
        CONTAINER_TYPE = 'Singularity'
        SINGULARITY_IMAGE = os.getenv('SINGULARITY_CONTAINER')
        if not SINGULARITY_IMAGE:
            msg = "SINGULARITY_CONTAINER is not defined in the launch environment.  container_type of Singularity is nonsensical"
            LOGGER.critical(msg)
            sys.exit(msg)
    else:
        sys.exit("container_type set to %s in config file.  Must be Docker or Singularity" %
                 config_dict['container_type'])


try:
    launchParametersAll: Dict[str, Any] = dict(config.items("%sParametersAll" % CLUSTER_TYPE))
except KeyError:
    launchParametersAll: Dict[str, Any] = {}

LOGGER.info("Command Line Arguments:\n%s", pprint.pformat(vars(args)))

launchParametersPathProx: Dict[str, Any] = copy.deepcopy(launchParametersAll)
launchParametersDDGMonomer: Dict[str, Any] = copy.deepcopy(launchParametersAll)
launchParametersUDNSequence: Dict[str, Any] = copy.deepcopy(launchParametersAll)
launchParametersDigenicAnalysis: Dict[str, Any] = copy.deepcopy(launchParametersAll)

# Initialize the slurm/LSF parameters for each of the 4 possible calculation sets
# Log the parameters (gathered from the config file
# Check the parameters that they have all required settings
slurm_required_settings = ['mem', 'ntasks', 'time']
LSF_required_settings = [
    ('G', 'user_group'),
    ('q', 'queue_name'),
    # ('a','application_name'),
    # ('J','job_name'), # Populated above
    ('u', 'mail_user'),
    ('W', 'limit'),
    ('R', 'resource'),
    ('n', 'tasks')
]
for launchParameterDictDesc, launchParameters in [
    ('ParametersPathProx', launchParametersPathProx),
    ('ParametersDDGMonomer', launchParametersDDGMonomer),
    ('ParametersUDNSequence', launchParametersUDNSequence),
    ('ParametersDigenicAnalysis', launchParametersDigenicAnalysis)
]:
    # Update each flavor-specific dictionaryof parameters
    # with the flavor-specific overrides in the config file
    slurm_or_lsf_section_name = CLUSTER_TYPE + launchParameterDictDesc
    config_items = {}
    if slurm_or_lsf_section_name in config:
        launchParameters.update(dict(config.items(slurm_or_lsf_section_name)))

    LOGGER.info("%s:\n%s", slurm_or_lsf_section_name, pprint.pformat(launchParameters))

    # Now make _sure_ that required parameters are in the (updated) launchParameters dictionary.
    if CLUSTER_TYPE == 'Slurm':
        for req in slurm_required_settings:
            if req not in launchParameters:
                LOGGER.error(
                    """\
                    Can't launch jobs because you have not provided the required\n%s setting [%s] for section %s
                    (or %sParametersAll)\n in either file %s or %s\n""", 
                    CLUSTER_TYPE, req, CLUSTER_TYPE, slurm_or_lsf_section_name, args.config, args.userconfig
                )
                sys.exit(1)
    elif CLUSTER_TYPE == 'LSF':
        for req in LSF_required_settings:
            if req[0] not in launchParameters and req[1] not in launchParameters:
                LOGGER.error(
                    """\
                    Can't launch jobs because you have not provided the required\n%s setting [%s] for section %s
                    (or %sParametersAll)\n in either file %s or %s\n""",
                    CLUSTER_TYPE, req, CLUSTER_TYPE, slurm_or_lsf_section_name, args.config, args.userconfig
                )
                sys.exit(1)
    else:
        sys.exit("cluster_type=%s not supported" % CLUSTER_TYPE)

# The collaboration_dir is the master directory for the case, i.e. for one patient
# Example: /dors/capra_lab/projects/psb_collab/UDN/UDN532183
udn_root_directory = os.path.join(config_dict['output_rootdir'], config_dict['collaboration'])

# Mounted directories refer to filesystem locations _outside_ the docker
# But which will be known to the Slurm and LSF environments prior to container start
mounted_UDN_directory = None
mounted_collaboration_directory = None

if 'mounted_udn_directory' in config_dict:
    mounted_UDN_directory = config_dict['mounted_udn_directory']

if oneMutationOnly:
    # Get the parent directory above the workplan
    collaboration_dir = os.path.dirname(os.path.dirname(os.path.abspath(args.projectORworkplan)))
else:
    collaboration_dir = os.path.join(udn_root_directory, args.projectORworkplan)
    if mounted_UDN_directory:
        mounted_collaboration_directory = os.path.join(mounted_UDN_directory, args.projectORworkplan)

psb_permissions = psb_perms.PsbPermissions(config_dict)
psb_permissions.makedirs(collaboration_dir)


class JobsLauncher:
    """
    Create .slurm or .bsub files for a dataframe of jobs.
    """

    def __init__(self, workplan_or_case_filename: str, df_all_jobs: pd.DataFrame, df_prior_exit0: pd.DataFrame) -> None:
        """
        @param workplan_or_case_filename: Retained for output to slurm/bsub comments section.
        @param df_all_jobs:     Dataframe from workplan.csv, created by psb_plan.py.  IT contains all jobs and params
        @param df_prior_exit0:  Dataframe from workstatus.csv, containing rows where previous jobs exited
                                with code 0 (success)
        """

        self._workplan_or_case_filename = workplan_or_case_filename
        self._df_all_jobs = df_all_jobs
        self._df_prior_exit0 = df_prior_exit0

        row0 = self._df_all_jobs.iloc[0]
        self._gene = row0['gene']
        self._refseq = row0['refseq']
        self._mutation = row0['mutation']

        self._geneRefseqMutation_OR_casewideString = 'casewide'
        # If gene is not casewide run, then fill it in specifically
        if self._gene != self._geneRefseqMutation_OR_casewideString:
            self._geneRefseqMutation_OR_casewideString = "%s_%s_%s" % (
                self._gene, self._refseq, self._mutation)  # This is usual case for the gene by gene mutation run sets

        self._slurm_or_bsub = 'bsub' if CLUSTER_TYPE == 'LSF' else 'slurm'

        # It is convenient to harvest the current working directory from the workplan/all_jobs_df
        # It should be the case that the 'cwd' column is the same across ALL jobs to be run.
        self._all_jobs_cwd = None

        # Whether by slurm, LSF or other, each job has s specific command line with arguments that must be run.
        self._launch_strings = {'PathProx': [], 'ddG_monomer': [], 'SequenceAnnotation': [], 'DigenicAnalysis': []}

        # After jobs are launched, we will create a 'submitted' file in the status directories for the jobs
        # Might collide with launched program though - so perhaps rethink a bit....
        self._status_dir_list = []

        # Dictionares to match each launched job unique key to it';s job id and array id
        self._jobids = {}
        self._arrayids = {}
        self._jobinfo = {}
        self._jobprogress = {}
        self._exitcodes = {}

    @property
    def mutation_dir(self):
        return os.path.join(collaboration_dir, self._geneRefseqMutation_OR_casewideString)

    def _build_a_launch_dir(self, flavor: str) -> Tuple[str,str]:
        launch_dir_inside_container = os.path.join(self.mutation_dir,
                            self._slurm_or_bsub,
                            "PathProx" if "_PP_" in flavor else flavor
                            )
        if mounted_collaboration_directory:
            launch_dir_outside_container = os.path.join(mounted_collaboration_directory,
                            self._geneRefseqMutation_OR_casewideString,
                            self._slurm_or_bsub,
                            "PathProx" if "_PP_" in flavor else flavor
                            )
        else:
            launch_dir_outside_container = None

        return launch_dir_inside_container, launch_dir_outside_container

    def _add_launch_directory_to_df_all_jobs(self) -> None:
        """
        Add a new column in the all_jobs dataframe, which
        will be named either 'slurm_directory' or 'bsub_directory'
        @return:  None
        """
        self._df_all_jobs['launch_directory'] = self._df_all_jobs.apply(
            lambda a_row: self._build_a_launch_dir(a_row['flavor'])[0],
            axis=1
        )

    def _add_gene_refseq_mutation_flavor_to_df_all_jobs(self):
        self._df_all_jobs['gene_refseq_mutation_flavor'] = self._df_all_jobs.apply(
            lambda a_row: a_row['flavor'] if a_row['gene'] == 'casewide' else
            "%s_%s_%s_%s" % (a_row['gene'], a_row['refseq'], a_row['mutation'], a_row['flavor']),
            axis=1
        )

    def _add_launch_filename_to_df_all_jobs(self):
        self._df_all_jobs['launch_file'] = self._df_all_jobs.apply(
            lambda a_row: os.path.join(a_row['launch_directory'], a_row['gene_refseq_mutation_flavor'] + '.' + self._slurm_or_bsub),
            axis=1
        )

    def _create_launch_strings(self):
        # Prepare directories for non-ddG jobs
        # Also create launch_strings for subsequence integration into slurm files
        # and launches
        self._all_jobs_cwd = None
        for _, row in self._df_all_jobs.iterrows():
            if (not self._df_prior_exit0.empty) and (row['uniquekey'] in self._df_prior_exit0.index):
                LOGGER.info("Job %s was successful previously.  %s file will not be recreated", row['uniquekey'], self._slurm_or_bsub)
            else:
                assert self._all_jobs_cwd is None or self._all_jobs_cwd == row['cwd']
                self._all_jobs_cwd = row['cwd']
                # Create the meat of the slurm script for this job
                # ddG monomer is simpler because it is not managed by the pipeline
                launch_string = "%(command)s %(options)s" % row
                LOGGER.info("%s %s", row['flavor'], launch_string)
                if 'ddG' not in row['flavor']:
                    launch_string += " --outdir %(outdir)s/%(flavor)s --uniquekey %(uniquekey)s" % row
                # else: << No just do this in the container
                #    launch_string = """\
                #export LD_LIBRARY_PATH=/psbadmin/rosetta3.7/main/source/build/src/release/linux/3.10/64/x86/gcc/4.8/default:\
                #/psbadmin/rosetta3.7/main/source/build/external/release/linux/3.10/64/x86/gcc/4.8/default; """ + launch_string

                uniquekey_launch_string = (row['uniquekey'], launch_string)
                # Launch PathProx Clinvar and COSMIC both with same parameters
                # out of same .slurm file
                if "PP_" in row['flavor']:
                    self._launch_strings['PathProx'].append(uniquekey_launch_string)
                elif row['flavor'] in self._launch_strings:  # Great for ddG, SequenceAnnocation, DigenicAnalysis
                    self._launch_strings[row['flavor']].append(uniquekey_launch_string)
                else:
                    print("I don't know this job flavor: ", row['flavor'])
                    sys.exit(1)

                # Before launching a job, make path to its
                # status directory and erase all files from that
                # directory if any there from previous run
                status_dir = os.path.join(row['outdir'], row['flavor'], "status")
                old_umask = os.umask(2)
                os.makedirs(status_dir, exist_ok=True)
                os.umask(old_umask)
                set_capra_group_sticky(status_dir)

                # Since we are launching anew, clear out ANY old junk files hanging around'
                # the status_dir of the jobs
                for the_file in os.listdir(status_dir):
                    file_path = os.path.join(status_dir, the_file)
                    LOGGER.info('Deleting old file %s', file_path)
                    try:
                        if os.path.isfile(file_path):
                            os.unlink(file_path,)
                    except OSError:
                        LOGGER.exception("Unable to delete %s in status directory", file_path)
                        sys.exit(1)

                # Well, we have not _actually_ submitted it _quite_ yet - but just save this dirname
                # to mark later
                self._status_dir_list.append(status_dir)

    def _update_job_statuses_from_prior_exit0s(self):
        for index, row in self._df_all_jobs.iterrows():
            if (not self._df_prior_exit0.empty) and (row['uniquekey'] in self._df_prior_exit0.index):
                LOGGER.info(
                    "Job %s was successful previously.  Prior results will be copied to new workstatus file", row[
                        'uniquekey'])
                prior_success = self._df_prior_exit0.loc[row['uniquekey']]
                self._jobinfo[row['uniquekey']] = prior_success['jobinfo']
                self._jobprogress[row['uniquekey']] = prior_success['jobprogress']
                self._jobids[row['uniquekey']] = prior_success['jobid']
                self._arrayids[row['uniquekey']] = prior_success['arrayid'] if ('arrayid' in prior_success) else 0
                self._exitcodes[row['uniquekey']] = prior_success['ExitCode']

    def _write_launch_independent_comments(self, launch_filename: str, slurm_f: TextIO, subdir: str):
        slurm_f.write("""\
#!/bin/bash
# Project        : PSB Pipeline
# Filename       : %s
# Generated By   : %s
# For work plan  : %s
# Organization   : Vanderbilt Genetics Institute,
#                : Program in Personalized Structural Biology,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University
# Generated on   : %s
# Description    : Runs a python script on a cluster to accomplish: %s
#===============================================================================
""" % (launch_filename, __file__, self._workplan_or_case_filename, str(datetime.now()), subdir))

    def _write_launch_independent_environment(self, launch_filename: str, slurm_f: TextIO, subdir: str):
        # In a container, our execution path should point to the container areas.  Annoyingly,
        # variations in Docker/Singlularity/LSF demand we explicitly restate the PATH here.
        # If NOT in a container, then init the PATH to be whatever it is outside the container.

        if CONTAINER_TYPE:
            path_statement="PATH=%s" % ':'.join([
                "/opt/conda/bin",
                "/ensembl/ensembl-git-tools/bin",
                "/psbadmin/bin",
                "/psbadmin/pdbmap",
                "/psbadmin/pathprox",
                "$PATH"])
            LOGGER.info("Setting PATH for container %s", CONTAINER_TYPE)
            LOGGER.info(path_statement)
            slurm_f.write(path_statement + '\n')
        else:
            path_statement = "PATH=%s:$PATH\n" % os.getenv('PATH')
            LOGGER.info("Setting PATH to match launch-time %s environment", __file__)
            LOGGER.info(path_statement)
            slurm_f.write(path_statement)

        if 'ensembl_registry' in config_dict:
            ensembl_registry = config_dict['ensembl_registry']
            LOGGER.info("Setting ENSEMBL_REGISTRY from .config: %s", ensembl_registry)
        else:
            ensembl_registry = os.getenv('ENSEMBL_REGISTRY')
            LOGGER.info("Setting ENSEMBL_REGISTRY from launch environment: %s", ensembl_registry)

        if ensembl_registry:
            slurm_f.write("""
#  Point the ENSEMBL PERL_API to the desired genome release
export ENSEMBL_REGISTRY=%s
echo ENSEMBL_REGISTRY=$ENSEMBL_REGISTRY
""" % ensembl_registry)
        else:
            LOGGER.warning("ENSEMBL_REGISTRY is not defined in .config or environment.  Omitting from %s",
                           launch_filename)
            slurm_f.write("echo WARNING: ENSEMBL_REGISTRY not defined in launch environment.  Omitted")

        slurm_f.write("""
# Prevent loading of an individual users ~/.local type python - really a mess for pipeline!
export PYTHONNOUSERSITE=x

# Don't buffer output, to get better errors if things go south
export PYTHONUNBUFFERED=x
""")


    def _write_launch_independent_cd(self, launch_filename: str, slurm_f: TextIO, subdir: str):
        slurm_f.write("""
cd %s
if [ $? != 0 ]; then
echo Failure at script launch: Unable to change to directory %s
exit 1
fi
""" % ( self._all_jobs_cwd,
        self._all_jobs_cwd))

    def _create_slurm_file(self, subdir: str, launch_filename: str, launch_stdout_directory: str):
        """
        @param subdir: The subdirectory under variant/slurm/ where action takes place  "PathProx", etc
        @param launch_filename: The name, with path, of the .slurm or .bsub file being created.
        @param launch_stdout_directory: "stdout" subdir of above.
        @return:
        """
        with open(launch_filename, 'w') as slurm_f:
            self._write_launch_independent_comments(launch_filename, slurm_f, subdir)

            slurm_f.write("\n# Slurm Parameters\n")
            # It is important to make shallow copies of the default dictionaries,
            # else mutations #7 can pick up dictionary entries set by mutation #4
            if "PathProx" in subdir:
                slurm_dict = dict(launchParametersPathProx)
            elif subdir == "ddG_monomer":
                slurm_dict = dict(launchParametersDDGMonomer)
            elif subdir == "SequenceAnnotation":
                slurm_dict = dict(launchParametersUDNSequence)
            elif subdir == "DigenicAnalysis":
                slurm_dict = dict(launchParametersDigenicAnalysis)
            else:
                print("I don't know what flavor(subdir) is: ", subdir)
                sys.exit(1)

            slurm_dict['output'] = "%s/%s.out" % (
                launch_stdout_directory, "%s_%%A_%%a" % self._geneRefseqMutation_OR_casewideString)
            slurm_dict['job-name'] = "%s_%s" % (self._geneRefseqMutation_OR_casewideString, subdir)
            job_count = len(self._launch_strings[subdir])

            # Build a slurm array file if jobCount > 1
            if job_count > 1:
                slurm_dict['array'] = "0-%d" % (job_count - 1)
            else:
                # We've had trobule with modifications to the default dictionaries getting passed forward.
                # Make sure this did not happen
                assert 'array' not in slurm_dict

            # We can set a max on the number of simultaneously running jobs in the array.  For now, don't do this
            # if jobCount > 10:
            #   slurmDict['array'] += "%%10"

            for slurm_field in ['job-name', 'mail-user', 'mail-type', 'ntasks', 'time', 'mem', 'account', 'output',
                                'array']:
                if slurm_field in slurm_dict and slurm_dict[slurm_field]:
                    slurm_f.write("#SBATCH --%s=%s\n" % (slurm_field, slurm_dict[slurm_field]))

            if job_count > 1:
                slurm_f.write("""\
echo "SLURM_ARRAY_TASKID="$SLURM_ARRAY_TASKID
""")
            slurm_f.write("""\
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
# echo "SLURMTMPDIR="$SLURMTMPDIR
echo "SLURM_SUBMIT_DIR = "$SLURM_SUBMIT_DIR
""")

            self._write_launch_independent_environment(launch_filename, slurm_f, subdir)

            self._write_launch_independent_cd(launch_filename, slurm_f, subdir)

            # Prepend the slurm file case/esac command with a 
            if CONTAINER_TYPE == 'Singularity':
                if 'singularity_bind' in config_dict:
                    slurm_f.write("export SINGULARITY_BIND=%s\n" % config_dict['singularity_bind'])
                else:
                    slurm_f.write("# SINGULARITY_BIND not written.  singularity_bind omitted from .config files\n")

                if 'singularity_pre_command' in config_dict:
                    slurm_f.write("%s\n" % config_dict['singularity_pre_command'])
                else:
                    slurm_f.write("# No singularity pre command written.  Add singularity_pre_command to .config if this is desired")

                container_exec_prefix = 'singularity exec %s ' % SINGULARITY_IMAGE
            else:
                container_exec_prefix = ''

            if job_count == 1:
                # No need to fiddle with the slurm case statement
                slurm_f.write(container_exec_prefix + self._launch_strings[subdir][0][1])
                slurm_f.write("\n")
            else:
                # Via the case/esac mechanism, select the specific job to run
                slurm_f.write("case $SLURM_ARRAY_TASK_ID in\n")
                slurm_array_id = 0
                for uniquekey_launchstring in self._launch_strings[subdir]:
                    # Write out the "case" tag
                    slurm_f.write("%d)\n" % slurm_array_id)
                    # Save this to put back later
                    self._arrayids[uniquekey_launchstring[0]] = slurm_array_id
                    slurm_array_id += 1
                    slurm_f.write(container_exec_prefix + uniquekey_launchstring[1])
                    # end the case tag
                    slurm_f.write("\n;;\n\n")
                slurm_f.write("esac\n")
        LOGGER.info("Slurm file created: %s", launch_filename)

    def _create_bsub_file(self, subdir: str, launch_filename: str, launch_stdout_directory: str) -> Tuple[str, int]:
        """
        @param subdir: The subdirectory under variant/bsub/ where action takes place  "PathProx", etc
        @param launch_filename: The name, with path, of the .bsub or .bsub file being created.
        @param launch_stdout_directory: "stdout" subdir of above.
        @return: job_name, job_count
        """

        _job_name = None
        _job_count = None
        with open(launch_filename, 'w') as bsub_f:
            self._write_launch_independent_comments(launch_filename, bsub_f, subdir)

            bsub_f.write("\n# LSF Parameters\n")
            # It is important to make shallow copies of the default dictionaries,
            # else mutations #7 can pick up dictionary entries set by mutation #4
            if "PathProx" in subdir:
                bsub_dict = dict(launchParametersPathProx)
            elif subdir == "ddG_monomer":
                bsub_dict = dict(launchParametersDDGMonomer)
            elif subdir == "SequenceAnnotation":
                bsub_dict = dict(launchParametersUDNSequence)
            elif subdir == "DigenicAnalysis":
                bsub_dict = dict(launchParametersDigenicAnalysis)
            else:
                print("I don't know what flavor(subdir) is: ", subdir)
                sys.exit(1)

            bsub_dict['output'] = "%s/%s.out" % (
                launch_stdout_directory, "%s_%%J_%%I" % self._geneRefseqMutation_OR_casewideString)
            _job_name = "%s_%s" % (self._geneRefseqMutation_OR_casewideString, subdir)
            _job_count = len(self._launch_strings[subdir])
            if _job_count > 1:
                _job_name += "[1-%d]" % _job_count
            bsub_dict['job_name'] = _job_name

            # Build a LSF array file if jobCount > 1
            if _job_count > 1:
                # Annoyingingly there is no mechanism inside a .bsub script to specify all array elements
                pass
            else:
                # We've had trobule with modifications to the default dictionaries getting passed forward.
                # Make sure this did not happen
                assert 'array' not in bsub_dict

            # min_tasks, max_tasks
            bsub_dict['n'] = 1

            # RESUME HERE WITH OTHER BSUB OPTIONS
            for bsub_field in [
                ('P', 'project_name'),
                ('G', 'user_group'),
                ('q', 'queue_name'),
                ('a', 'application_name'),
                ('J', 'job_name'),  # Populated above
                ('u', 'mail_user'),
                ('g', 'job_group_name'),
                ('W', 'limit'),
                ('R', 'resource'),
                ('n', 'tasks'),
                ('o', 'output')     # Populated above
             ]:
                bsub_value = None
                # Accept either LSF short description or sanity preserving long description
                # from CONFIG file
                if bsub_field[0] in bsub_dict:
                    bsub_value = bsub_dict[bsub_field[0]]
                elif bsub_field[1] in bsub_dict:
                    bsub_value = bsub_dict[bsub_field[1]]

                # For BSUB itself, use the short value
                # Comment it with the long value
                if bsub_value:
                    if bsub_field[0] == 'u':  # This is the email - so add a -N to ensure notify
                        bsub_f.write("#BSUB -N\n")
                    bsub_f.write("#     %s\n" % bsub_field[1])
                    bsub_f.write("#BSUB -%s %s\n" % (bsub_field[0], bsub_value))
                else:
                    bsub_f.write("#     %s/%s not provided in .config file(s)\n" % (bsub_field[1], bsub_field[0]))

            if _job_count > 1:
                bsub_f.write("""\
echo "LSB_JOBINDEX="$LSB_JOBINDEX
""")

            for lsf_environment_variable in [
                'LSB_JOBNAME', 'LSB_HOSTS', 'LSB_JOBPIDS', 'LSB_OUTDIR', 'LSB_PROJECT_NAME',
                'LSB_QUEUE', 'LSB_BIND_CPU_LIST', 'LS_SUBCWD'
            ]:

                bsub_f.write('echo "{0}="${0}\n\n'.format(lsf_environment_variable))

            if 'lsf_docker_volumes' in config_dict:
                bsub_f.write("export LSF_DOCKER_VOLUMES=%s\n" % config_dict['lsf_docker_volumes'])
            else:
                bsub_f.write("# LSF_DOCKER_VOLUMES not written.  lsf_docker_volumes omitted from .config files\n")

            self._write_launch_independent_environment(launch_filename, bsub_f, subdir)

            self._write_launch_independent_cd(launch_filename, bsub_f, subdir)

            if _job_count == 1:
                # No need to fiddle with the slurm case statement
                bsub_f.write(self._launch_strings[subdir][0][1])
                bsub_f.write("\n")
            else:
                # Via the case/esac mechanism, select the specific job to run
                bsub_f.write("case $LSB_JOBINDEX in\n")
                # Unlike slurm arrays, bsubs are indexed from 1
                bsub_array_id = 1
                for uniquekey_launchstring in self._launch_strings[subdir]:
                    # Write out the "case" tag
                    bsub_f.write("%d)\n" % bsub_array_id)
                    # Save this to put back later
                    self._arrayids[uniquekey_launchstring[0]] = bsub_array_id
                    bsub_array_id += 1
                    bsub_f.write(uniquekey_launchstring[1])
                    # end the case tag
                    bsub_f.write("\n;;\n\n")
                bsub_f.write("esac\n")
        LOGGER.info("BSUB file created: %s", launch_filename)
        return _job_name, _job_count

    def _df_updated_workstatus(self) -> pd.DataFrame:
        """
        Add new columns for jobid, arrayid, jobinfo, jobprogress, and ExitCode to the workstatus dataframe,
        which will be written back to disk by the mainline called
        @return: dataframe containing the workplan, extended with workstatus information from process launches.
        """
        df_updated_workstatus = self._df_all_jobs.merge(pd.DataFrame(
            list(self._jobids.items()),
            columns=['uniquekey', 'jobid']),
            on='uniquekey', how='left')
        df_updated_workstatus = df_updated_workstatus.merge(
            pd.DataFrame(list(self._arrayids.items()), columns=['uniquekey', 'arrayid']),
            on='uniquekey', how='left')
        df_updated_workstatus = df_updated_workstatus.merge(
            pd.DataFrame(list(self._jobinfo.items()), columns=['uniquekey', 'jobinfo']),
            on='uniquekey', how='left')
        df_updated_workstatus = df_updated_workstatus.merge(
            pd.DataFrame(list(self._jobprogress.items()), columns=['uniquekey', 'jobprogress']), on='uniquekey',
            how='left')
        df_updated_workstatus = df_updated_workstatus.merge(
            pd.DataFrame(list(self._exitcodes.items()), columns=['uniquekey', 'ExitCode']),
            on='uniquekey', how='left')

        return df_updated_workstatus.set_index(['uniquekey', 'jobid'], inplace=False)

    def doit(self) -> Tuple[pd.DataFrame, List[str]]:
        """
        @return: Workstatus dataframe (workplan with additional stats columns)
                 List of .slurm or .bsub files that may need to be launched outside the container.
        """

        self._add_launch_directory_to_df_all_jobs()
        self._add_gene_refseq_mutation_flavor_to_df_all_jobs()
        self._add_launch_filename_to_df_all_jobs()

        self._create_launch_strings()

        self._update_job_statuses_from_prior_exit0s()

        # For creation of a launch script external to the container
        # retail the list of all the .slurm or .bsub files getting launched.
        launch_filenames = []

        for subdir in ['PathProx', 'ddG_monomer', 'SequenceAnnotation', 'DigenicAnalysis']:
            if len(self._launch_strings[subdir]) == 0:
                # It's noteworth if we have a normal gene entry (not casewide) and a gene-related job is not running
                if (self._gene == 'casewide' and subdir == 'DigenicAnalysis') or (
                        self._gene != 'casewide' and subdir != 'DigenicAnalysis'):
                    LOGGER.info(
                        "No %s jobs will be run for %s %s %s", subdir, self._gene, self._refseq, self._mutation)
            else:
                launch_directory, launch_directory_outside_container = self._build_a_launch_dir(subdir)
                launch_stdout_directory = os.path.join(launch_directory, "stdout")

                old_umask = os.umask(2)
                os.makedirs(launch_stdout_directory, exist_ok=True)
                os.umask(old_umask)

                set_capra_group_sticky(launch_stdout_directory)

                launch_filename = os.path.join(launch_directory,
                                               "%s_%s.%s" % (self._geneRefseqMutation_OR_casewideString,
                                                             subdir,
                                                             self._slurm_or_bsub))
                LOGGER.info("Creating: %s", launch_filename)

                launch_filename_outside_container = launch_filename
                if launch_directory_outside_container and launch_directory_outside_container != launch_directory:
                    launch_filename_outside_container = os.path.join(
                        launch_directory_outside_container,
                        "%s_%s.%s" % (self._geneRefseqMutation_OR_casewideString,
                        subdir,
                        self._slurm_or_bsub))


                job_id = None

                # Repoint the -o output directory in case we are launching a container
                if launch_directory_outside_container:
                    launch_stdout_directory = os.path.join(launch_directory_outside_container, "stdout")

                launch_filenames.append(launch_filename_outside_container)

                if CLUSTER_TYPE == 'LSF':
                    # job_name, job_count = \
                    self._create_bsub_file(subdir, launch_filename, launch_stdout_directory)
                    if args.nolaunch:
                        LOGGER.warning("bsub will not be called to launch %s because args.nolaunch=%s",
                                       launch_filename,
                                       args.nolaunch)
                    else:
                        # bsub_command_list = ['bsub', '-i', launch_filename]
                        # Extension of job name to handle array now in .bsub file itself
                        # if job_count > 1:
                        #    bsub_command_list.extend(['-J', "%s[1-%d]" % (job_name, job_count)])
                        job_id = bsub_submit(launch_filename)
                else:
                    self._create_slurm_file(subdir, launch_filename, launch_stdout_directory)
                    if args.nolaunch:
                        LOGGER.warning("sbatch will not be called to launch %s because args.nolaunch=%s",
                                       launch_filename,
                                       args.nolaunch)
                    else:
                        job_id = slurm_submit(['sbatch', launch_filename])
                        # job_id = 999 # slurm_submit(['sbatch',slurm_file])

                if not args.nolaunch:
                    for uniquekey_launchstring in self._launch_strings[subdir]:
                        self._jobids[uniquekey_launchstring[0]] = job_id
                        self._jobinfo[uniquekey_launchstring[0]] = 'Submitted'
                        self._jobprogress[uniquekey_launchstring[0]] = 'In queue'
                        self._exitcodes[uniquekey_launchstring[0]] = None

            # Now, finally, record "Submitted" in all the status info files
            for status_dir in self._status_dir_list:
                with open('%s/info' % status_dir, 'w') as f:
                    f.write('Submitted')

        return self._df_updated_workstatus(), launch_filenames


# collaboration_log_dir = os.path.join(collaboration_dir,"log")
# if not os.path.exists(collaboration_log_dir):  # python 3 has exist_ok parameter...
#   LOGGER.error("%s not found.  It should have been created by psb_plan.py"%collaboration_log_dir)
#  sys.exit(1)

mutation_dir = ''


def launch_one_mutation(workplan_or_case_filename: str) -> Tuple[pd.DataFrame, str, List[str]]:
    # Load the schedule of jobs that was created by psb_plan.py
    df_all_jobs = pd.read_csv(workplan_or_case_filename, sep='\t', keep_default_na=False, na_filter=False)

    if len(df_all_jobs) < 1:
        LOGGER.warning("No rows in work plan file %s.  No jobs will be launched.", workplan_or_case_filename)
        return pd.DataFrame(), workplan_or_case_filename, []

    LOGGER.info("%d rows read from work plan file %s", len(df_all_jobs), workplan_or_case_filename)

    # Let's take care to not relaunch jobs that are already running, or complete, unless --relaunch flag is active
    workstatus_filename = workplan_or_case_filename.replace('workplan', 'workstatus')

    df_prior_exit0 = pd.DataFrame()

    try:
        # import pdb; pdb.set_trace()
        df_workstatus = pd.read_csv(workstatus_filename, sep='\t', keep_default_na=False, na_filter=False,
                                    dtype=str)

    except FileNotFoundError:
        # It's A-OK to not have a workstatus
        pass

    else:
        df_prior_exit0 = df_workstatus[df_workstatus['ExitCode'] == '0'].copy()
        df_prior_exit0.set_index('uniquekey', inplace=True)

    if args.relaunch and len(df_prior_exit0) > 0:
        LOGGER.info("--relaunch flag was passed, including %d with prior ExitCode=0", len(df_prior_exit0))

    jobs_launcher = JobsLauncher(workplan_or_case_filename, df_all_jobs, df_prior_exit0)
    old_umask = os.umask(2)
    os.makedirs(jobs_launcher.mutation_dir, exist_ok=True)
    os.umask(old_umask)

    set_capra_group_sticky(mutation_dir)

    # pw_name = pwd.getpwuid( os.getuid() ).pw_name # example jsheehaj or mothcw
    log_filename = os.path.join(mutation_dir, "psb_launch.log")

    if oneMutationOnly:
        sys.stderr.write("psb_launch log file is %s\n" % log_filename)
    else:
        LOGGER.info("additionally logging to file %s", log_filename)

    need_roll = os.path.isfile(log_filename)

    local_fh = RotatingFileHandler(log_filename, backupCount=7)
    formatter = logging.Formatter('%(asctime)s %(levelname)-4s [%(filename)20s:%(lineno)d] %(message)s',
                                  datefmt="%H:%M:%S")
    local_fh.setFormatter(formatter)
    local_fh.setLevel(logging.INFO)
    LOGGER.addHandler(local_fh)

    if need_roll:
        local_fh.doRollover()

    df_updated_workstatus,launch_filenames = jobs_launcher.doit()

    # Now create a new workstatus file
    print("Recording %d jobids in %s" % (len(df_updated_workstatus), workstatus_filename))

    df_updated_workstatus.to_csv(workstatus_filename, sep='\t')

    # DO NOT FORGET TO TEAR DOWN THE EXTRA LOGGER

    LOGGER.removeHandler(local_fh)
    # local_fh = None
    return df_updated_workstatus, workstatus_filename, launch_filenames


def main():
    global mutation_dir
    # Main logic here.  A period in the argument means the user wants to launch one mutation only,
    # directly from a single mutation output file of psb_plan.py


    # custom_launcher_filename = os.path.join(
    #     collaboration_dir,
    #     ('launch_%s.py' % args.projectORworkplan) if not oneMutationOnly else "launch_custom.py")

    custom_launcher_filename = ('launch_%s.py' % args.projectORworkplan) if not oneMutationOnly else "launch_custom.py"

    LOGGER.info("%s will be created as an external launch utility", custom_launcher_filename)

    # We create a super-simple (few dependencies) launcher program that can be run _outside_ the container
    launcher_f = open(custom_launcher_filename,'w')

    # Make the customer launcher .py file executable by user and group
    mode = os.fstat(launcher_f.fileno()).st_mode
    mode |= stat.S_IXUSR | stat.S_IXGRP
    os.fchmod(launcher_f.fileno(), stat.S_IMODE(mode))

    # We care not whether python 2 or python 3 in this small app....
    launcher_f.write("""\
#!/usr/bin/env python3\n""")

    # Add a docstring to the launcher program
    launcher_f.write('"""\n')
    launcher_f.write("%s: Launcher utility for PSB Pipeline\n" % custom_launcher_filename)
    launcher_f.write("""\
Created %s by command line:
%s
""" % (datetime.now(), ' '.join(sys.argv)))

    launcher_f.write('\n"""\n')

    launcher_f.write("""\

import subprocess as sp
import sys
import re
import csv
import logging

logging.basicConfig(level=logging.INFO)
""")

    bsub_submit_function_source_code = inspect.getsource(bsub_submit)

    launcher_f.write("""\

# bsub_submit() function copied from bsub.py
""")

    launcher_f.write(bsub_submit_function_source_code)

    launcher_f.write("""\

# slurm_submit() function copied from slurm.py
""")

    slurm_submit_function_source_code = inspect.getsource(slurm_submit)
    launcher_f.write(slurm_submit_function_source_code)
    launcher_f.write("\n")

    launcher_f.write("""\
# Read the workstatus.csv file for a variant
# And load it into a list of dictionaries
# one per row.  The Dictionary keys for each row are the 
# row1 header elements
def read_workstatus_csv(workstatus_filename):
    workstatus_rows = []
    with open(workstatus_filename) as f:
        reader = csv.DictReader(f,delimiter='\\t')
        for csv_row in reader:
            workstatus_rows.append(csv_row)
    return workstatus_rows

# Write the workstatus rows, with the 'jobid' now populated    
def write_workstatus_csv(workstatus_filename,workstatus_rows):
    with open(workstatus_filename,'w') as f:
        writer = csv.DictWriter(f,fieldnames=workstatus_rows[0].keys(),delimiter='\\t')
        writer.writeheader()
        for csv_row in workstatus_rows:
            writer.writerow(csv_row)

# For every workstatus file row that is launched by the sbatch/bsub file
# Update the jobid
def add_jobid_to_workstatus(workstatus_csv, launch_filename, job_id):
    for csv_row in workstatus_csv:
        if 'launch_file' in  csv_row and csv_row['launch_file'] == launch_filename:
            csv_row['jobid'] = job_id
""")

    if oneMutationOnly:
        # The argument is a complete workplan filename
        df_updated_workstatus, workstatus_filename, launch_filenames = launch_one_mutation(args.projectORworkplan)
    else:
        # The argument is an entire project UDN124356
        udn_csv_filename = os.path.join(collaboration_dir,
                                        "%s_missense.csv" % args.projectORworkplan)
        print("Retrieving project mutations from %s" % udn_csv_filename)
        df_all_mutations = pd.read_csv(udn_csv_filename, sep=',', index_col=None,
                                       keep_default_na=False, encoding='utf8',
                                       comment='#', skipinitialspace=True)
        df_all_mutations.fillna('NA', inplace=True)
        print("Launching all jobs for %d mutations" % len(df_all_mutations))

        # For each variant in the case, create workstatus files, and gather job ids if
        # launching now.  Build up the custom launcher launcher_f script for a later launch
        for index, row in df_all_mutations.iterrows():
            print("Launching %-10s %-10s %-6s" % (row['gene'], row['refseq'], row['mutation']))
            if 'GeneOnly' in row['refseq']:
                row['refseq'] = 'NA'
            mutation_dir = os.path.join(collaboration_dir, "%s_%s_%s" % (row['gene'], row['refseq'], row['mutation']))
            if not os.path.exists(mutation_dir):  # python 3 has exist_ok parameter...
                LOGGER.critical(
                    "The specific mutation directory %s should have been created by psb_plan.py.  Fatal problem.",
                    mutation_dir)
                sys.exit(1)

            workplan_filename = "%s/%s_%s_%s_workplan.csv" % (mutation_dir, row['gene'], row['refseq'], row['mutation'])
            # print workplan_filename
            # The argument is a complete workplan filename
            df_updated_workstatus, workstatus_filename, launch_filenames = launch_one_mutation(workplan_filename)

            # Very cheesy - but the idea is that we need the external launcher to use the mounted
            # Directory - NOT the container-visible directory:
            if mounted_collaboration_directory and mounted_collaboration_directory != collaboration_dir:
                workstatus_filename = workstatus_filename.replace(collaboration_dir, mounted_collaboration_directory, 1)

            launcher_f.write("""

# Launching all jobs for %s
""" % "%s_%s_%s" % (row['gene'], row['refseq'], row['mutation']))

            launcher_f.write("workstatus_csv = read_workstatus_csv('%s')\n" % workstatus_filename)
            launcher_f.write("launch_filenames = [\n    '%s']\n" % "',\n    '".join(launch_filenames))
            launcher_f.write("for launch_filename in launch_filenames:\n")

            # Write code to call bsub_submit or slurm_submit included in the .py above
            if CLUSTER_TYPE == 'LSF':
                launcher_f.write("    job_id = bsub_submit(launch_filename)\n")
            else:
                launcher_f.write("    job_id = slurm_submit(['sbatch',launch_filename])\n")

            launcher_f.write("    add_jobid_to_workstatus(workstatus_csv, launch_filename, job_id)\n")

            launcher_f.write("write_workstatus_csv('%s',workstatus_csv)\n" % workstatus_filename)

        # The casewide launch is a separate case from the above
        casewide_workplan_filename = "casewide/casewide_workplan.csv"
        if os.path.exists(casewide_workplan_filename):
            print("Launching casewide job(s) listed in %s" % casewide_workplan_filename)
            df_updated_workstatus, workstatus_filename, launch_filenames = launch_one_mutation(
                casewide_workplan_filename)

            # Very cheesy - but the idea is that we need the external launcher to use the mounted
            # Directory - NOT the container-visible directory:
            if mounted_collaboration_directory and mounted_collaboration_directory != collaboration_dir:
                workstatus_filename = workstatus_filename.replace(collaboration_dir, mounted_collaboration_directory, 1)

            launcher_f.write("""

# Launching casewide jobs
""")
            launcher_f.write("workstatus_csv = read_workstatus_csv('%s')\n" % workstatus_filename)
            launcher_f.write("launch_filenames = [\n    '%s']\n" % "',\n    '".join(launch_filenames))
            launcher_f.write("for launch_filename in launch_filenames:\n")
            # Write code to call bsub_submit or slurm_submit included in the .py above
            if CLUSTER_TYPE == 'LSF':
                launcher_f.write("    job_id = bsub_submit(launch_filename)\n")
            else:
                launcher_f.write("    job_id = slurm_submit(['sbatch',launch_filename])\n")

            launcher_f.write("    add_jobid_to_workstatus(workstatus_csv, launch_filename, job_id)\n")

            launcher_f.write("write_workstatus_csv('%s',workstatus_csv)\n" % workstatus_filename)
        else:
            print("No casewide file (%s) was found.  No casewide jobs will be started" % casewide_workplan_filename)

    user_launch_message = None
    if CONTAINER_TYPE:
        user_launch_message = "./%s should be run outside %s container to launch jobs" % (
            custom_launcher_filename, CONTAINER_TYPE)
    elif args.nolaunch:
        user_launch_message = "Run ./%s to launch jobs" % custom_launcher_filename

    if user_launch_message:
        LOGGER.info(user_launch_message)
        print("\n\n" + user_launch_message)


if __name__ == '__main__':
    main()
sys.exit(0)
