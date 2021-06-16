#!/usr/bin/env python
"""
class DDG_base underpins both DDG_monomer and DDG_cartesian with
a variety of methods shared by bost.
"""
import os
import time
import datetime
import subprocess
import lzma
import pprint
from io import StringIO
from psb_shared.ddg_repo import DDG_repo
from typing import List, Tuple, Union

import logging

LOGGER = logging.getLogger(__name__)

class DDG_base(object):
    def refresh_ddg_repo_and_mutations(self, ddg_repo: DDG_repo, mutations: Union[str, List[str]]):
        """
        When retrieving results from the repository, one should be able to update the ddg_repo
        pointer to a new variant without repeating all the other checks
        """

        self._ddg_repo = ddg_repo

        if type(mutations) == str:
            self._mutations = [mutations]  # Nasty to switch out like this - watch closely
        else:
            assert type(self._mutations) == list
            # Sort the mutations by the integer part of their pdb code
            # and the PDB insertion code if mutation[-2] isalpha(!)
            self._mutations = sorted(mutations, key=lambda mutation: (
                int(mutation[1:-2 if mutation[-2].isalpha() else -1]),
                mutation[-2] if mutation[-2].isalpha() else ' '
            ))

        # Not sure we really want this potentially long string in filenames
        # Keep thinking about it.
        self._mutationtext = ",".join(self._mutations)

    def __init__(self, ddg_repo: DDG_repo, mutations: Union[str, List[str]]):
        """
        Initialize self._ddg_repo, self._mutations from caller's data.

        Generally called from derived DDG_monomer or DDG_cartesian function to handle shared
        initializtion needs.
        :param ddg_repo:
        :param mutations: S123AF format (or list of these) that indicate AA, res#, insert code, new AA.
        """

        self._mutations = None
        self._ddg_repo = None
        self._mutationtext = None
        self.refresh_ddg_repo_and_mutations(ddg_repo, mutations)

        self._to_pdb = []

        LOGGER.info("ddg_base.__init__ completed after setting:\n%s" % pprint.pformat(self.__dict__))

        # timestamp_suffix = datetime.datetime.now().strftime("%d%m%y%H%M%S")
        # isoresult, isopdbname, to_pdb = isolate_chain(self._pdb, self._chain)
        # if not isoresult:
        #    return self.failure(None, to_pdb)
        # else:
        #    self._to_pdb = to_pdb

    def mutations_to_rosetta_poses(self):
        """
        Give the list of self._mutations in A123{I}B form where I is optional insertion code,
        convert these to their rosetta-numbered 1..N formats in 3 ways:
            rosetta_mutations: ['A123X',...]  List of Rosetta-reNumbered AA res# AA entries
            mutation_resids: [(' ',123,'I'),... ] biopython residue IDs not renumbered
            rosetta_residue_number_to_residue_xref: A cross-reference of the rosetta renumberings to original residues
        """
        # Convert mutation to pose numbering for internal consistency
        rosetta_mutations = []  # List of mutations input to ddg

        mutation_resids = []
        # Capture the inner bytes,123A of the S123AW format, as Biopython resid tuple
        for mutation in self._mutations:
            if mutation[-2].isalpha():  # Then Mutation is like ex: S123AF (Ser->Phe at res 123A)
                mutation_resid = (' ', int(mutation[1:-2]), mutation[-2])
            else:  # The much more common case of format S456T
                mutation_resid = (' ', int(mutation[1:-1]), ' ')
            assert mutation_resid in self._ddg_repo.residue_to_clean_xref, (
                    "mutation %s is not in source structure" % str(mutation_resid))
            mutation_resids.append(mutation_resid)

        for mutation, mutation_resid in zip(self._mutations, mutation_resids):
            rosetta_residue_number = self._ddg_repo.residue_to_clean_xref[mutation_resid][1]
            LOGGER.info("Res %s converted to rosetta cleaned res # %d" % (mutation_resid, rosetta_residue_number))
            # Add to list of rosetta mutations in form A123S
            rosetta_mutations.append(mutation[0] + str(rosetta_residue_number) + mutation[-1])

        if not rosetta_mutations:
            self._ddg_repo.psb_status_manager.sys_exit_failure("No mutations found in parameter self._mutations %s" % self._mutations)

        if len(rosetta_mutations) != 1:
            self._ddg_repo.psb_status_manager.sys_exit_failure("Only 1 mutation can be setisfied in parameter self._mutations %s" % self._mutations)

        # Create a reverse lookup dictionary, mapping the 1..N rosetta resdiue numbers
        # back to source structure residue IDs (which might have insertion codes)
        rosetta_residue_number_to_residue_xref = {}
        for residue, clean in self._ddg_repo.residue_to_clean_xref.items():
            try:
                rosetta_residue_number = int(clean[1])
                rosetta_residue_number_to_residue_xref[rosetta_residue_number] = residue
                LOGGER.debug("%s %s"%(rosetta_residue_number,residue))
            except:
                self._ddg_repo.psb_status_manager.sys_exit_failure(
                    "Creating rosetta->residue reverse xref Failure to deal with residue=%s clean=%s"%(residue,clean))

        return mutation_resids,rosetta_mutations,rosetta_residue_number_to_residue_xref


    @staticmethod
    def evaluate_swiss(swiss_modelid, swiss_remark3_metrics):
        """Swiss models are assumed to be of reasonable fundamental quality.
           We just  check for sequence id >= 40 and move on

           :param swiss_modelid:         Unique swiss model identifier
           :param swiss_remark3_metrics:  Dictionary preloaded via PDBMapSwiss call

           :return True/False for OK to continue or NOT
        """
        sid_threshold = 30.0
        if 'sid' in swiss_remark3_metrics:
            swiss_sid = float(swiss_remark3_metrics['sid'])
            sequence_identity_acceptable = (swiss_sid >= sid_threshold)
            if sequence_identity_acceptable:
                LOGGER.info("Swiss Model %s seq identity=%f  >= Minimum Acceptable: %f" % (
                    swiss_modelid, swiss_sid, sid_threshold))
            else:
                LOGGER.warning("Swiss Model %s seq identity=%f  <= Minimum Acceptable: %f" % (
                    swiss_modelid, swiss_sid, sid_threshold))
            return sequence_identity_acceptable

        LOGGER.warning("Swiss Model %s lacks REMARK for sid (seq identity)" % swiss_modelid)
        return False

    @staticmethod
    def evaluate_modbase(modbase_fullpath_or_fin):
        """Modbase models are evaluated for whether or not they are of high enough quality
        To run in ddg monomer.
        Comuted quality metrics must all be sufficient, including
        modpipe: modpipe quality score (cumulative overall score from modbase pipeline)
        ztsvmod: RMSD of model to template.
        seqid: sequence id of model seq and template seq

        param: modbase_fullpath_or_fin  Full pathname of modbase.gz file or a StringIO

        """

        modpipe_acceptable = True  # Removed filter for low modpipe scores (seqid and rmsd to template should be enough)
        tsvmod_acceptable = False
        sequence_identity_acceptable = False

        sid_threshold = 30.0
        rmsd_threshold = 4.0
        qual_threshold = 1.1

        LOGGER.info('Checking quality of %s', modbase_fullpath_or_fin)
        with modbase_fullpath_or_fin if type(modbase_fullpath_or_fin) == StringIO else \
                lzma.open(modbase_fullpath_or_fin, 'rt') as infile:
            for line in infile:
                if line.startswith('REMARK 220 SEQUENCE IDENTITY'):
                    try:
                        modbase_sid = float(line.strip().split()[4])
                    except IndexError:  # Sometimes this line is entirely blank in modbase
                        modbase_sid = 0.0
                    except TypeError:
                        modbase_sid = 0.0
                    except ValueError:
                        modbase_sid = 0.0
                    if modbase_sid >= sid_threshold:
                        sequence_identity_acceptable = True
                        LOGGER.info("Modbase Sequence identity of %0.1f acceptable as > threshold %0.1f", modbase_sid,
                                    sid_threshold)
                    else:
                        sequence_identity_acceptable = False
                        LOGGER.warning("Modbase Sequence Identity: %0.1f < threshold %0.1f", modbase_sid, sid_threshold)

                elif line.startswith('REMARK 220 TSVMOD RMSD'):
                    try:
                        rmsd = float(line.strip().split()[4])
                    except IndexError:  # Sometimes this line is entirely blank in modebase
                        rmsd = 1000.0
                    except TypeError:
                        rmsd = 1000.0
                    except ValueError:
                        rmsd = 1000.0
                    tsvmod_acceptable = rmsd <= rmsd_threshold
                    if tsvmod_acceptable:
                        LOGGER.info("TSVMOD RMSD: %0.1f OK ( <= threshold %0.1f)", rmsd, rmsd_threshold)
                    else:
                        LOGGER.info("TSVMOD RMSD: %0.1f TOO HIGH, threshold %0.1f", rmsd, rmsd_threshold)
                elif line.startswith('REMARK 220 MODPIPE QUALITY SCORE'):
                    try:
                        qual = float(line.strip().split()[5])
                    except IndexError:  # Sometimes this line is entirely blank in modbase
                        qual = 0.0
                    except TypeError:
                        qual = 0.0
                    except ValueError:
                        qual = 0.0

                    if qual >= 1.1:
                        modpipe = True  # Well - does not do much given low filter removal ...
                    LOGGER.info("Modpipe quality score (not checked) %0.2f" % qual)

        modbase_ok = modpipe_acceptable and tsvmod_acceptable and sequence_identity_acceptable

        LOGGER.log(logging.INFO if modbase_ok else logging.WARNING,
                   "Modbase modpipe: %s, tsvmod: %s, template identity: %s" % (
                       str(modpipe_acceptable), str(tsvmod_acceptable), str(sequence_identity_acceptable)))

        return modbase_ok

    @staticmethod
    def _move_prior_file_if_exists(filename):
        if os.path.exists(filename):
            archive_dir = "./archive"
            os.makedirs(archive_dir, exist_ok=True)
            # append the modification time of the file to the name
            # and archive in ./archive subdir
            archive_filename = os.path.join(
                archive_dir,
                # precede the filename with time stamp of the file's last modification.
                datetime.datetime.fromtimestamp(os.path.getmtime(filename)).strftime(
                    "%Y%M%d%H%M%S") + "_%s" % filename
            )
            try:
                os.replace(filename, archive_filename)
                LOGGER.info("Saved prior run %s as %s" % (filename, archive_filename))
            except OSError:
                LOGGER.exception("Failed to save prior run %s as %s: " % (filename, archive_filename))

    @staticmethod
    def _command_ran_previously(binary_program_basename: str) -> Tuple[int, str, str]:
        """
        If this program ran earlier, don't re-run it.
        Instead return the results of the last run

        param: binary_program_basename:  The binary program file to check for previous completion
        :return (
            exit_code: integer exit code we can test for 0
            str: stdout from command
            str: stderr from command
            )

        """
        stdout_filename, stderr_filename, exitcd_filename = DDG_base._command_result_filenames(binary_program_basename)
        previous_exit, previous_exit_int = DDG_base._get_previous_exit(exitcd_filename)

        previous_stdout = ""
        if os.path.isfile(stdout_filename):
            with open(stdout_filename, 'r') as stdout_record:
                previous_stdout = stdout_record.read()  # This stdout from previous successful run

        previous_stderr = ""
        if os.path.isfile(stderr_filename):
            with open(stderr_filename, 'r') as stderr_record:
                previous_stderr = stderr_record.read()  # Read stderr from previous successful run

        if previous_exit_int == 0:
            LOGGER.info("Previous successful run found\n.  NOT re-running:\n%s" % (
                binary_program_basename))
        else:
            LOGGER.warning("Previous run did not exit with code 0\n.  Re-running:\n%s" % (
                binary_program_basename))

        return (previous_exit_int,
                previous_stdout,
                previous_stderr)  # runtime 0 may or may not be quite right

    def _run_command_line_terminate_on_nonzero_exit(self,
                                                    command_line_list: List[str],
                                                   additional_files_to_archive: List[str] = [],
                                                   force_rerun=False) -> Tuple[int, str, str, float]:
        """
        Launch a shell to run the command line.

        Definitely run if force_rerun==True (unusual)

        Typically, we do not re-run if we find a .exit file with a zero exit value
        In that case, we return the saved .stdout and .stderr files, and save
        the cpu cycles

        Terminates hard if non-zero return value from command line

        :param command_line_list: list of command line components

        :return (
            int: return/exit code from command invocation
            str: stdout from command
            str: stderr from command
            float: run time

        """

        # Example: minimzer.linuxgccrelease
        binary_program_basename = os.path.basename(command_line_list[0])
        stdout_filename, stderr_filename, exitcd_filename = DDG_base._command_result_filenames(binary_program_basename)

        # We've been playing with letting the command list have embedded spaces which have to be
        # split back apart...
        # We may get bitten at some point because some of our arguments are actually 2 arguments
        # separated by a space.  For now...
        submitted_command_line_list = []
        for arg_with_spaces in command_line_list:
            submitted_command_line_list.extend(arg_with_spaces.split(' '))

        # Archive any additional output files once we commit to re-launching
        for additional_file in additional_files_to_archive:
            DDG_base._move_prior_file_if_exists(additional_file)

            # If there is already a file out there with exit code 0, then return the old status

        # For sanity, capture precicsely a reproducible shell command at work for future analysis
        commandline_record_filename = binary_program_basename + ".commandline_record.sh"
        DDG_base._move_prior_file_if_exists(commandline_record_filename)
        with open(commandline_record_filename, 'w') as commandline_record:
            commandline_record.write("#!/bin/bash\n")
            commandline_record.write("# Record of command invocation\n")
            commandline_record.write("# %s\n" % time.ctime(time.time()))
            commandline_record.write(" \\\n".join(command_line_list))
            # End the command with redirects to the stderr/stdout files
            commandline_record.write(" \\\n> %s 2> %s\n" % (stdout_filename, stderr_filename))
            commandline_record.write("echo Command exited with $?\n")

        command_start_time = time.perf_counter()
        LOGGER.info("Running below command line from cwd=%s:\n%s" % (os.getcwd(),' '.join(submitted_command_line_list)))
        # run using new API call,
        completed_process = subprocess.run(submitted_command_line_list, shell=False, text=True, capture_output=True)
        fail_message = None
        if completed_process.returncode == 0:
            LOGGER.info("%s completed successfully (exit 0)", binary_program_basename)
        else:
            fail_message = "%s failed with exit %d" % (binary_program_basename, completed_process.returncode)
            LOGGER.critical(fail_message)

        # Record all stdout and stderr and exit code
        DDG_base._move_prior_file_if_exists(stdout_filename)
        with open(stdout_filename, 'w') as stdout_record:
            stdout_record.write(completed_process.stdout)  # This writes the 'str'

        DDG_base._move_prior_file_if_exists(stderr_filename)
        with open(stderr_filename, 'w') as stderr_record:
            stderr_record.write(completed_process.stderr)

        DDG_base._move_prior_file_if_exists(exitcd_filename)
        with open(exitcd_filename, 'w') as exitcd_record:
            exitcd_record.write(str(completed_process.returncode))
        LOGGER.info(
            "stdout/stderr/exit saved in files %s/%s/%s" % (stdout_filename, stderr_filename, exitcd_filename))

        if fail_message:
            self._ddg_repo.psb_status_manager.sys_exit_failure(fail_message)

        # Convert output byte streams to manageable unicode strings.
        # Return (stdout,stderr,elapsed_time)
        return (completed_process.returncode,
                completed_process.stdout,
                completed_process.stderr,
                time.perf_counter() - command_start_time)


    def final_results_filename(self) -> str:
        return "%s_%s_%s.csv" % (
            self._ddg_repo.structure_id,
            self._ddg_repo.chain_id,
            self._mutationtext
        )

    @staticmethod
    def _command_result_filenames(binary_program_basename: str):
        return (
            binary_program_basename + ".stdout",
            binary_program_basename + ".stderr",
            binary_program_basename + ".exit")

    @staticmethod
    def _get_previous_exit(exitcd_filename: str) -> Tuple[str, int]:
        previous_exit = None
        if os.path.isfile(exitcd_filename):
            with open(exitcd_filename, 'r') as exitcd_record:
                previous_exit = exitcd_record.read()
                if previous_exit:
                    previous_exit = previous_exit.split()[0]

        previous_exit_int = 1
        if previous_exit:
            try:
                previous_exit_int = int(previous_exit)
            except TypeError:
                pass
            except ValueError:
                pass

        return previous_exit, previous_exit_int

