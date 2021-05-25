#!/usr/bin/env python
"""
class DDG_monomer manages all aspects of Rosetta ddg_monomer calculations 
inside the directory heirarchy supplied by a ddg_repo object
"""
import os
import datetime
import gzip
import lzma
import pprint
import subprocess
import tempfile
import time
import pandas as pd
from io import StringIO
from psb_shared.ddg_repo import DDG_repo
from typing import Dict, List, Tuple, Union
# from sys import argv, stderr, stdout
# from os import popen, system
# from os.path import exists, basename
from psb_shared.psb_progress import PsbStatusManager

import logging

LOGGER = logging.getLogger(__name__)


class DDG_monomer(object):

    def refresh_ddg_repo_and_mutations(self,ddg_repo: DDG_repo, mutations: Union[str, List[str]]):
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


    def _verify_applications_available(self):
        for application in [
                self._minimize_application_filename,
                self._per_residues_application_filename,
                self._ddg_monomer_application_filename,
                self._score_jd2_application_filename]:
            application_fullpath = os.path.join(self._ddg_repo.rosetta_bin_dir, application)
            assert os.path.exists(application_fullpath), (
                    "%s is a required application for ddg_monomer, but not found" % application_fullpath)

    def __init__(self, ddg_repo: DDG_repo, mutations: Union[str, List[str]]):
        """
        :param ddg_repo:
        :param mutations: S123AF format (or list of these) that indicate AA, res#, insert code, new AA.
        """

        self.refresh_ddg_repo_and_mutations(ddg_repo, mutations)

        #        self._root = root_path
        # Set the binary filenames in one place and make sure we have them before we start
        self._minimize_application_filename =    'minimize_with_cst.default.linuxgccrelease'
        self._per_residues_application_filename ='per_residue_energies.linuxgccrelease'
        self._ddg_monomer_application_filename = 'ddg_monomer.linuxgccrelease'
        self._score_jd2_application_filename =   'score_jd2.linuxgccrelease'

        self._application_timers = ['minimize', 'rescore', 'ddg_monomer']
        self._to_pdb = []

        LOGGER.info("ddg_monomer.__init__ completed after setting:\n%s" % pprint.pformat(self.__dict__))

    @staticmethod
    def evaluate_swiss(swiss_modelid,swiss_remark3_metrics):
        """Swiss models are assumed to be of reasonable fundamental quality.
           We just  check for sequence id >= 40 and move on

           :param swiss_model_id:         Unique swiss model identifier
           :param swiss_remark3_metrics:  Dictionary preloaded via PDBMapSwiss call
           
           :return True/False for OK to continue or NOT
        """
        sid_threshold = 30.0
        if 'sid' in swiss_remark3_metrics:
            swiss_sid = float(swiss_remark3_metrics['sid'])
            sequence_identity_acceptable = (swiss_sid >= sid_threshold)
            if sequence_identity_acceptable:
                LOGGER.info("Swiss Model %s seq identity=%f  >= Minimum Acceptable: %f"%(
                    swiss_modelid,swiss_sid,sid_threshold))
            else:
                LOGGER.warning("Swiss Model %s seq identity=%f  <= Minimum Acceptable: %f"%(
                    swiss_modelid,swiss_sid,sid_threshold))
            return sequence_identity_acceptable

        LOGGER.warning("Swiss Model %s lacks REMARK for sid (seq identity)"%modelid)
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

        modpipe_acceptable = True #Removed filter for low modpipe scores (seqid and rmsd to template should be enough)
        tsvmod_acceptable  = False
        sequence_identity_acceptable = False

        sid_threshold = 30.0
        rmsd_threshold = 4.0
        qual_threshold = 1.1

        LOGGER.info('Checking quality of %s',modbase_fullpath_or_fin)
        with modbase_fullpath_or_fin if type(modbase_fullpath_or_fin) == StringIO else \
             lzma.open(modbase_fullpath_or_fin,'rt') as infile:
            for line in infile:
                if line.startswith('REMARK 220 SEQUENCE IDENTITY'):
                    try:
                        modbase_sid = float(line.strip().split()[4])
                    except IndexError: # Sometimes this line is entirely blank in modbase
                        modbase_sid = 0.0
                    except TypeError:
                        modbase_sid = 0.0
                    except ValueError:
                        modbase_sid = 0.0
                    if modbase_sid>=sid_threshold:
                        sequence_identity_acceptable = True
                        LOGGER.info("Modbase Sequence identity of %0.1f acceptable as > threshold %0.1f",modbase_sid,sid_threshold)
                    else:
                        sequence_identity_acceptable = False
                        LOGGER.warning("Modbase Sequence Identity: %0.1f < threshold %0.1f",modbase_sid,sid_threshold)

                elif line.startswith('REMARK 220 TSVMOD RMSD'):
                    try:
                        rmsd = float(line.strip().split()[4])
                    except IndexError: # Sometimes this line is entirely blank in modebase
                        rmsd = 1000.0
                    except TypeError:
                        rmsd = 1000.0
                    except ValueError:
                        rmsd = 1000.0
                    tsvmod_acceptable = rmsd<=rmsd_threshold
                    if tsvmod_acceptable:
                        LOGGER.info("TSVMOD RMSD: %0.1f OK ( <= threshold %0.1f)",rmsd,rmsd_threshold)
                    else:
                        LOGGER.info("TSVMOD RMSD: %0.1f TOO HIGH, threshold %0.1f",rmsd,rmsd_threshold)
                elif line.startswith('REMARK 220 MODPIPE QUALITY SCORE'):
                    try:
                        qual = float(line.strip().split()[5])
                    except IndexError: # Sometimes this line is entirely blank in modbase
                        qual = 0.0
                    except TypeError:
                        qual = 0.0
                    except ValueError:
                        qual = 0.0

                    if qual>=1.1:
                        modpipe = True # Well - does not do much given low filter removal ...
                    LOGGER.info("Modpipe quality score (not checked) %0.2f"%qual)

        modbase_ok = modpipe_acceptable and tsvmod_acceptable and sequence_identity_acceptable

        LOGGER.log(logging.INFO if modbase_ok else logging.WARNING, "Modbase modpipe: %s, tsvmod: %s, template identity: %s"%(
            str(modpipe_acceptable),str(tsvmod_acceptable),str(sequence_identity_acceptable)))

        return modbase_ok

    def final_results_filename(self) -> str:
        return "%s_%s_%s.csv" % (
            self._ddg_repo.structure_id,
            self._ddg_repo.chain_id,
            self._mutationtext
        )

    def _command_result_filenames(self,binary_program_basename:str):
        return (
            binary_program_basename + ".stdout",
            binary_program_basename + ".stderr",
            binary_program_basename + ".exit")

    def _get_previous_exit(self,exitcd_filename:str) -> Tuple[str,int]:
        previous_exit = None
        if os.path.isfile(exitcd_filename):
            with open(exitcd_filename,'r') as exitcd_record:
                previous_exit=exitcd_record.read()
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

        return previous_exit,previous_exit_int
        



    #        with open(self._pdb) as infile:
    #            self._pdblines = infile.readlines()

    # def log_success(self, times):
    #   mutationtext = ",".join("".join(str(x) for x in y) for y in self._mutations)
    #   nres = str(len(self._to_pdb) - 1)
    #   curdate = datetime.datetime.now().strftime("%d_%m_%y")
    #   curtime = datetime.datetime.now().strftime("%H_%M")
    #   wd = os.getcwd()
    #   if not os.path.isfile(self._success):
    #       with open(self._success, 'w') as outfile:
    #           header = ['date', 'time', 'pdbfile', 'chain', 'nres', 'mutation'] + self._timers + ['path']
    #           outfile.write("\t".join(header))
    #           outfile.write("\n")
    #   with open(self._success, 'a') as outfile:
    #       line = [curdate, curtime, self._pdb, self._chain, nres, mutationtext] + \
    #              ["%.2f" % (times[x] / 60.0) if x in times else 'NaN' for x in self._timers] + [wd]
    #       outfile.write("\t".join(line))
    #       outfile.write("\n")

    # def log_failure(self, message):
    #   wd = os.getcwd()
    #   mutationtext = ",".join("".join(str(x) for x in y) for y in self._mutations)
    #   curdate = datetime.datetime.now().strftime("%d_%m_%y")
    #   curtime = datetime.datetime.now().strftime("%H_%M")
    #   if not os.path.isfile(self._fail):
    #       with open(self._fail, 'w') as outfile:
    #           header = ['date', 'time', 'pdbfile', 'chain', 'mutation', 'message', 'path']
    #           outfile.write("\t".join(header))
    #           outfile.write("\n")
    #   with open(self._fail, 'a') as outfile:
    #       line = [curdate, curtime, self._pdb, self._chain, mutationtext, message, wd]
    #       outfile.write("\t".join(line))
    #       outfile.write("\n")

    # def failure(self, dumps, message):
    #   LOGGER.warning(message)
    #   if dumps is not None:
    #       for item in dumps:
    #           LOGGER.warning("Dumping %s to %s", item[0], item[1])
    #           with open(item[1], 'w') as outfile:
    #               outfile.write(item[2])
    #   # self.log_failure(message)
    #   return [False, message]

    """
    Old code that seems no longer used "non-standard ddg"
    def get_ddg(self, quality, silent=False):
        timestamp_suffix = datetime.datetime.now().strftime("%d%m%y%H%M%S")

        with open("ddg_predictions.out") as infile:
            muts = [x.strip().split()[1] for x in infile.readlines()[1:] if len(x.strip().split()) > 1]
        if not os.path.exists('top_models'):
            os.makedirs('top_models')
        top_three = {'wt': []}

        if not silent:
            wt_files = glob.glob("repacked*.pdb")
            mut_files = glob.glob("mut*.pdb")
            modellist = "models_%s.ls" % timestamp_suffix
            rescorefile = "ddg_%s.sc" % timestamp_suffix
            with open(modellist, 'w') as outfile:
                for x in wt_files:
                    outfile.write("\n".join(wt_files + mut_files))
                    outfile.write("\n")
            command = "score_jd2.linuxgccrelease -l %s -score:weights talaris2014 -out:file:scorefile %s" % (
            modellist, rescorefile)
            LOGGER.info("Scoring all ddg models with command:\n%s/%s", self._rosetta, command)
            runscore = subprocess.Popen("%s/%s" % (self._rosetta, command), shell=True, stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
            (scout, scerr) = runscore.communicate()
            scout = out.decode('latin')
            scerr = scerr.decode('latin')
            try:
                assert os.path.isfile(rescorefile)
            except AssertionError:
                scoutfile = "modscorefail%s.out" % timestamp_suffix
                scerrfile = "modscorefail%s.err" % timestamp_suffix
                return self.failure([["rescore output", scoutfile, scout], ["rescore stderr", scerrfile, scerr]],
                                    "Error scoring ddg models, scorefile is missing.")

            with open(rescorefile) as infile:
                raw_rescore = [x.strip().split() for x in infile.readlines() if x.strip() != '']

            try:
                assert len(raw_rescore) > 0
            except AssertionError:
                scoutfile = "modscorefail%s.out" % timestamp_suffix
                scerrfile = "modscorefail%s.err" % timestamp_suffix
                return self.failure([["rescore output", scoutfile, scout], ["rescore stderr", scerrfile, scerr]],
                                    "Error scoring ddg models, scorefile is empty.")
        else:
            raw_rescore = []
            with open("wt_.out") as infile:
                raw_rescore += [x.strip().split() for x in infile.readlines() if
                                x.strip().split()[-1] != 'description' and x.strip().split()[0] == 'SCORE:']
            for item in muts:
                with open("mut_%s.out" % item) as infile:
                    raw_rescore += [x.strip().split() for x in infile.readlines() if
                                    x.strip().split()[-1] != 'description' and x.strip().split()[0] == 'SCORE:']

        if not silent:
            try:
                score_ix = raw_rescore[0].index('total_score')
                name_ix = raw_rescore[0].index('description')
                headend = 0
            except ValueError:
                score_ix = raw_rescore[1].index('total_score')
                name_ix = raw_rescore[1].index('description')
                headend = 1
        else:
            score_ix = 1
            name_ix = -1
            headend = -1

        for line in raw_rescore[headend + 1:]:
            if silent:
                curmodel = line[name_ix]
            else:
                curmodel = "%s.pdb" % line[name_ix][:-5]
            curscore = float(line[score_ix])
            curkey = None
            if curmodel[:11] == 'repacked_wt':
                curkey = 'wt'
            else:
                for mut in muts:
                    if curmodel.startswith("mut_%s" % mut):
                        curkey = mut
            if curkey is None:
                continue
            if curkey not in top_three:
                top_three[curkey] = [[curscore, curmodel]]
            elif len(top_three[curkey]) < 3:
                top_three[curkey].append([curscore, curmodel])
            elif curscore < sorted(top_three[curkey])[-1][0]:
                LOGGER.info(sorted(top_three[curkey])[-1][0])
                top_three[curkey][-1] == [curscore, curmodel]
        allmeans = {}
        for keys in top_three:
            topmod = sorted(top_three[keys])[0][1]
            if silent:
                if keys == 'wt':
                    outfile = "wt_.out"
                else:
                    outfile = "mut_%s.out" % keys
                command = "extract_pdbs.linuxgccrelease -in:file:silent %s -in:file:tags %s" % (outfile, topmod)
                topmod += ".pdb"
                extcom = subprocess.Popen("%s/%s" % (self._rosetta, command), shell=True, stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE)
                (cstout, csterr) = extcom.communicate()
                cstout = cstout.decode('latin')
                csterr = csterr.decode('latin')
            copyfile(topmod, 'top_models/%s' % topmod)
            curmean = mean([x[0] for x in top_three[keys]])
            allmeans[keys] = curmean
        scores = []
        for keys in allmeans:
            if keys == 'wt': continue
            scores.append("%s\t%.3f" % (keys, (allmeans[keys] - allmeans['wt'])))
        return [True, scores]
    """

    def run(self,
            iterations: int = 50,
            standard_ddg: bool = True,
            high_resolution: bool = True,
            silent=True) -> Tuple[bool, pd.DataFrame]:
        """
        Perform a rosetta ddg_monomer calculation per the guidelines at
        https://www.rosettacommons.org/docs/latest/application_documentation/analysis/ddg-monomer

        :param iterations: 50 recommended by Rosetta.  Might reduce for quicker debugging
        :param standard_ddg:  # Seems always true.
        :param high_resolution: Boolean picks one of two ddg_monomer algoritms True->Row16; False->Row3
        :param silent:
        :return: (True if success,  dataframe of final results)
        """

        self._ddg_repo.make_calculation_directory_heirarchy()

        # Mirror the log - with INFO level details - to a local file
        
        # commandfile = "commands.txt"
        # Always use silent file/standard ddg for low quality models
        if not high_resolution:
            silent = True
            standard_ddg = True

        LOGGER.info('DDG_monomer.run called with\nsilent: %s\n,standard_ddg: %s\n,Row16: %s' % (
            silent, standard_ddg, high_resolution))

        # timestamp_suffix = datetime.datetime.now().strftime("%d%m%y%H%M%S")
        # isoresult, isopdbname, to_pdb = isolate_chain(self._pdb, self._chain)
        # if not isoresult:
        #    return self.failure(None, to_pdb)
        # else:
        #    self._to_pdb = to_pdb

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
         

        save_current_directory = os.getcwd()
        os.chdir(self._ddg_repo.calculation_dir)

        def move_prior_file_if_exists(filename):
            if os.path.exists(filename):
                archive_dir = "./archive"
                os.makedirs(archive_dir, exist_ok=True)
                # append the modification time of the file to the name
                # and archive in ./archive subdir
                archive_filename = os.path.join(
                    archive_dir,
                    # precede the filename with time stamp of the file's last modification.
                    datetime.datetime.fromtimestamp(os.path.getmtime(filename)).strftime("%Y%M%d%H%M%S") + "_%s"%filename
                )
                try:
                    os.replace(filename, archive_filename)
                    LOGGER.info("Saved prior run %s as %s" % (filename, archive_filename))
                except OSError:
                    LOGGER.exception("Failed to save prior run %s as %s: " % (filename, archive_filename))


        def command_ran_previously(binary_program_basename:str) -> Tuple[int,str,str]:
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
            stdout_filename,stderr_filename,exitcd_filename = self._command_result_filenames(binary_program_basename)
            previous_exit,previous_exit_int = self._get_previous_exit(exitcd_filename)

            previous_stdout = ""                        
            if os.path.isfile(stdout_filename):
                with open(stdout_filename, 'r') as stdout_record:
                    previous_stdout = stdout_record.read() # This stdout from previous successful run

            previous_stderr = ""                        
            if os.path.isfile(stderr_filename):
                with open(stderr_filename, 'r') as stderr_record:
                    previous_stderr = stderr_record.read() # Read stderr from previous successful run

            if previous_exit_int == 0:
                LOGGER.info("Previous successful run found\n.  NOT re-running:\n%s" % (
                            binary_program_basename))
            else:
                LOGGER.warning("Previous run did not exit with code 0\n.  Re-running:\n%s" % (
                            binary_program_basename))

            return (previous_exit_int,
                    previous_stdout,
                    previous_stderr) # runtime 0 may or may not be quite right


             


        def run_command_line_terminate_on_nonzero_exit(command_line_list: List[str],
                additional_files_to_archive: List[str]=[],
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
            stdout_filename,stderr_filename,exitcd_filename = self._command_result_filenames(binary_program_basename)

            # We've been playing with letting the command list have embedded spaces which have to be 
            # split back apart...
            # We may get bitten at some point because some of our arguments are actually 2 arguments
            # separated by a space.  For now...
            submitted_command_line_list = []
            for arg_with_spaces in command_line_list:
                submitted_command_line_list.extend(arg_with_spaces.split(' '))

            # Archive any additional output files once we commit to re-launching
            for additional_file in additional_files_to_archive: 
                move_prior_file_if_exists(additional_file) 



            # If there is already a file out there with exit code 0, then return the old status

            # For sanity, capture precicsely a reproducible shell command at work for future analysis
            commandline_record_filename = binary_program_basename + ".commandline_record.sh"
            move_prior_file_if_exists(commandline_record_filename)
            with open(commandline_record_filename, 'w') as commandline_record:
                commandline_record.write("#!/bin/bash\n")
                commandline_record.write("# Record of command invocation\n")
                commandline_record.write("# %s\n" % time.ctime(time.time()))
                commandline_record.write(" \\\n".join(command_line_list))
                # End the command with redirects to the stderr/stdout files
                commandline_record.write(" \\\n> %s 2> %s\n" % (stdout_filename, stderr_filename))
                commandline_record.write("echo Command exited with $?\n")

            command_start_time = time.perf_counter()
            LOGGER.info("Running via commands:\n%s" % ' '.join(submitted_command_line_list))
            # run using new API call, 
            completed_process = subprocess.run(submitted_command_line_list, shell=False, text=True, capture_output=True)
            fail_message = None
            if completed_process.returncode == 0:
                LOGGER.info("%s completed successfully (exit 0)", binary_program_basename)
            else:
                fail_message = "%s failed with exit %d" % (binary_program_basename, completed_process.returncode)
                LOGGER.critical(fail_message)

            # Record all stdout and stderr and exit code
            move_prior_file_if_exists(stdout_filename)
            with open(stdout_filename, 'w') as stdout_record:
                stdout_record.write(completed_process.stdout) # This writes the 'str'

            move_prior_file_if_exists(stderr_filename)
            with open(stderr_filename, 'w') as stderr_record:
                stderr_record.write(completed_process.stderr)

            move_prior_file_if_exists(exitcd_filename)
            with open(exitcd_filename, 'w') as exitcd_record:
                exitcd_record.write(str(completed_process.returncode))
            LOGGER.info("stdout/stderr/exit saved in files %s/%s/%s" % (stdout_filename, stderr_filename,exitcd_filename))

            if fail_message:
                self._ddg_repo.psb_status_manager.sys_exit_failure(fail_message)

            # Convert output byte streams to manageable unicode strings.
            # Return (stdout,stderr,elapsed_time)
            return (completed_process.returncode,
                    completed_process.stdout,
                    completed_process.stderr,
                    time.perf_counter() - command_start_time)

        #############################################################################################
        # Part 1 of 4        Preminimize the cleaned and renumbered PDB
        #############################################################################################

        # We will generate file of harmonic restraints - an input to the part 3 calculation
        minimize_directory = os.path.join(self._ddg_repo.structure_dir,"minimize")

        ddg_monomer_atom_pair_constraints_basename = "%s_%s.cst" % (
            self._ddg_repo.structure_id, self._ddg_repo.chain_id)
        ddg_monomer_atom_pair_constraints_fullpath = os.path.join(
            minimize_directory,ddg_monomer_atom_pair_constraints_basename)

        def run_part1_minimizer() -> Tuple[str, List[str]]:
            """
            :return minimized_pdb filename and orignal atom pair distnaces for constraints or halt on errors:
            """
            # 2020-June-28 Chris Moth Command copied from
            # https://www.rosettacommons.org/docs/latest/application_documentation/analysis/ddg-monomer
            self._verify_applications_available()

            save_curwd = os.getcwd()
            previous_exit_code = 1
            stderr = ""
            stdout = ""
            # If the minimize directory is not there
            # Or somehow there without good exit code (should be impossible)
            # Then rerun
            def load_from_prior_minimize():
                LOGGER.info("os.chdir('%s')"%minimize_directory)
                os.chdir(minimize_directory)
                previous_exit_code,stdout,stdin = command_ran_previously(self._minimize_application_filename)
                assert previous_exit_code  == 0,\
                    "%s directory lacks a 0 exit code recorded.  This should never happen"%\
                    minimize_directory
                assert os.path.isfile(ddg_monomer_atom_pair_constraints_basename),\
                    "%s directory lacks atom_pair constraints file %s"%(
                        minimize_directory,ddg_monomer_atom_pair_constraints_basename)
                LOGGER.info("Returning via os.chdir('%s')"%save_curwd)
                os.chdir(save_curwd)
                return previous_exit_code,stdout,stdin

            if os.path.isdir(minimize_directory):
                previous_exit_code,stdout,stderr = load_from_prior_minimize()

            minimized_pdb_filename = "minimized.%s_0001.pdb" % (
                os.path.basename(self._ddg_repo.cleaned_structure_filename).split(".")[0])

            if previous_exit_code == 0: 
                LOGGER.info("Skipping minimization.  Using results from prior calculation")
                returncode = previous_exit_code
                runtime = 0.0
            else:
                # Then re-run the minimization in a subdirectory of the variant current directory
                # If all goes well, this directory will be moved to remove the minimize_directory
                tmp_directory = 'tmp_minimize'
                os.makedirs(tmp_directory,mode=0o770,exist_ok=True)
                LOGGER.info("os.chdir('%s')"%tmp_directory)
                os.chdir(tmp_directory)

                # Premin only takes a list of structs even if you only have 1
                structure_list_filename = 'filename.input'
                with open(structure_list_filename,mode="w") as structure_list_fp:
                    structure_list_filename = structure_list_fp.name
                    structure_list_fp.write(self._ddg_repo.cleaned_structure_filename)

                # WAS  minimize_with_cst.default.linuxgccrelease
                minimize_with_cst_command = [
                    # Commented out parameters are in current ddG_monomer documentation
                    # but incompatible with current rosetta version.  Chris Moth
                    # discussed with Rocco - and we are staying with our old parameters
                    # even though they are likely 'wrong'er than they need to be
                    # 
                    # OLD command = minimize_with_cst.default.linuxgccrelease -in:file:fullatom -fa_max_dis 9.0 -ddg:harmonic_ca_tether 0.5 -ddg::constraint_weight 1.0 -ddg::out_pdb_prefix minimized -ddg::sc_min_only false -score:weights talaris2014 -in:file:l %s -overwrite -ignore_zero_occupancy false % (templistfile)
                    #



                    os.path.join(self._ddg_repo.rosetta_bin_dir, self._minimize_application_filename),
                    '-in:file:l {0}'.format(structure_list_filename),
                    '-in:file:fullatom',
                    '-ignore_unrecognized_res',
                    '-ignore_zero_occupancy false',
                    '-fa_max_dis 9.0',
                    # '-database {0}'.format(self._ddg_repo.rosetta_database_dir),  # the full oath to the database is required
                    '-ddg:harmonic_ca_tether 0.5',
                    '-ddg:out_pdb_prefix minimized',
                    # '-score:weights standard',
                    '-score:weights talaris2014',
                    '-ddg:constraint_weight 1.0',
                    # '-ddg:out_pdb_prefix minimized',  # Makes no sense to have file name have recommended min_cst_0.5
                    '-ddg:sc_min_only false'
                    # '-score:patch {0}'.format(os.path.join(self._ddg_repo.rosetta_database_dir,
                    #                                        "scoring", "weights", "score12.wts_patch"))
                    # '> mincst.log'
                    ]



                returncode, stdout, stderr, runtime = run_command_line_terminate_on_nonzero_exit(
                    minimize_with_cst_command,
                    additional_files_to_archive=[minimized_pdb_filename])

                # Collect data from successful run and create file for part 3

                # Grab the pairwise C-alpha distances (pre-minimized) as constraints to input to next
                original_atom_pair_distances_constraints = []
                for line in stdout.split("\n"):
                    if len(line) > 7 and line[:7] == "c-alpha":
                        cur = line.split()
                        original_atom_pair_distances_constraints.append(
                            "AtomPair CA %s CA %s HARMONIC %s %s" % (cur[5], cur[7], cur[9], cur[12]))

                if not original_atom_pair_distances_constraints:
                    self._ddg_repo.psb_status_manager.sys_exit_failure(
                        "No AtomPair contraints returned from %s" % minimize_application_filename)

                with open(ddg_monomer_atom_pair_constraints_basename, 'w') as outfile:
                    outfile.write("\n".join(original_atom_pair_distances_constraints))
                    outfile.write("\n")

                LOGGER.info("%d atom pair constraints for ddg_monomer written to %s"%(
                    len(original_atom_pair_distances_constraints),
                    ddg_monomer_atom_pair_constraints_fullpath))

                LOGGER.info("Returning via os.chdir(%s)"%save_curwd)
                os.chdir(save_curwd)

                # Move our work to become the new "minimize" directory
                # so that other tasks do not have to repeat this work
                #
                # If OTHER tasks beat us to installing their work directory, then 
                # no worries!  Load that as the data source instead
                # and discard our hard work
                LOGGER.info("Attempting os.rename(%s,%s)"%(tmp_directory,minimize_directory))
                rename_succeeded = False
                if returncode == 0:
                    try:
                        os.rename(tmp_directory,minimize_directory)
                        rename_succeeded = True
                    except OSError:
                        # The final minimize directory already exists there
                        pass

                if not rename_succeeded:
                    # Then some other task beat us...
                    LOGGER.info("%s already installed by another process."%minimize_directory)
                    LOGGER.info("Discarding current calculation and loading prior results")
                    # Load their outputs so all the ddgs are on the same page
                    returncode, stdout,stderr = load_from_prior_minimize()
                    assert returncode == 0 


            _minimized_pdb_relpath = os.path.join('..','..',minimize_directory,minimized_pdb_filename)
            if not os.path.isfile(_minimized_pdb_relpath):
                self._ddg_repo.psb_status_manager.sys_exit_failure('Minimized pdb file %s was not created by %s' % (
                    os.path.abspath(_minimized_pdb_relpath), minimize_application_filename ))

            return _minimized_pdb_relpath

        self._ddg_repo.psb_status_manager.write_info("Part 1: Minimization")
        minimized_pdb_relpath = run_part1_minimizer()

        LOGGER.info("Minimized structure %s available for part 2", 
            minimized_pdb_relpath)

        #############################################################################################
        # Part 2 of 4        Rescore minimized model and get residue per-residue score  (quasi energy)
        #############################################################################################
        def run_part2_rescore():
            self._verify_applications_available()
            min_residues_filename = 'min_residues.sc'  # Not sure why we'd want timestamp_suffix...
            per_residues_directory = os.path.join(self._ddg_repo.structure_dir,"per_residue_energies")
            save_curwd = os.getcwd()
            previous_exit_code = 1
            stderr = ""
            stdout = ""


            #######################################################################
            def load_raw_per_residue_scores():
                raw_per_residue_scores = []
                with open(min_residues_filename) as min_residues_f:
                    for line in min_residues_f.readlines():
                        line_components = line.strip().split()
                        if len(line_components) > 2:
                            if line_components[-1] != 'description':
                                raw_per_residue_scores.append(line_components[-2])

                # Make sure that we have a score for each residue
                if len(raw_per_residue_scores) != len(self._ddg_repo.residue_to_clean_xref):
                    self._ddg_repo.psb_status_manager.sys_exit_failure("%s yielded %d raw_scores.  %d were expected" % (
                        self._per_residues_application_filename,
                        len(raw_per_residue_scores),
                        len(self._ddg_repo.residue_to_clean_xref)))
                return raw_per_residue_scores

            # If the per_residues directory is not there
            # Or somehow there without good exit code (should be impossible)
            # Then rerun
            def load_from_prior_per_residues_run():
                LOGGER.info("os.chdir('%s')"%per_residues_directory)
                os.chdir(per_residues_directory)
                previous_exit_code,stdout,stderr = command_ran_previously(self._per_residues_application_filename)
                assert previous_exit_code == 0,\
                    "%s directory lacks a 0 exit code recorded.  This should never happen"%\
                    per_residues_directory

                raw_scores = load_raw_per_residue_scores()
                LOGGER.info("Returning via os.chdir('%s')"%save_curwd)
                os.chdir(save_curwd)
                return previous_exit_code,stdout,stderr,raw_scores
            #######################################################################

            if os.path.isdir(per_residues_directory):
                previous_exit_code,stdout,stderr,raw_scores = load_from_prior_per_residues_run()
                returncode = previous_exit_code
                runtime = 0.0
            else:
                # Then re-run the minimization in a subdirectory of the variant current directory
                # If all goes well, this directory will be moved to remove the minimize_directory
                tmp_directory = tempfile.mkdtemp(prefix='tmp_per_residues',dir='.')
                os.makedirs(tmp_directory,mode=0o770,exist_ok=True)
                LOGGER.info("os.chdir('%s')"%tmp_directory)
                os.chdir(tmp_directory)

                # Then re-run the minimization in a subdirectory of the variant current directory
                rescore_command = [
                    os.path.join(self._ddg_repo.rosetta_bin_dir, self._per_residues_application_filename),
                    '-s ' + minimized_pdb_relpath,
                    '-score:weights talaris2014',
                    '-out:file:silent %s' % min_residues_filename
                ]

                # stdout, stderr, runtime
                returncode, _, _, _ = run_command_line_terminate_on_nonzero_exit(
                    rescore_command,
                    additional_files_to_archive=[min_residues_filename])

                raw_scores = load_raw_per_residue_scores()

                LOGGER.info("Returning via os.chdir(%s)"%save_curwd)
                os.chdir(save_curwd)

                # Move our work to become the new "per_residues" directory
                # so that other tasks do not have to repeat this work
                #
                # If OTHER tasks beat us to installing their work directory, then 
                # no worries!  Load that as the data source instead
                # and discard our hard work
                LOGGER.info("Attempting os.rename('%s','%s')"%(tmp_directory,per_residues_directory))
                rename_succeeded = False
                if returncode == 0:
                    try:
                        os.rename(tmp_directory,per_residues_directory)
                        rename_succeeded = True
                    except OSError:
                        # The final per_residues directory already exists there
                        pass

                if not rename_succeeded:
                    # Then some other task beat us...
                    LOGGER.info("%s already installed by another process."%per_residues_directory)
                    LOGGER.info("Discarding current calculation and loading prior results")
                    # Load their outputs so all the ddgs are on the same page
                    returncode, stdout,stderr,raw_scores = load_from_prior_per_residues_run()
                    assert returncode == 0 


            mutation_scores = {}
            for original_residue_id,rosetta_residue_id in self._ddg_repo.residue_to_clean_xref.items():
                # Extract the residue # from biopython resid tuple
                rosetta_residue_number = int(rosetta_residue_id[1])
                # Rosetta residues start with 1.  Scores indexed from 0
                raw_score = raw_scores[rosetta_residue_number - 1]
                if not raw_score:
                    self._ddg_repo.psb_status_manager.sys_exit_failure("Mutation residue %s result is None - investigate" % mutation)

                LOGGER.debug("Mutation score for %d %s=%s"%(
                    rosetta_residue_number,
                    original_residue_id,
                    raw_score))
                mutation_scores[original_residue_id] = float(raw_score)

            return mutation_scores

        self._ddg_repo.psb_status_manager.write_info("Part 2: Rescore all residues")
        mutation_scores = run_part2_rescore()

        #############################################################################################
        # Part 3 of 4        Run ddg_mononomer
        # Unlike parts 1 and 2, this run is unique to each mutation
        # We assume that use of the parent directory here will be unique
        #############################################################################################


        def run_part3_ddg():
            self._verify_applications_available()
            # As other steps, don't re-run ddg monomer if we already have run it successfully
            previous_exit_code,stdout,stdin = command_ran_previously(self._ddg_monomer_application_filename)
            ddg_predictions_filename = 'ddg_predictions.out'

            # Does it look like we exited with no error before?
          
            previous_run_acceptable = False 
            if previous_exit_code == 0:
                previous_run_acceptable = True

                if not os.path.isfile(ddg_predictions_filename):
                    previous_run_acceptable = False
                    LOGGER.warning("Prior ddg_monomer run missing %s.  Restarting"%\
                        ddg_predictions_filename)

                if not os.path.isfile(self.final_results_filename()):
                    previous_run_acceptable = False
                    LOGGER.warning("Prior ddg_monomer run missing %s.  Restarting"%\
                        self.final_results_filename())

            if previous_run_acceptable:
                # report back from prior run
                returncode = previous_exit_code
                final_results_df = pd.read_csv(self.final_results_filename(), sep='\t', dtype={'ddG': float, 'WT_Res_Score': float})
            else:
                # Start a new run of ddg_monomer
                # Get rid of the predictions output, else ddg appends to it - problems!
                try:
                    os.unlink(ddg_predictions_filename)
                except OSError:
                    pass

                # Create a list of mutations to be analyzed by ddg_monomer, using Rosetta 1..N numberings.
                # Documented https://www.rosettacommons.org/docs/latest/application_documentation/analysis/ddg-monomer
                mutation_list_filename = "%s_%s_%s.mut" % (
                    self._ddg_repo.structure_id, self._ddg_repo.chain_id, self._mutationtext)
                with open(mutation_list_filename, 'w') as mutation_list_f:
                    mutation_list_f.write("total %d\n" % len(rosetta_mutations))
                    for rosetta_mutation in rosetta_mutations:
                        mutation_list_f.write("1\n")
                        mutation_list_f.write("%s %d %s\n" % (
                                              rosetta_mutation[0],
                                              int(rosetta_mutation[1:-1]), 
                                              rosetta_mutation[-1]))

                #
                # https://www.rosettacommons.org/docs/latest/application_documentation/analysis/ddg-monomer
                # STart with obtions common to both protocols
                ddg_options_filename = "%s_%s_%s.options" % (
                    self._ddg_repo.structure_id,
                    self._ddg_repo.chain_id,
                    self._mutationtext)

                ddg_options = [
                    '-in:file:s ' + minimized_pdb_relpath,  # PDB file of structure on which point mutations should be made
                    '-in::file::fullatom',  # read the input PDB file as a fullatom structure
                    '-ddg::mut_file ' + mutation_list_filename,  # the list of point mutations to consider in this run
                    '-fa_max_dis 9.0',  # optional -- if not given, the default value of 9.0 Angstroms is used.
                    '-ddg::dump_pdbs true',  # write one PDB for the wildtype and one for the pointmutant for each iteration
                    '-database ' + self._ddg_repo.rosetta_database_dir,  # the full oath to the database is required
                    '-ddg::iterations %d' % iterations,  # Typically 50 iterations of the algorithm
                    '-mute all',  # silence all of the log-file / stdout output generated by this protocol
                    '-ignore_unrecognized_res',  # ignore Rosetta-foreign res instead of quitting with error
                    '-ignore_zero_occupancy false', # Plenty of good CryoEM structures publish with 0 occupancy
                    '-ddg::suppress_checkpointing true',  # don't checkpoint LIZ DOES CHECKPOINTING WORK AT ALL?
                    '-ddg::output_silent true',  # write output to a silent file
                    '-ddg:weight_file soft_rep_design'  # soft-repulsive weights for initial sidechain optimization stage'
                ]

                if high_resolution:
                    protocol = "Row16"
                    comment1 = "#High resolution protocol"
                    comment2 = "#Based on kellogg et al 2011 table 1 row 16"
                    # Add Row 16-specific options
                    ddg_options.extend([
                        # optional -- the weights file to use, if not given, then
                        # "score12" will be used (score12 = standard.wts + score12.wts_patch)
                        # '-ddg:minimization_scorefunction <weights file>',

                        # optional -- the weight-patch file to apply to the weight file; does not have to be given"
                        # -ddg::minimization_patch <weights patch file>

                        # recommended: local optimization restricts the sidechain optimization to only the
                        # 8 A neighborhood of the mutation (equivalent to row 13)
                        '-ddg::local_opt_only false',

                        '-ddg::min_cst true',  # use distance restraints (constraints) during backbone minimization phase
                        # the set of constraints gleaned from the original (non-pre-relaxed/minimized) structure
                        '-constraints::cst_file ' + ddg_monomer_atom_pair_constraints_fullpath,
                        '-ddg::mean false',  # do not report the mean energy
                        '-ddg::min true',  # report the minimum energy
                        '-ddg::sc_min_only false',  # don't minimize only the backbone during backbone minimization phase

                        # perform three rounds of minimization (and not just the default 1 round)
                        # where the weight on the repulsive term is increased from 10% to 33% to 100%
                        '-ddg::ramp_repulsive true',
                        '-unmute core.optimization.LineMinimizer'  # optional -- unsilence a particular tracer

                        # Options we USED to use:
                        # -ddg:minimization_scorefunction talaris2014
                        # - ignore_zero_occupancy false
                    ])
                    LOGGER.info("Running high quality (row 16) ddG monomer protocol: %d iterations" % iterations)
                else:
                    protocol = "Row3"
                    comment1 = "#Low resolution protocol"
                    comment2 = "#Based on kellogg et al 2011 table 1 row 3"
                    # Add Row 3-specific options
                    ddg_options.extend([
                        # repack the residues in an 8 Angstrom shell around the site of the point mutation
                        '-ddg::local_opt_only true',
                        '-ddg::mean true # do not report the mean energy',
                        '-ddg::min false # report the minimum energy'
                        # Options we USED to use
                        # -ignore_zero_occupancy false
                        # -ddg:minimization_scorefunction talaris2014
                        # Worried about -constraints file
                    ])

                    LOGGER.info("Running low quality (row 3) ddG monomer protocol: %d iterations" % iterations)

                move_prior_file_if_exists(ddg_options_filename)
                with open(ddg_options_filename, 'w') as ddg_options_f:
                    ddg_options_f.write("%s\n" % comment1)
                    ddg_options_f.write("%s\n" % comment2)
                    for ddg_option in ddg_options:
                        ddg_options_f.write("%s\n" % ddg_option)

                # Run actual ddg_monomer application

                # stdout, stderr, runtime =
                returncode, _, _, _ = run_command_line_terminate_on_nonzero_exit([
                    os.path.join(self._ddg_repo.rosetta_bin_dir, self._ddg_monomer_application_filename),
                    "@" + ddg_options_filename],
                    additional_files_to_archive=[ddg_predictions_filename])
                 
                ######################################################3
                # Make sure all expected files have been generated
                # From what I can tell, this branch is abandoned - but wire it in just in case it might be resurrected soon
                ######################################################3
                if not silent:
                    # Cool - looks like one could get all 50 (#iterations) pdbs
                    # But this code just looks way wrong.
                    expected_files = [
                            "repacked_wt_round_%d.pdb" % (x + 1) for x in range(iterations)] + \
                                     [ddg_predictions_filename]
                    for mut in rosetta_mutations:
                        expected_files += ["mut_%s_round_%d.pdb" % ("".join(str(x) for x in mut), x + 1) for x in
                                           range(iterations)]
                else:
                    # What we normally get
                    expected_files = [ddg_predictions_filename]

                for expected_file in expected_files:
                    if not os.path.exists(expected_file):
                        self._ddg_repo.psb_status_manager.sys_exit_failure(
                            "Something went wrong with ddg_monomer, expected file %s is missing!"%expected_file)

                # These files are always empty - double check
                if os.path.isfile("wt_traj"):
                   if os.stat("wt_traj").st_size == 0:
                       os.remove("wt_traj")

                # for mut in mutation_pose:
                #    if os.path.isfile("mutant_traj%s" % "".join(str(x) for x in mut)):
                #        os.remove("mutant_traj%s" % "".join(str(x) for x in mut))

                ddg_predictions_df = pd.read_csv(
                    ddg_predictions_filename,
                    delim_whitespace=True)

                # Ultimate target is final_results file with everything we "want"
                final_results_df = pd.DataFrame(
                    columns=["File", "Chain", "Residue", "WT_Res_Score", "Mutation", "ddG", "Protocol"])

                expected = len(self._mutations)

                if ddg_predictions_df.shape[0] != expected:
                    self._ddg_repo.psb_status_manager.sys_exit_failure(
                        "Something went wrong with ddg_monomer, expected %d pridictions, got only %d!" % (
                            expected, ddg_predictions_df.shape[0])
                    )

                mutation_ddg = {}
                # temp = mutation_pose[:]
                # temp.sort(key=lambda x: x[1])
                # tempraw = "".join("".join(str(x) for x in y) for y in sorted(self._mutations, key=lambda z: z[1]))
                # Standard ddg uses the output ddG score as provided by application
                # Nonstandard ddg uses difference in pose scores between the mean score of the top 3 models per condition.
                if standard_ddg:
                    for index,row in ddg_predictions_df.iterrows():
                        rosetta_mutation = row['description']

                        # Populate the ddgs back to the original structure positions (perhaps with insertion codes)

                        original_residue = rosetta_residue_number_to_residue_xref[int(rosetta_mutation[1:-1])]
                        original_residue_insertion_code = original_residue[2] if original_residue[2].isalpha() else ''
                        original_mutation = (
                                rosetta_mutation[0] +
                                str(original_residue[1]) +
                                original_residue_insertion_code +
                                rosetta_mutation[-1])

                        mutation_ddg[original_mutation] = float(row['total'])

                        final_results_df = final_results_df.append({
                            "File": self._ddg_repo.structure_id,
                            "Chain": self._ddg_repo.chain_id,
                            "Residue": str(original_residue[1]) + original_residue_insertion_code,
                            "WT_Res_Score": mutation_scores[original_residue],
                            "Mutation": original_mutation,
                            "ddG": mutation_ddg[original_mutation],
                            "Protocol": protocol
                        },ignore_index=True)

                else:
                    self._ddg_repo.psb_status_manager.sys_exit_failure(
                        "NOn Standard DDG is not supported in 2020 see ddg_monomer.py")
                    # scout, scres = self.get_ddg(high_resolution, silent)
                    # if not scout:
                    #     return [False, scres]
                    # for line in scres:
                    #    mut, dg = line.split("\t")
                    #    try:
                    #        curmut = pose_to_raw[mut]
                    #    except KeyError:
                    #        curmut = tempraw
                    #    mutation_ddg[curmut] = dg

                # full_final_results_filename = os.path.join(save_current_directory, self.final_results_filename())

                final_results_df.to_csv(self.final_results_filename(), sep='\t')
                LOGGER.info("Results saved to %s" % os.path.abspath(self.final_results_filename()))
                # End of ddg_monomer calculation
            return True, final_results_df

        self._ddg_repo.psb_status_manager.write_info("Part 3: DDG monomer")
        success,final_results_df = run_part3_ddg()
        if not success: 
            self._ddg_repo.psb_status_manager.sys_exit_failure("run_part3_ddg() failed - not sure why it got to this point")

        self._ddg_repo.psb_status_manager.mark_complete()

        # self._log_close()
        os.chdir(save_current_directory)
        return True, final_results_df

    def jobstatus(self) -> Dict[str,str]:
        jobstatus_info = {}
        save_current_directory = os.getcwd()
        try:
            os.chdir(self._ddg_repo.calculation_dir)
        except FileNotFoundError:
            jobstatus_info['ExitCode'] = None
            jobstatus_info['jobinfo'] = 'ddg_monomer has not created %s'%self._ddg_repo.calculation_dir
            jobstatus_info['jobprogress'] = jobstatus_info['jobinfo']
            return jobstatus_info

        _,_,exitcd_filename = self._command_result_filenames(self._ddg_monomer_application_filename)
        previous_exit,previous_exit_int = self._get_previous_exit(exitcd_filename)
        if previous_exit_int == 0:
            jobstatus_info['ExitCode'] = '0'
            jobstatus_info['jobprogress'] = 'Completed' # For now, we will overwrite memory on happy exit - need to think
            jobstatus_info['jobinfo'] = 'Completed'
        else:
            jobstatus_info['ExitCode'] = previous_exit
            info,progress = self._ddg_repo.psb_status_manager.read_info_progress()
            jobstatus_info['jobinfo'] = info
            jobstatus_info['jobprogress'] = progress

        os.chdir(save_current_directory)
        return jobstatus_info

    def retrieve_result(self) -> pd.Series:
        """Return the calculated ddG results as a pandas series, or None"""
        _,_,exitcd_filename = self._command_result_filenames(self._ddg_monomer_application_filename)
        save_current_directory = os.getcwd()
        try:
            os.chdir(self._ddg_repo.calculation_dir)
        except FileNotFoundError:
            LOGGER.info("No results directory %s",self._ddg_repo.calculation_dir)
            return None

        previous_exit,previous_exit_int = self._get_previous_exit(exitcd_filename)

        final_results_fullpath = os.path.join(self._ddg_repo.calculation_dir,self.final_results_filename())

        if previous_exit_int == 0:
            final_results_df = pd.read_csv(final_results_fullpath, sep='\t', dtype={'ddG': float, 'WT_Res_Score': float})
            assert len(final_results_df) == 1,"Format error: single row of data not found in %s - halting"%final_results_fullpath
            final_results_df['RefAA'] = final_results_df['Mutation'].str[0]
            final_results_df['AltAA'] = final_results_df['Mutation'].str[-1]
            LOGGER.debug("ddg of %f read from %s"%(final_results_df['ddG'],final_results_fullpath))
        else:
            LOGGER.info("ddg results file %s not found"%(final_results_fullpath))
            final_results_df = None

        os.chdir(save_current_directory)
        return final_results_df.iloc[0] if final_results_df is not None else None
        
