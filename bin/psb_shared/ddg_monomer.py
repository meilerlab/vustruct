#!/usr/bin/env python
"""
class DDG_monomer manages all aspects of Rosetta ddg_monomer calculations 
inside the directory heirarchy supplied by a ddg_repo object
"""
import os
import datetime
# import gzip
import pprint
import subprocess
import tempfile
import time
import pandas as pd
from bin.psb_shared.ddg_repo import DDG_repo
from typing import Dict, List, Tuple, Union
# from sys import argv, stderr, stdout
# from os import popen, system
# from os.path import exists, basename
from bin.psb_shared.psb_progress import PsbStatusManager

import logging
from logging.handlers import RotatingFileHandler

LOGGER = logging.getLogger(__name__)


class DDG_monomer(object):
    def __init__(self, ddg_config_dict: Dict[str, str], ddg_repo: DDG_repo, mutations: Union[str, List[str]]):
        """
        :param ddg_config_dict:
        :param ddg_repo:
        :param mutations: S123AF format (or list of these) that indicate AA, res#, insert code, new AA.
        """

        self._rosetta_bin_dir = ddg_config_dict['rosetta_bin_dir']
        self._rosetta_database_dir = ddg_config_dict['rosetta_database_dir']
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

        #        self._root = root_path
        for application in ['minimize_with_cst.linuxgccrelease', 'per_residue_energies.linuxgccrelease',
                            'ddg_monomer.linuxgccrelease', 'score_jd2.linuxgccrelease']:
            application_fullpath = os.path.join(self._rosetta_bin_dir, application)
            assert os.path.exists(application_fullpath), (
                    "%s is a required application for ddg_monomer, but not found" % application_fullpath)

        self._application_timers = ['minimize', 'rescore', 'ddg_monomer']
        self._to_pdb = []

        need_roll = os.path.isfile(ddg_repo.log_filename)

        self._local_fh = RotatingFileHandler(ddg_repo.log_filename, backupCount=5)

        # self._local_fh.setFormatter(formatter)
        self._local_fh.setLevel(logging.INFO)
        LOGGER.addHandler(self._local_fh)

        if need_roll:
            self._local_fh.doRollover()

        LOGGER.info("ddg_mmonomer.__init__ completed after setting:\n%s" % pprint.pformat(self.__dict__))

    def final_results_filename(self) -> str:
        return "%s_%s_%s" % (
            self._ddg_repo.structure_id,
            self._ddg_repo.chain_id,
            self._mutationtext
        )

    def _log_close(self):
        if self._local_fh is not None:
            self._local_fh.flush()
            self._local_fh.close()
            LOGGER.removeHandler(self._local_fh)
        self._local_fh = None

    #        with open(self._pdb) as infile:
    #            self._pdblines = infile.readlines()

    def __del__(self):
        if self._local_fh:
            # LOGGER.warning("ddg_monomer.__del__ called with logfile still open")
            self._log_close()

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
            psb_status_manager: PsbStatusManager,
            iterations: int = 50,
            standard_ddg: bool = True,
            high_resolution: bool = True,
            silent=True) -> Tuple[bool, pd.DataFrame]:
        """
        Perform a rosetta ddg_monomer calculation per the guidelines at
        https://www.rosettacommons.org/docs/latest/application_documentation/analysis/ddg-monomer

        :param psb_status_manager: Utility class for recording status to the compute cluster, and failure messages
        :param iterations: 50 recommended by Rosetta.  Might reduce for quicker debugging
        :param standard_ddg:  # Seems always true.
        :param high_resolution: Boolean picks one of two ddg_monomer algoritms True->Row16; False->Row3
        :param silent:
        :return: (True if success,  dataframe of final results)
        """
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
            psb_status_manager.sys_exit_failure("No mutations found in parameter self._mutations %s" % self._mutations)

        rosetta_residue_number_to_residue_xref = {
            int(clean[1]): residue for residue, clean in self._ddg_repo.residue_to_clean_xref
        }

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
                    # add the time stamp f the file's last modification.
                    filename + datetime.datetime.fromtimestamp(os.path.getmtime(filename)).strftime("%Y%M%d%H%M%S")
                )
                try:
                    os.replace(filename, archive_filename)
                    LOGGER.info("Saved %s as %s" % (filename, archive_filename))
                except OSError:
                    LOGGER.exception("Failed to save %s as %s: " % (filename, archive_filename))

        def run_command_line(command_line_list: List[str]) -> Tuple[str, str, float]:
            """
            Launch a shell to run the command line.
            Terminate if non-zero return value

            :param command_line_list: list of command line components
            :return (
                str: stdout from command
                str: stderr from command
                float: run time
            """

            # Example: minimzer.linuxgccrelease
            binary_program_basename = os.path.basename(command_line_list[0])
            stdout_filename = binary_program_basename + ".stdout"
            stderr_filename = binary_program_basename + ".stderr"

            # For sanity, capture precicsely a reproducible shell command at work for future analysis
            commandline_record_filename = binary_program_basename + ".commandline_record.sh"
            move_prior_file_if_exists(commandline_record_filename)
            with open(commandline_record_filename, 'w') as commandline_record:
                commandline_record.write("#!/bin/bash\n")
                commandline_record.write("# Record of command invocation")
                commandline_record.write("# %s\n" % time.ctime(time.time()))
                commandline_record.write(" \\\n".join(command_line_list))
                # End the command with redirects to the stderr/stdout files
                commandline_record.write(" \\\n> %s >> %s\n" % (stdout_filename, stderr_filename))

            command_start_time = time.perf_counter()
            # We may get bitten at some point because some of our arguments are actually 2 arguments
            # separated by a space.  For now...
            run_process = subprocess.Popen(command_line_list, shell=True, stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)
            LOGGER.info("Running via shell command:\n%s" % '\n'.join(run_process.args))
            (stdout_binary, stderr_binary) = run_process.communicate()
            if run_process.returncode == 0:
                LOGGER.info("%s completed successfully (exit 0)", binary_program_basename)
            else:
                message = "%s failed with exit %d" % (binary_program_basename, run_process.returncode)
                psb_status_manager.sys_exit_failure(message)

            stdout_str = stdout_binary.decode('latin')
            stderr_str = stderr_binary.decode('latin')

            # Record all stdout and stderr
            move_prior_file_if_exists(stdout_filename)
            with open(stdout_filename, 'w') as stdout_record:
                stdout_record.write(stdout_str)

            move_prior_file_if_exists(stderr_filename)
            with open(stderr_filename, 'w') as stderr_record:
                stderr_record.write(stderr_str)

            LOGGER.info("stdout/err saved in %s/%s" % (stdout_filename, stderr_filename))

            # Convert output byte streams to manageable unicode strings.
            # Return (stdout,stderr,elapsed_time)
            return (stdout_str,
                    stderr_str,
                    time.perf_counter() - command_start_time)

        #############################################################################################
        # Part 1 of 4        Preminimize the cleaned and renumbered PDB
        #############################################################################################

        def run_part1_minimizer() -> Tuple[str, List[str]]:
            """
            :return minimized_pdb filename and orignal atom pair distnaces for constraints or halt on errors:
            """
            # Premin only takes a list of structs even if you only have 1
            temp_structure_list_filename = None
            with tempfile.NamedTemporaryFile(delete=False, dir=".", mode="w") as temp_structure_list_fp:
                temp_structure_list_filename = temp_structure_list_fp.name
                temp_structure_list_fp.write(self._ddg_repo.cleaned_structure_filename)

            # 2020-June-28 Chris Moth Command copied from
            # https://www.rosettacommons.org/docs/latest/application_documentation/analysis/ddg-monomer
            minimize_application_filename = "minimize_with_cst.linuxgccrelease"
            # WAS  minimize_with_cst.default.linuxgccrelease
            minimize_with_cst_command = [
                os.path.join(self._rosetta_bin_dir, minimize_application_filename),
                '-in:file:l {0}'.format(temp_structure_list_filename),
                '-in:file:fullatom',
                '-ignore_unrecognized_res',
                '-fa_max_dis 9.0',
                '-database {0}'.format(self._rosetta_database_dir),  # the full oath to the database is required
                '-ddg::harmonic_ca_tether 0.5',
                '-score:weights standard',
                '-ddg::constraint_weight 1.0',
                '-ddg::out_pdb_prefix minimized',  # Makes no sense to have file name have recommended min_cst_0.5
                '-ddg::sc_min_only false',
                '-score:patch {0}'.format(os.path.join(self._rosetta_database_dir,
                                                       "scoring", "weights", "score12.wts_patch"))
                # '> mincst.log'
                # Prior different options included:
                # -score:weights talaris2014
                ]

            stdout, stderr, runtime = run_command_line(minimize_with_cst_command)

            # Grab the pairwise C-alpha distances (pre-minimized) as constraints to input to next step.
            _original_atom_pair_distances_constraints = []
            for line in stdout.split("\n"):
                if len(line) > 7 and line[:7] == "c-alpha":
                    cur = line.split()
                    _original_atom_pair_distances_constraints.append(
                        "AtomPair CA %s CA %s HARMONIC %s %s" % (cur[5], cur[7], cur[9], cur[12]))

            if not _original_atom_pair_distances_constraints:
                psb_status_manager.sys_exit_failure(
                    "No AtomPair contraints returned from %s" % minimize_application_filename)

            os.remove(temp_structure_list_filename)
            minimized_pdb_filename = "minimized.%s_0001.pdb" % (
                os.path.basename(self._ddg_repo.cleaned_structure_filename).split(".")[0])

            if not os.path.isfile(minimized_pdb_filename):
                psb_status_manager.sys_exit_failure('Minimized pdb file %s not created by %s' % (
                    minimized_pdb_filename, minimize_application_filename ))

            return minimized_pdb_filename, _original_atom_pair_distances_constraints

        minimized_pdb_filename, original_atom_pair_distances_constraints = run_part1_minimizer()

        LOGGER.info("%s available for part 2", minimized_pdb_filename)

        #############################################################################################
        # Part 2 of 4        Rescore model and get residue score
        #############################################################################################
        def run_part2_rescore():
            min_residues_filename = 'min_residues.sc'  # Not sure why we'd want timestamp_suffix...
            per_residues_application_filename = 'per_residue_energies.linuxgccrelease'
            rescore_command = [
                os.path.join(self._rosetta_bin_dir, per_residues_application_filename),
                '-s ' + minimized_pdb_filename,
                '-score:weights talaris2014',
                '-out:file:silent %s' % min_residues_filename
            ]

            # stdout, stderr, runtime
            _, _, _ = run_command_line(rescore_command)

            raw_scores = []
            with open(min_residues_filename) as min_residues_f:
                for line in min_residues_f.readlines():
                    line_components = line.strip().split()
                    if len(line_components) > 2:
                        if line_components[-1] != 'description':
                            raw_scores.append(line_components[-2])

            # Make sure that we have a score for each residue
            if len(raw_scores) != len(self._ddg_repo.residue_to_clean_xref):
                psb_status_manager.sys_exit_failure("%s yielded %d raw_scores.  %d were expected" % (
                    per_residues_application_filename,
                    len(raw_scores),
                    len(self._ddg_repo.residue_to_clean_xref)))

            mutation_scores = {}

            for mutation, mutation_resid in zip(self._mutations, mutation_resids):
                rosetta_residue = self._ddg_repo.residue_to_clean_xref[mutation_resid]
                # Extract the residue # from biopython tuple
                rosetta_residue_number = int(rosetta_residue[1])
                raw_score = raw_scores[rosetta_residue_number - 1][0]
                if not raw_score:
                    psb_status_manager.sys_exit_failure("Mutation residue %s result is None - investigate" % mutation)
                mutation_scores[mutation] = raw_score

                return mutation_scores

        mutation_scores = run_part2_rescore()

        #############################################################################################
        # Part 3 of 4        Run ddg_mononomer
        #############################################################################################

        # Generate file of harmonic restraints
        atom_pair_constraints_filename = "%s_%s_%s.cst" % (
            self._ddg_repo.structure_id, self._ddg_repo.chain_id, self._mutationtext)
        with open(atom_pair_constraints_filename, 'w') as outfile:
            outfile.write("\n".join(original_atom_pair_distances_constraints))
            outfile.write("\n")

        def run_part3_ddg():
            # Create a list of mutations to be analyzed by ddg_monomer, using Rosetta 1..N numberings.
            # Documented https://www.rosettacommons.org/docs/latest/application_documentation/analysis/ddg-monomer
            mutation_list_filename = "%s_%s_%s.mut" % (
                self._ddg_repo.structure_id, self._ddg_repo.chain_id, self._mutationtext)
            with open(mutation_list_filename, 'w') as mutation_list_f:
                mutation_list_f.write("total %d\n" % len(rosetta_mutations))
                for rosetta_mutation in rosetta_mutations:
                    mutation_list_f.write("1\n")
                    mutation_list_f.write("%d %s %s\n" % (
                                          int(rosetta_mutation[1:-1]), rosetta_mutation[0], rosetta_mutation[-1]))
                # if len(self._mutations) > 1:
                #    outfile.write("%d\n" % len(mutation_pose))
                #    for mut in mutation_pose:
                #        outfile.write(" ".join(str(x) for x in mut))
                #        outfile.write("\n")

            #
            # https://www.rosettacommons.org/docs/latest/application_documentation/analysis/ddg-monomer
            # STart with obtions common to both protocols
            ddg_options_filename = "%s_%s_%s.options" % (
                self._ddg_repo.structure_id,
                self._ddg_repo.chain_id,
                self._mutationtext)

            ddg_options = [
                '-in:file:s ' + minimized_pdb_filename,  # PDB file of structure on which point mutations should be made
                '-in::file::fullatom',  # read the input PDB file as a fullatom structure
                '-ddg::mut_file ' + mutation_list_filename,  # the list of point mutations to consider in this run
                '-fa_max_dis 9.0',  # optional -- if not given, the default value of 9.0 Angstroms is used.
                '-ddg::dump_pdbs true',  # write one PDB for the wildtype and one for the pointmutant for each iteration
                '-database ' + self._rosetta_database_dir,  # the full oath to the database is required
                '-ddg::iterations %d' % iterations,  # Typically 50 iterations of the algorithm
                '-mute all',  # silence all of the log-file / stdout output generated by this protocol
                '-ignore_unrecognized_res',  # ignore Rosetta-foreign res instead of quitting with error
                '-ddg::suppress_checkpointing true',  # don't checkpoint LIZ DOES CHECKPOINTING WORK AT ALL?
                '-ddg::output_silent true',  # write output to a silent file
                '-ddg:weight_file soft_rep_design'  # soft-repulsive weights for initial sidechain optimization stage'
            ]

            if high_resolution:
                protocol = "Row16"
                comment = "#Based on kellogg et al 2011 table 1 row 16"
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
                    '- constraints::cst_file ' + atom_pair_constraints_filename,
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
                comment = "#Based on kellogg et al 2011 table 1 row 3"
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
                ddg_options_f.write("%s\n" % comment)
                for ddg_option in ddg_options:
                    ddg_options_f.write("%s\n" % ddg_option)

            # Run actual ddg_monomer application
            move_prior_file_if_exists("ddg_predictions.out")

            # stdout, stderr, runtime =
            _, _, _ = run_command_line([
                os.path.join(self._rosetta_bin_dir, "ddg_monomer.linuxgccrelease"),
                "@" + ddg_options_filename])

            # Make sure all expected files have been generated
            # From what I can tell, this branch is abandoned - but wire it in just in case it might be resurrected soon
            if not silent:
                # Cool - looks like one could get all 50 (#iterations) pdbs
                # But this code just looks way wrong.
                expected_files = [
                        "repacked_wt_round_%d.pdb" % (x + 1) for x in range(iterations)] + \
                                 ["ddg_predictions.out"]
                for mut in rosetta_mutations:
                    expected_files += ["mut_%s_round_%d.pdb" % ("".join(str(x) for x in mut), x + 1) for x in
                                       range(iterations)]
            else:
                # What we normally get
                expected_files = ["ddg_predictions.out"]

            for expected_file in expected_files:
                if not os.path.exists(expected_file):
                    psb_status_manager.sys_exit_failure(
                        "Something went wrong with ddg_monomer, expected file %s is missing!")

            # These files are always empty
            # if os.path.isfile("wt_traj"):
            #    os.remove("wt_traj")

            # for mut in mutation_pose:
            #    if os.path.isfile("mutant_traj%s" % "".join(str(x) for x in mut)):
            #        os.remove("mutant_traj%s" % "".join(str(x) for x in mut))

            ddg_predictions_df = pd.read_csv(
                'ddg_prodictions.out',
                delim_whitespace=True)

            # Ultimate target is final_results file with everything we "want"
            final_results_df = pd.DataFrame(
                columns=["File", "Chain", "Residue", "WT_Res_Score", "Mutation", "ddG", "Protocol"])

            expected = len(self._mutations)

            if ddg_predictions_df.shape[0] != expected:
                psb_status_manager.sys_exit_failure(
                    "Something went wrong with ddg_monomer, expected %d preidctions, got only %d!" % (
                        expected, ddg_predictions_df.shape[0])
                )

            mutation_ddg = {}
            # temp = mutation_pose[:]
            # temp.sort(key=lambda x: x[1])
            # tempraw = "".join("".join(str(x) for x in y) for y in sorted(self._mutations, key=lambda z: z[1]))
            # Standard ddg uses the output ddG score as provided by application
            # Nonstandard ddg uses difference in pose scores between the mean score of the top 3 models per condition.
            if standard_ddg:
                for row in ddg_predictions_df.iterrows():
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
                        "WT_Res_Score": float(mutation_scores[original_residue]),
                        "Mutation": original_mutation,
                        "ddG": mutation_ddg[original_mutation],
                        "Protocol": protocol
                    })

            else:
                assert "NOn Standard DDG is not supported in 2020 see ddg_monomer.py"
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

            full_final_results_filename = os.path.join(save_current_directory, self.final_results_filename())
            final_results_df.to_csv(self.final_results_filename, sep='\t')
            LOGGER.info("Results saved to %s" % full_final_results_filename)
            return True, final_results_df

        success,final_results_df = run_part3_ddg()
        assert success,"run_par3_ddg() failed - not sure why it got to this point"

        self._log_close()
        os.chdir(save_current_directory)
        return True, final_results_df
