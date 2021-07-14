#!/usr/bin/env python
"""
class DDG_cartesian manages all aspects of Rosetta ddg_cartesian calculations 
inside the directory heirarchy supplied by a ddg_repo object

Perform a rosetta ddg_cartesian calculation per the guidelines at
https://new.rosettacommons.org/docs/latest/cartesian-ddG
supplemented by 2 critically important papers.

Paper 1 of 2: 2016
https://doi.org/10.1021/acs.jctc.6b00819
Hahnbeom Park, et al. David Baker, and Frank DiMaio
Simultaneous optimization of biomolecular energy function on features from small molecules and macromolecules"
JCTC.

Paper 2 of 2: 2020 - Critical improvements:
https://doi.org/10.3389/fbioe.2020.558247
Frenz, Frank DiMaio, et al. Yifan Song. 2020.
Prediction of Protein Mutational Free Energy:
Benchmark and Sampling Improvements Increase Classification Accuracy.
Frontiers in Bioengineering and Biotechnology 8 (October): 558247.

:return: (True if success,  dataframe of final results)
"""
import os
import sys
import datetime
import lzma
import pprint
from collections import OrderedDict
import subprocess
import tempfile
import time
import pandas as pd
import numpy as np
from io import StringIO
from psb_shared.ddg_repo import DDG_repo
from psb_shared.ddg_base import DDG_base
from typing import Dict, List, Tuple, Union
from psb_shared.psb_os import pushdir

import logging

LOGGER = logging.getLogger(__name__)


class DDG_cartesian(DDG_base):

    # DDG_base defines
    # def refresh_ddg_repo_and_mutations(self, ddg_repo: DDG_repo, mutations: Union[str, List[str]]):

    def _verify_applications_available(self):
        for application in [
            self._relax_application_filename,
            self._ddg_cartesian_application_filename]:
            # self._score_jd2_application_filename]:
            application_fullpath = os.path.join(self._ddg_repo.rosetta_bin_dir, application)
            assert os.path.exists(application_fullpath), (
                    "%s is a required application for ddg_cartesian, but not found" % application_fullpath)

    def __init__(self, ddg_repo: DDG_repo, mutations: Union[str, List[str]]):
        """
        :param ddg_repo:
        :param mutations: S123AF format (or list of these) that indicate AA, res#, insert code, new AA.
        """

        # Let the base class init the repository self._ddg_repo and self._mutations
        super(DDG_cartesian, self).__init__(ddg_repo,mutations)

        self.refresh_ddg_repo_and_mutations(ddg_repo, mutations)

        #        self._root = root_path
        # Set the binary filenames in one place and make sure we have them before we start
        # self._relax_application_filename = 'relax.default.linuxgccrelease'
        self._relax_application_filename = 'rosetta_scripts.default.linuxgccrelease'
        self._ddg_cartesian_application_filename = 'cartesian_ddg.default.linuxgccrelease'
        # self._ddg_cartesian_application_filename = 'rosetta_scripts.default.linuxgccrelease'

        self._to_pdb = []

        # Convert mutation to pose numbering for internal consistency
        # rosetta_mutations = []  # List of mutations input to ddg
        self._mutation_resids, self._rosetta_mutations, self._rosetta_residue_number_to_residue_xref = \
            self.mutations_to_rosetta_poses()

        if not self._rosetta_mutations:
           self._ddg_repo.psb_status_manager.sys_exit_failure(
               "No rosetta numbered mutations found to match your request of %s" % self._mutations)

        if len(self._rosetta_mutations) != 1:
            self._ddg_repo.psb_status_manager.sys_exit_failure(
                "Only 1 mutation can be specified in __init__ parameter mutations.  You supplied: %s" % self._mutations)


        # Dump a LOT of stuff.
        LOGGER.debug("ddg_cartesian.__init__ completed after setting:\n%s" % pprint.pformat(self.__dict__))

    @property
    def _ddg_predictions_filename(self):
        filename_prefix = "%s_%s_%s" % (
            self._ddg_repo.structure_id, self._ddg_repo.chain_id, self._mutationtext)
        return filename_prefix + ".ddg"

    @property
    def _mutation_list_filename(self):
        filename_prefix = "%s_%s_%s" % (
            self._ddg_repo.structure_id, self._ddg_repo.chain_id, self._mutationtext)
        return filename_prefix + ".mut"

    @property
    def _ddg_cartesian_final_directory(self):
        """
        return the ....variant/ddg_cartesian where final results should be found, if available
        """
        return os.path.join(self._ddg_repo.variant_dir,"ddg_cartesian")

    @property
    def analyze_cartesian_ddg(self) -> pd.DataFrame:
        """
        Closely copied from Bian Li's analyse_cartesian.py script
        This class function opens the weird Rosetta ddg_cartesian output file, parses the energy/value pairs
        and returns the ddgs
        """

        # Frenz et al 2020:
        # In the analysis step, we changed the number of mutant
        # models generated using the following convergence criterion: the
        # lowest energy 2 structures must converge to within 1 Rosetta
        # Energy Unit, or take the best of 5 models, whichever comes first."""

        # parse the cartesian_ddg results
        # The format of this file is quite annoying.  It combines WT and MUT energy scores in one file.
        """  Below is first 80 columns of a .ddg file output\
COMPLEX:   Round1: WT:  -623.944  fa_atr: -1146.164 fa_rep:   147.805 fa_sol:.......   
COMPLEX:   Round2: WT:  -623.951  fa_atr: -1146.164 fa_rep:   147.756 fa_sol:   
COMPLEX:   Round3: WT:  -623.946  fa_atr: -1146.181 fa_rep:   147.800 fa_sol:   
COMPLEX:   Round1: MUT_106VAL:  -616.753  fa_atr: -1141.595 fa_rep:   148.003 fa
COMPLEX:   Round2: MUT_106VAL:  -616.754  fa_atr: -1141.586 fa_rep:   147.922 fa
COMPLEX:   Round3: MUT_106VAL:  -615.593  fa_atr: -1141.217 fa_rep:   147.892 fa
        """


        wt_energy_dictionaries = []
        mut_energy_dictionaries = []

        LOGGER.info("analyze_cartesian_ddg occuring in directory %s", self._ddg_cartesian_final_directory)

        ddg_predictions_filename_fullpath = os.path.join(
            self._ddg_cartesian_final_directory,self._ddg_predictions_filename
        )
        if not os.path.exists(ddg_predictions_filename_fullpath):
            LOGGER.critical("Unable to open %s - no results available",ddg_predictions_filename_fullpath)
            return None

        with open(ddg_predictions_filename_fullpath, 'rt') as ddg_pf:
            for line in ddg_pf:
                if not line.startswith("COMPLEX: "):
                    message = "ddg_cartesian output file %s:%s lacks required COMPLEX: start" % (
                        self._ddg_predictions_filename,
                        line)
                    LOGGER.critial(message)
                    sys.exit(message)

                # But we make sure that they are formatted properly
                energy_fields = line.strip().split()
                if len(energy_fields) < 4:
                    if not line.startswith("COMPLEX: "):
                        message = "ddg_cartesian output file %s:%s has insufficiennt columns" % (
                            self._ddg_predictions_filename,
                            line)

                    LOGGER.critial(message)
                    sys.exit(message)

                # We now split the WT and MUT..
                # split out the energy description words, and their associated energie floating point values
                # We ignore the first two columns which are 0=COMPLEX: and 1=RoundN:
                energy_words = [energy_description.strip(':') for energy_description in energy_fields[2::2]]
                try:
                    energy_values = [float(energy_value) for energy_value in energy_fields[3::2]]
                except ValueError:
                    message = "ddg_cartesian output file %s:%s has energy values that are not float" % (
                        self._ddg_predictions_filename,
                        line)

                    LOGGER.critial(message)
                    sys.exit(message)


                wt_or_mut_energy_dictionaries = None
                if energy_words[0].startswith('WT_'):  # Add to the wild type (starting) energies dataframe
                    wt_or_mut_energy_dictionaries = wt_energy_dictionaries
                elif energy_words[0].startswith('MUT_'):  # Add to the variant (end point) energies dataframe
                    wt_or_mut_energy_dictionaries = mut_energy_dictionaries
                else:
                    message = "ddg_cartesian output file %s:%s has neither WT nor MUT_ energies" % (
                        self._ddg_predictions_filename,
                        line)
                    LOGGER.critical(message)
                    sys.exit(message)

                if len(energy_words) != len(energy_values):
                    message = "ddg_cartesian output file %s:%s lacks poairing of energy words and values" % (
                        self._ddg_predictions_filename,
                        line)
                    LOGGER.critial(message)
                    sys.exit(message)

                # Replace the WT or MUT_..: with "total" in the new dataframes
                energy_words[0] = 'total'
                wt_or_mut_energy_dictionaries.append(OrderedDict(zip(energy_words, energy_values)))

        LOGGER.info("Parsed %d WT rows and %d MUT rows from %s", len(wt_energy_dictionaries),len(mut_energy_dictionaries),self._ddg_predictions_filename)

        # calculate the mean and standard deviations for each energy term
        wt_energy_df = pd.DataFrame(wt_energy_dictionaries)
        wt_energy_means = wt_energy_df.mean(axis=0)
        wt_energy_stds = wt_energy_df.std(axis=0)

        mut_energy_df = pd.DataFrame(mut_energy_dictionaries)
        mut_energy_2_lowest = mut_energy_df.nsmallest(2,'total') # Do NOT use 'all' as you can get more than 2 rows!!

        lowest_2_total_energies = mut_energy_2_lowest['total'].to_numpy()
        abs_mut_energy_2_lowest_diff = abs(lowest_2_total_energies[0] - lowest_2_total_energies[1])

        if abs_mut_energy_2_lowest_diff > 1:
            LOGGER.warning("DDG Cartesian min energy mutant calculations did not converge to within 1 Rosetta Energy Unit.  Difference=%f",
                           abs_mut_energy_2_lowest_diff)
        else:
            LOGGER.warning("DDG Cartesian min energy mutant calculations converged to %f Rosetta Enerfy Units",
                           abs_mut_energy_2_lowest_diff)

        # mut_energy_means = mut_energy_2_lowest.mean(axis=0)

        wt_2_rows_matching_mut_2_lowest_df = wt_energy_df.iloc[mut_energy_2_lowest.index]

        # mut_energy_stds = mut_energy_2_lowest.std(axis=0)

        ddgs_df = pd.DataFrame((mut_energy_2_lowest - wt_2_rows_matching_mut_2_lowest_df).mean()).transpose()

        LOGGER.debug("DDGs Calculated from %s\n%s",self._ddg_predictions_filename,str(ddgs_df))

        return ddgs_df

    def run(self) -> Tuple[bool, pd.DataFrame]:

        self._ddg_repo.make_variant_directory_heirarchy()

        # Mirror the log - with INFO level details - to a local file

        # commandfile = "commands.txt"
        # Always use silent file/standard ddg for low quality models
        #if not high_resolution:
        #    silent = True
        #    standard_ddg = True

        LOGGER.info('DDG_cartesian.run() starting')
        # called with\nsilent: %s\n,standard_ddg: %s\n,Row16: %s' % (
        #    silent, standard_ddg, high_resolution))

        timestamp_suffix = datetime.datetime.now().strftime("%d%m%y%H%M%S")
        # isoresult, isopdbname, to_pdb = isolate_chain(self._pdb, self._chain)
        # if not isoresult:
        #    return self.failure(None, to_pdb)
        # else:
        #    self._to_pdb = to_pdb

        save_current_directory = os.getcwd()
        os.chdir(self._ddg_repo.variant_dir)

        #############################################################################################
        # Part 1 of 2        Relax the cleaned a PDB
        #############################################################################################

        relax_directory = os.path.join(self._ddg_repo.structure_dir, "relax")

        def run_part1_relax() -> Tuple[str, List[str]]:
            """
            Relax creates a number of candidate 'relaxed' .pdb files.  In part 2, we select
            lowest energy conformation for ddg cartesian calculation
            """
            # See 2020 Frenz et al for guidance
            self._verify_applications_available()

            save_curwd = os.getcwd()
            previous_exit_code = 1

            # If the minimize directory is not there
            # Or somehow there without good exit code (should be impossible)
            # Then rerun
            def load_from_prior_relax():
                with pushdir(relax_directory):
                    previous_exit_code, stdout, stdin = self._command_ran_previously(self._relax_application_filename)
                    assert previous_exit_code == 0, \
                        "%s directory lacks a 0 exit code recorded.  This should never happen" % \
                        relax_directory
                    return previous_exit_code, stdout, stdin

            if os.path.isdir(relax_directory):
                previous_exit_code, stdout, stderr = load_from_prior_relax()

            target_score_filename = "target.score"

            if previous_exit_code == 0:
                LOGGER.info("Skipping relax.  Using results from prior calculation")
                returncode = previous_exit_code
                runtime = 0.0
            else:
                # Then re-run the minimization in a subdirectory of the variant current directory
                # If all goes well, this directory will be moved to remove the relax_directory
                tmp_directory = tempfile.mkdtemp(prefix='tmp_relax',dir='.')

                os.makedirs(tmp_directory, mode=0o770, exist_ok=True)
                LOGGER.info("os.chdir('%s')" % tmp_directory)
                os.chdir(tmp_directory)
 
                
                # relax_script_filename = 'cartesian_relax.script'
                # with open(relax_script_filename, mode='w') as relax_script_fp:
                #    relax_script_fp.write(
                """\
switch:cartesian
repeat 2
ramp_repack_min 0.02  0.01     1.0  50
ramp_repack_min 0.250 0.01     0.5  50
ramp_repack_min 0.550 0.01     0.0 100
ramp_repack_min 1     0.00001  0.0 200
accept_to_best
endrepeat"""
                #)

                # Taken directly from 2020 Frenz et all supplement information
                with open("cartesianrelaxprep.xml", mode='w') as relax_script_fp:
                    relax_script_fp.write("""\
<ROSETTASCRIPTS>
  <SCOREFXNS>
    <ScoreFunction name="fullatom" weights="ref2015_cart" symmetric="0" />
  </SCOREFXNS>
  <MOVERS>
    <FastRelax name="fastrelax" scorefxn="fullatom" cartesian="1" repeats="4" />
  </MOVERS>
  <PROTOCOLS>
    <Add mover="fastrelax"/>
  </PROTOCOLS>
<OUTPUT scorefxn="fullatom"/>
</ROSETTASCRIPTS>
""")

                not_used_old_idea_relax_command = [
                    os.path.join(self._ddg_repo.rosetta_bin_dir, self._relax_application_filename),
                    "-in:file:s %s"%format(self._ddg_repo.cleaned_structure_filename),
                    # Not sure about this one: "-in:file:native %s"%format(self._ddg_repo.cleaned_structure_filename),
                    "-in:file:fullatom",
                    "-nstruct 20",
                    "-ignore_unrecognized_res",
                    "-ignore_zero_occupancy false",
                    "-parser:protocol cartesianrelaxprep.xml",
                    "-out:file:scorefile target.score",
                    "-score:set_weights cart_bonded 0.5 pro_close 0",
                    "-relax:cartesian true",
                    "-multiple_processes_writing_to_one_directory",
                    "-no_color"]

                    # "-relax:cartesian-score:weights ref2015_cart",
                    # "-relax:min_type lbfgs_armijo_nonmonotone",
                    # "-relax_script card2.script",
                    # Note from colleages: # modify fa_atr and fa_sol
                    # behavior, really important for protein stability (default: 6).
                    # This flag needs to match what is used in the cartesian ddg options below.
                    # "-fa_max_dis 9.0"
                    #]

                # Try the new approach using Mover as described in 2020 paper
                relax_command = [
                    os.path.join(self._ddg_repo.rosetta_bin_dir, self._relax_application_filename),
                    "-in:file:s %s"%format(self._ddg_repo.cleaned_structure_filename),
                    # Not sure about this one: "-in:file:native %s"%format(self._ddg_repo.cleaned_structure_filename),
                    "-in:file:fullatom",
                    "-default_max_cycles 200",
                    "-ignore_unrecognized_res",
                    "-ignore_zero_occupancy false",
                    "-parser:protocol cartesianrelaxprep.xml",
                    "-out:file:scorefile target.score",
                    "-missing_density_to_jump",
                    "-nstruct 20",
                    "-fa_max_dis 9",
                    "-relax:cartesian true",
                    "-multiple_processes_writing_to_one_directory"]


                returncode, stdout, stderr, runtime = self._run_command_line_terminate_on_nonzero_exit(
                    relax_command)
                    # additional_files_to_archive=[minimized_pdb_filename])

                LOGGER.info("Returning via os.chdir(%s)" % save_curwd)
                os.chdir(save_curwd)

                # Move our work to become the new "minimize" directory
                # so that other tasks do not have to repeat this work
                #
                # If OTHER tasks beat us to installing their work directory, then 
                # no worries!  Load that as the data source instead
                # and discard our hard work
                LOGGER.info("Attempting os.rename(%s,%s)" % (tmp_directory, relax_directory))
                rename_succeeded = False
                if returncode == 0:
                    try:
                        os.rename(tmp_directory, relax_directory)
                        rename_succeeded = True
                    except OSError:
                        # The final minimize directory already exists there
                        pass

                if not rename_succeeded:
                    # Then some other task beat us...
                    LOGGER.info("%s already installed by another process." % relax_directory)
                    LOGGER.info("Discarding current calculation and loading prior results")
                    # Load their outputs so all the ddgs are on the same page
                    returncode, stdout, stderr = load_from_prior_relax()
                    assert returncode == 0

            _target_score_relpath = os.path.join('..', '..', relax_directory, target_score_filename)
            if not os.path.isfile(_target_score_relpath):
                self._ddg_repo.psb_status_manager.sys_exit_failure('file %s was not created by %s' % (
                    os.path.abspath(_target_score_relpath), self._relax_application_filename))

            return _target_score_relpath

        self._ddg_repo.psb_status_manager.write_info("Part 1: Relax")
        target_score_relpath = run_part1_relax()
        target_score_df = pd.read_csv(target_score_relpath,skiprows=1,sep='\s+')


        LOGGER.info("%s with %d scores available for part 2",target_score_relpath,len(target_score_df))

        #############################################################################################
        # Part 2 of 2       Determine the lowest energy relax structure and run Cartesian
        #############################################################################################
        def run_part2_ddg_cartesian():
            previous_exit_code, stdout, stdin = self._command_ran_previously(self._ddg_cartesian_application_filename)


            # Step 2.1 is to identify the lowest energy structure from the relax operation.
            # Grab the .pdb filename based on the minimum score in teh dataframe.
            min_total_score_rows = target_score_df[target_score_df.total_score == target_score_df.total_score.min()]
            min_scoring_relaxed_pdb = os.path.join("../../../relax/",min_total_score_rows.iloc[0].description + ".pdb")

            self._verify_applications_available()
            save_curwd = os.getcwd()
            previous_exit_code = 1
            stderr = ""
            stdout = ""



            # If the ddg_cartesian directory is not there
            # Or somehow there without good exit code (should be impossible)
            # Then rerun
            def load_from_prior_ddg_cartesian_run():
                with pushdir(self._ddg_cartesian_final_directory):
                    previous_exit_code, stdout, stderr = self._command_ran_previously(self._ddg_cartesian_application_filename)
                    assert previous_exit_code == 0, \
                        "%s directory lacks a 0 exit code recorded.  This should never happen" % \
                        self._ddg_cartesian_final_directory

                    return previous_exit_code, stdout, stderr

            #######################################################################

            if os.path.isdir(self._ddg_cartesian_final_directory):
                previous_exit_code, stdout, stderr = load_from_prior_ddg_cartesian_run()
                returncode = previous_exit_code
                runtime = 0.0
            else:
                # We need to run the cartesian calculation.
                tmp_directory = tempfile.mkdtemp(prefix='tmp_ddg_cartesian', dir='.')
                os.makedirs(tmp_directory, mode=0o770, exist_ok=True)
                LOGGER.info("os.chdir('%s')" % tmp_directory)
                os.chdir(tmp_directory)

               
                # Create a list of mutations to be analyzed by ddg_monomer, .
                # Documented https://www.rosettacommons.org/docs/latest/application_documentation/analysis/ddg-monomer
                # Cartesian uses same format as monomer for th e".mut" file
                optimize_proline_flag = 'false'
                with open(self._mutation_list_filename, 'w') as mutation_list_f:
                    mutation_list_f.write("total %d\n" % len(self._mutation_resids))
                    for rosetta_mutation in self._rosetta_mutations:
                        mutation_list_f.write("1\n")
                        rosetta_resno = int(rosetta_mutation[1:-1])

                        mutation_list_f.write("%s %d %s\n" % (
                            rosetta_mutation[0],
                            rosetta_resno,
                            rosetta_mutation[-1]))

                        if rosetta_mutation[0] == 'P' or rosetta_mutation[-1] == 'P':
                            optimize_proline_flag = 'true'


                # ddg_cartesian_command = [
                #     os.path.join(self._ddg_repo.rosetta_bin_dir, self._ddg_cartesian_application_filename),
                #     "-database %s"%self._ddg_repo.rosetta_database_dir,
                #     "-in:file:s %s"%min_scoring_relaxed_pdb,       # The lowest scoring pdb produced by relax
                #     '-ddg::mut_file ' + self._mutation_list_filename,  # the list of point mutations to consider in this run
                #     "-score:weights ref2015_cart",
                #     "-ddg:iterations 10",
                #     "-fa_max_dis 9.0",
                #     "-ddg:dump_pdbs true",
                #     "-ddg:cartesian",
                #     "-ddg:bbnbrs 1",
                #     "-detect_disulf false",
                #     "-fa_max_dis 9.0"
                #     ]

                # min_scoring_relaxed_pdb = "/dors/capra_lab/users/mothcw/ddg_cartesian/FrenzReconstruction/pdbs/%s_%s.pdb"%(
                #    self._ddg_repo.structure_id.upper(),self._ddg_repo.chain_id)
                # LOGGER.warning("DO NOT CHECK THIS ABOVE IN OVERRIDE OF STRUCT TO %s",min_scoring_relaxed_pdb)
                ddg_cartesian_command = [
                    os.path.join(self._ddg_repo.rosetta_bin_dir, self._ddg_cartesian_application_filename),
                    "-database %s"%self._ddg_repo.rosetta_database_dir,
                    "-s %s"%min_scoring_relaxed_pdb,       # The lowest scoring pdb produced by relax
                    "-ddg::iterations 3",
                    "-ddg::score_cutoff 1",
                    "-ddg:dump_pdbs true",
                    "-ddg::bbnbrs 1",
                    "-score:weights ref2015_cart",
                    '-ddg::mut_file ' + self._mutation_list_filename,  # the list of point mutations to consider in this run
                    "-ddg:frag_nbrs 2",
                    "-ignore_zero_occupancy false",
                    "-missing_density_to_jump",
                    "-ddg:flex_bb false",
                    "-ddg::force_iterations false",
                    "-fa_max_dis 9.0",
                    # "-ddg::json true",
                    "-ddg:optimize_wt true",
                    "-ddg:optimize_proline %s"%optimize_proline_flag, # Set above if WT or MUT res is Pro
                    "-ddg:legacy false"
                    # "-ddg:cartesian",
                    # "-detect_disulf false",
                    ]







                # stdout, stderr, runtime
                returncode, _, _, _ = self._run_command_line_terminate_on_nonzero_exit(
                    ddg_cartesian_command,
                    additional_files_to_archive=["Need to list there here"])

                LOGGER.info("Returning via os.chdir(%s)" % save_curwd)
                os.chdir(save_curwd)

                # Move our work to become the new "ddg_cartesian" directory
                # so that other tasks do not have to repeat this work
                #
                # If OTHER tasks beat us to installing their work directory, then 
                # no worries!  Load that as the data source instead
                # and discard our hard work
                LOGGER.info("Attempting os.rename('%s','%s')" % (tmp_directory, self._ddg_cartesian_final_directory))
                rename_succeeded = False
                if returncode == 0:
                    try:
                        os.rename(tmp_directory, self._ddg_cartesian_final_directory)
                        rename_succeeded = True
                    except OSError:
                        # The final ddg_cartesian directory already exists there
                        pass

                if not rename_succeeded:
                    # Then some other task beat us...
                    LOGGER.info("%s already installed by another process." % self._ddg_cartesian_final_directory)
                    LOGGER.info("Discarding current calculation and loading prior results")
                    # Load their outputs so all the ddgs are on the same page
                    returncode, stdout, stderr, raw_scores = load_from_prior_ddg_cartesian_run()
                    assert returncode == 0

            # end of part 2 - nothing is returned


        self._ddg_repo.psb_status_manager.write_info("Part 2: Cartesian ddg")
        run_part2_ddg_cartesian()

        return 0,self.analyze_cartesian_ddg;

    def jobstatus(self) -> Dict[str, str]:
        jobstatus_info = {}
        save_current_directory = os.getcwd()
        try:
            os.chdir(self._ddg_cartesian_final_directory)
        except FileNotFoundError:
            jobstatus_info['ExitCode'] = None
            jobstatus_info['jobinfo'] = 'ddg_cartesian has not created %s' % self._ddg_cartesian_final_directory
            jobstatus_info['jobprogress'] = jobstatus_info['jobinfo']
            return jobstatus_info

        _, _, exitcd_filename = self._command_result_filenames(self._ddg_cartesian_application_filename)
        previous_exit, previous_exit_int = self._get_previous_exit(exitcd_filename)
        if previous_exit_int == 0:
            jobstatus_info['ExitCode'] = '0'
            jobstatus_info[
                'jobprogress'] = 'Completed'  # For now, we will overwrite memory on happy exit - need to think
            jobstatus_info['jobinfo'] = 'Completed'
        else:
            jobstatus_info['ExitCode'] = previous_exit
            info, progress = self._ddg_repo.psb_status_manager.read_info_progress()
            jobstatus_info['jobinfo'] = info
            jobstatus_info['jobprogress'] = progress

        os.chdir(save_current_directory)
        return jobstatus_info

    def retrieve_result(self) -> pd.Series:
        """
        Return the ddg_cartesian results details for the one variant as a Pandas Series
        """
        save_current_directory = os.getcwd()
        try:
            os.chdir(self._ddg_repo.variant_dir)
        except FileNotFoundError:
            LOGGER.info("No results directory %s",self._ddg_repo.variant_dir)
            return None

        ddg_cartesian_result_df = None

        _,_,exitcd_filename = self._command_result_filenames(self._ddg_cartesian_application_filename)
        exitcd_fullpath = os.path.join(self._ddg_cartesian_final_directory,exitcd_filename)
        previous_exit,previous_exit_int = self._get_previous_exit(exitcd_fullpath)

        if previous_exit is None:
            LOGGER.info("ddg Results exit status file %s not found" % exitcd_fullpath)
        elif previous_exit_int == 0:
            ddg_cartesian_result_df = self.analyze_cartesian_ddg
            # Above returns None if results could not be loaded.  Put variant at front of series
            if isinstance(ddg_cartesian_result_df,pd.DataFrame):
                ddg_cartesian_result_df['structure_id'] = self._ddg_repo.structure_id
                ddg_cartesian_result_df['chain'] = self._ddg_repo.chain_id
                ddg_cartesian_result_df['variant'] = self._mutationtext
        else:
            LOGGER.info("ddg Results exit status file %s contains non-zero exit failure: %s" % (exitcd_fullpath,previous_exit))

        os.chdir(save_current_directory)

        return ddg_cartesian_result_df.iloc[0] if ddg_cartesian_result_df is not None else None
