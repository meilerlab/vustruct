#!/usr/bin/env python
"""
class DDG_monomer manages all aspects of Rosetta ddg_monomer calculations 
inside the directory heirarchy supplied by a ddg_repo object
"""
import os
import pprint
import pandas as pd
from io import StringIO
from psb_shared.ddg_repo import DDG_repo
from psb_shared.ddg_base import DDG_base
from typing import Dict, List, Tuple, Union

import logging

LOGGER = logging.getLogger(__name__)


class DDG_monomer(DDG_base):
    @property
    def _mutation_list_filename(self):
        filename_prefix = "%s_%s_%s" % (
            self._ddg_repo.structure_id, self._ddg_repo.chain_id, self._mutationtext)
        return filename_prefix + ".mut"

    def _verify_applications_available(self):
        for application in [
            self._minimize_application_filename,
            self._per_residues_application_filename,
            self._ddg_monomer_application_filename]:
            # self._score_jd2_application_filename]:
            application_fullpath = os.path.join(self._ddg_repo.rosetta_bin_dir, application)
            assert os.path.exists(application_fullpath), (
                    "%s is a required application for ddg_monomer, but not found" % application_fullpath)

    def __init__(self, ddg_repo: DDG_repo, mutations: Union[str, List[str]]):
        """
        :param ddg_repo:
        :param mutations: S123AF format (or list of these) that indicate AA, res#, insert code, new AA.
        """

        self._mutations = None
        self._ddg_repo = None
        self.refresh_ddg_repo_and_mutations(ddg_repo, mutations)

        # Convert mutation to pose numbering for internal consistency
        # rosetta_mutations = []  # List of mutations input to ddg
        self._mutation_resids, self._rosetta_mutations, self._rosetta_residue_number_to_residue_xref = \
            self.mutations_to_rosetta_poses()

        #        self._root = root_path
        # Set the binary filenames in one place and make sure we have them before we start
        self._minimize_application_filename = 'minimize_with_cst.default.linuxgccrelease'
        self._per_residues_application_filename = 'per_residue_energies.default.linuxgccrelease'
        self._ddg_monomer_application_filename = 'ddg_monomer.default.linuxgccrelease'
        # self._score_jd2_application_filename =   'score_jd2.linuxgccrelease'

        self._application_timers = ['minimize', 'rescore', 'ddg_monomer']
        self._to_pdb = []

        LOGGER.debug("ddg_monomer.__init__ completed after setting:\n%s" % pprint.pformat(self.__dict__))

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

        self._ddg_repo.make_variant_directory_heirarchy()

        # Mirror the log - with INFO level details - to a local file

        # commandfile = "commands.txt"
        # Always use silent file/standard ddg for low quality models
        if not high_resolution:
            silent = True
            standard_ddg = True

        LOGGER.info('DDG_monomer.run called with\nsilent: %s\n,standard_ddg: %s\n,Row16: %s' % (
            silent, standard_ddg, high_resolution))



        save_current_directory = os.getcwd()
        os.chdir(self._ddg_repo.variant_dir)


        #############################################################################################
        # Part 1 of 4        Preminimize the cleaned and renumbered PDB
        #############################################################################################

        # We will generate file of harmonic restraints - an input to the part 3 calculation
        minimize_directory = os.path.join(self._ddg_repo.structure_dir, "minimize")

        ddg_monomer_atom_pair_constraints_basename = "%s_%s.cst" % (
            self._ddg_repo.structure_id, self._ddg_repo.chain_id)
        ddg_monomer_atom_pair_constraints_fullpath = os.path.join(
            minimize_directory, ddg_monomer_atom_pair_constraints_basename)

        def run_part1_minimizer() -> Tuple[str, List[str]]:
            """
            :return minimized_pdb filename and orignal atom pair distnaces for constraints or halt on errors:
            """
            # 2020-June-28 Chris Moth Command copied from
            # https://www.rosettacommons.org/docs/latest/application_documentation/analysis/ddg-monomer
            self._verify_applications_available()

            save_curwd = os.getcwd()
            previous_exit_code = 1

            # If the minimize directory is not there
            # Or somehow there without good exit code (should be impossible)
            # Then rerun
            def load_from_prior_minimize():
                LOGGER.info("os.chdir('%s')" % minimize_directory)
                os.chdir(minimize_directory)
                previous_exit_code, stdout, stdin = DDG_base._command_ran_previously(self._minimize_application_filename)
                assert previous_exit_code == 0, \
                    "%s directory lacks a 0 exit code recorded.  This should never happen" % \
                    minimize_directory
                assert os.path.isfile(ddg_monomer_atom_pair_constraints_basename), \
                    "%s directory lacks atom_pair constraints file %s" % (
                        minimize_directory, ddg_monomer_atom_pair_constraints_basename)
                LOGGER.info("Returning via os.chdir('%s')" % save_curwd)
                os.chdir(save_curwd)
                return previous_exit_code, stdout, stdin

            if os.path.isdir(minimize_directory):
                previous_exit_code, stdout, stderr = load_from_prior_minimize()

            minimized_pdb_filename = "minimized.%s_0001.pdb" % (
                os.path.basename(self._ddg_repo.cleaned_structure_filename).split(".")[0])

            if previous_exit_code == 0:
                LOGGER.info("Skipping minimization.  Using results from prior calculation")
                return_code = previous_exit_code
                runtime = 0.0
            else:
                # Then re-run the minimization in a subdirectory of the variant current directory
                # If all goes well, this directory will be moved to remove the minimize_directory
                tmp_directory = self._ddg_repo.mkdtemp(prefix='tmp_minimize',dir='.')
                LOGGER.info("os.chdir('%s')" % tmp_directory)
                os.chdir(tmp_directory)

                with os.fdopen(self._ddg_repo.os_open('minimize.log','w'),'w') as minimize_log_f:
                    minimize_log= logging.StreamHandler(minimize_log_f)
                    # minimize_log.setFormatter(logging.Formatter())
                    # logging.getLogger(__name__).addHandler(minimize_log)
                    minimize_log.setFormatter(logging.Formatter(DDG_base.log_format_string, DDG_base.date_format_string))
                    logging.getLogger().addHandler(minimize_log)


                    # Premin only takes a list of structs even if you only have 1
                    structure_list_filename = 'filename.input'
                    with os.fdopen(self._ddg_repo.os_open(structure_list_filename,'w'),'w') as structure_list_fp:
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

                    return_code, stdout, stderr, runtime = self._run_command_line_terminate_on_nonzero_exit(
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
                            "No AtomPair contraints returned from %s" % self._minimize_application_filename)

                    with os.fdopen(self._ddg_repo.os_open(ddg_monomer_atom_pair_constraints_basename, 'w'),'w') as outfile:
                        outfile.write("\n".join(original_atom_pair_distances_constraints))
                        outfile.write("\n")

                    LOGGER.info("%d atom pair constraints for ddg_monomer written to %s" % (
                        len(original_atom_pair_distances_constraints),
                        ddg_monomer_atom_pair_constraints_fullpath) )

                    logging.getLogger().removeHandler(minimize_log)

                LOGGER.info("Returning via os.chdir(%s)" % save_curwd)
                os.chdir(save_curwd)

                # Move our work to become the new "minimize" directory
                # so that other tasks do not have to repeat this work
                #
                # If OTHER tasks beat us to installing their work directory, then 
                # no worries!  Load that as the data source instead
                # and discard our hard work
                if return_code == 0:
                    rename_succeeded = self._attempt_directory_rename_delete_if_fail(tmp_directory, minimize_directory)

                if not rename_succeeded:
                    # Then some other task beat us...
                    # Load their outputs so all the ddgs are on the same page
                    return_code, stdout, stderr = load_from_prior_minimize()
                    assert return_code == 0

            _minimized_pdb_relpath = os.path.join('..', '..', minimize_directory, minimized_pdb_filename)
            if not os.path.isfile(_minimized_pdb_relpath):
                self._ddg_repo.psb_status_manager.sys_exit_failure('Minimized pdb file %s was not created by %s' % (
                    os.path.abspath(_minimized_pdb_relpath), self._minimize_application_filename))

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
            per_residues_directory = os.path.join(self._ddg_repo.structure_dir, "per_residue_energies")
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
                LOGGER.info("os.chdir('%s')" % per_residues_directory)
                os.chdir(per_residues_directory)
                previous_exit_code, stdout, stderr = DDG_base._command_ran_previously(self._per_residues_application_filename)
                assert previous_exit_code == 0, \
                    "%s directory lacks a 0 exit code recorded.  This should never happen" % \
                    per_residues_directory

                raw_scores = load_raw_per_residue_scores()
                LOGGER.info("Returning via os.chdir('%s')" % save_curwd)
                os.chdir(save_curwd)
                return previous_exit_code, stdout, stderr, raw_scores

            #######################################################################

            if os.path.isdir(per_residues_directory):
                previous_exit_code, stdout, stderr, raw_scores = load_from_prior_per_residues_run()
                returncode = previous_exit_code
                runtime = 0.0
            else:
                # Then re-run the minimization in a subdirectory of the variant current directory
                # If all goes well, this directory will be moved to remove the minimize_directory
                tmp_directory = self._ddg_repo.mkdtemp(prefix='tmp_per_residues', dir='.')
                os.chdir(tmp_directory)

                # Then re-run the minimization in a subdirectory of the variant current directory
                rescore_command = [
                    os.path.join(self._ddg_repo.rosetta_bin_dir, self._per_residues_application_filename),
                    '-s ' + minimized_pdb_relpath,
                    '-score:weights talaris2014',
                    '-out:file:silent %s' % min_residues_filename
                ]

                # stdout, stderr, runtime
                returncode, _, _, _ = self._run_command_line_terminate_on_nonzero_exit(
                    rescore_command,
                    additional_files_to_archive=[min_residues_filename])

                raw_scores = load_raw_per_residue_scores()

                LOGGER.info("Returning via os.chdir(%s)" % save_curwd)
                os.chdir(save_curwd)

                # Move our work to become the new "per_residues" directory
                # so that other tasks do not have to repeat this work
                #
                # If OTHER tasks beat us to installing their work directory, then 
                # no worries!  Load that as the data source instead
                # and discard our hard work
                LOGGER.info("Attempting os.rename('%s','%s')" % (tmp_directory, per_residues_directory))
                rename_succeeded = False
                if returncode == 0:
                    try:
                        os.rename(tmp_directory, per_residues_directory)
                        rename_succeeded = True
                    except OSError:
                        # The final per_residues directory already exists there
                        pass

                if not rename_succeeded:
                    # Then some other task beat us...
                    LOGGER.info("%s already installed by another process." % per_residues_directory)
                    LOGGER.info("Discarding current calculation and loading prior results")
                    # Load their outputs so all the ddgs are on the same page
                    returncode, stdout, stderr, raw_scores = load_from_prior_per_residues_run()
                    assert returncode == 0

            mutation_scores = {}
            for original_residue_id, rosetta_residue_id in self._ddg_repo.residue_to_clean_xref.items():
                # Extract the residue # from biopython resid tuple
                rosetta_residue_number = int(rosetta_residue_id[1])
                # Rosetta residues start with 1.  Scores indexed from 0
                raw_score = raw_scores[rosetta_residue_number - 1]
                if not raw_score:
                    self._ddg_repo.psb_status_manager.sys_exit_failure(
                        "Mutation residue %s result is None - investigate" % mutation)

                LOGGER.debug("Mutation score for %d %s=%s" % (
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
            previous_exit_code, stdout, stdin = DDG_base._command_ran_previously(self._ddg_monomer_application_filename)
            ddg_predictions_filename = 'ddg_predictions.out'

            # Does it look like we exited with no error before?

            previous_run_acceptable = False
            if previous_exit_code == 0:
                previous_run_acceptable = True

                if not os.path.isfile(ddg_predictions_filename):
                    previous_run_acceptable = False
                    LOGGER.warning("Prior ddg_monomer run missing %s.  Restarting" % \
                                   ddg_predictions_filename)

                if not os.path.isfile(self.final_results_filename()):
                    previous_run_acceptable = False
                    LOGGER.warning("Prior ddg_monomer run missing %s.  Restarting" % \
                                   self.final_results_filename())

            if previous_run_acceptable:
                # report back from prior run
                returncode = previous_exit_code
                final_results_df = pd.read_csv(self.final_results_filename(), sep='\t',
                                               dtype={'ddG': float, 'WT_Res_Score': float})
            else:
                # Start a new run of ddg_monomer
                # Get rid of the predictions output, else ddg appends to it - problems!
                try:
                    os.unlink(ddg_predictions_filename)
                except OSError:
                    pass

                with os.fdopen(self._ddg_repo.os_open(self._mutation_list_filename, 'w'),'w') as mutation_list_f:
                    mutation_list_f.write("total %d\n" % len(self._mutation_resids))
                    for rosetta_mutation in self._rosetta_mutations:
                        mutation_list_f.write("1\n")
                        rosetta_resno = int(rosetta_mutation[1:-1])

                        mutation_list_f.write("%s %d %s\n" % (
                            rosetta_mutation[0],
                            rosetta_resno,
                            rosetta_mutation[-1]))

                #
                # https://www.rosettacommons.org/docs/latest/application_documentation/analysis/ddg-monomer
                # STart with obtions common to both protocols
                ddg_options_filename = "%s_%s_%s.options" % (
                    self._ddg_repo.structure_id,
                    self._ddg_repo.chain_id,
                    self._mutationtext)

                ddg_options = [
                    '-in:file:s ' + minimized_pdb_relpath,
                    # PDB file of structure on which point mutations should be made
                    '-in::file::fullatom',  # read the input PDB file as a fullatom structure
                    '-ddg::mut_file ' + self._mutation_list_filename,  # the list of point mutations to consider in this run
                    '-fa_max_dis 9.0',  # optional -- if not given, the default value of 9.0 Angstroms is used.
                    '-ddg::dump_pdbs true',
                    # write one PDB for the wildtype and one for the pointmutant for each iteration
                    '-database ' + self._ddg_repo.rosetta_database_dir,  # the full oath to the database is required
                    '-ddg::iterations %d' % iterations,  # Typically 50 iterations of the algorithm
                    '-mute all',  # silence all of the log-file / stdout output generated by this protocol
                    '-ignore_unrecognized_res',  # ignore Rosetta-foreign res instead of quitting with error
                    '-ignore_zero_occupancy false', # Plenty of good CryoEM structures publish with 0 occupancy
                    '-ddg::suppress_checkpointing true',  # don't checkpoint LIZ DOES CHECKPOINTING WORK AT ALL?
                    '-ddg::output_silent true',  # write output to a silent file
                    '-ddg:weight_file soft_rep_design'
                    # soft-repulsive weights for initial sidechain optimization stage'
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

                        '-ddg::min_cst true',
                        # use distance restraints (constraints) during backbone minimization phase
                        # the set of constraints gleaned from the original (non-pre-relaxed/minimized) structure
                        '-constraints::cst_file ' + ddg_monomer_atom_pair_constraints_fullpath,
                        '-ddg::mean false',  # do not report the mean energy
                        '-ddg::min true',  # report the minimum energy
                        '-ddg::sc_min_only false',
                        # don't minimize only the backbone during backbone minimization phase

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

                self._move_prior_file_if_exists(ddg_options_filename)
                with os.fdopen(self._ddg_repo.os_open(ddg_options_filename, 'w'),'w') as ddg_options_f:
                    ddg_options_f.write("%s\n" % comment1)
                    ddg_options_f.write("%s\n" % comment2)
                    for ddg_option in ddg_options:
                        ddg_options_f.write("%s\n" % ddg_option)

                # Run actual ddg_monomer application

                # stdout, stderr, runtime =
                returncode, _, _, _ = self._run_command_line_terminate_on_nonzero_exit([
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
                    for mut in self._rosetta_mutations:
                        expected_files += ["mut_%s_round_%d.pdb" % ("".join(str(x) for x in mut), x + 1) for x in
                                           range(iterations)]
                else:
                    # What we normally get
                    expected_files = [ddg_predictions_filename]

                for expected_file in expected_files:
                    if not os.path.exists(expected_file):
                        self._ddg_repo.psb_status_manager.sys_exit_failure(
                            "Something went wrong with ddg_monomer, expected file %s is missing!" % expected_file)

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
                    for index, row in ddg_predictions_df.iterrows():
                        rosetta_mutation = row['description']

                        # Populate the ddgs back to the original structure positions (perhaps with insertion codes)

                        original_residue = self._rosetta_residue_number_to_residue_xref[int(rosetta_mutation[1:-1])]
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
                        }, ignore_index=True)

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

                with os.fdopen(self._ddg_repo.os_open(self.final_results_filename(),'w'),'w') as frf:
                    final_results_df.to_csv(frf, sep='\t')
                LOGGER.info("Results saved to %s" % os.path.abspath(self.final_results_filename()))
                # End of ddg_monomer calculation
            return True, final_results_df

        self._ddg_repo.psb_status_manager.write_info("Part 3: DDG monomer")
        success, final_results_df = run_part3_ddg()
        if not success:
            self._ddg_repo.psb_status_manager.sys_exit_failure(
                "run_part3_ddg() failed - not sure why it got to this point")

        self._ddg_repo.psb_status_manager.mark_complete()

        # self._log_close()
        os.chdir(save_current_directory)
        return True, final_results_df

    def jobstatus(self) -> Dict[str, str]:
        jobstatus_info = {}
        save_current_directory = os.getcwd()
        try:
            os.chdir(self._ddg_repo.variant_dir)
        except FileNotFoundError:
            jobstatus_info['ExitCode'] = None
            jobstatus_info['jobinfo'] = 'ddg_monomer has not created %s'%self._ddg_repo.variant_dir
            jobstatus_info['jobprogress'] = jobstatus_info['jobinfo']
            return jobstatus_info

        _, _, exitcd_filename = self._command_result_filenames(self._ddg_monomer_application_filename)
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
        """Return the calculated ddG results as a pandas series, or None"""
        _,_,exitcd_filename = self._command_result_filenames(self._ddg_monomer_application_filename)
        save_current_directory = os.getcwd()
        try:
            os.chdir(self._ddg_repo.variant_dir)
        except FileNotFoundError:
            LOGGER.info("No results directory %s",self._ddg_repo.variant_dir)
            return None

        previous_exit, previous_exit_int = self._get_previous_exit(exitcd_filename)

        final_results_fullpath = os.path.join(self._ddg_repo.variant_dir,self.final_results_filename())

        if previous_exit_int == 0 and os.path.exists(final_results_fullpath):
            final_results_df = pd.read_csv(final_results_fullpath, sep='\t', dtype={'ddG': float, 'WT_Res_Score': float})
            assert len(final_results_df) == 1,"Format error: single row of data not found in %s - halting"%final_results_fullpath
            final_results_df['RefAA'] = final_results_df['Mutation'].str[0]
            final_results_df['AltAA'] = final_results_df['Mutation'].str[-1]
            LOGGER.debug("ddg of %f read from %s" % (final_results_df['ddG'], final_results_fullpath))
        else:
            LOGGER.info("ddg results file %s not found" % (final_results_fullpath))
            final_results_df = None

        os.chdir(save_current_directory)
        return final_results_df.iloc[0] if final_results_df is not None else None
