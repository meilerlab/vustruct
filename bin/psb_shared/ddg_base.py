#!/usr/bin/env python
"""
class DDG_base underpins both DDG_monomer and DDG_cartesian with
a variety of methods shared by bost.
"""
import os
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

    def final_results_filename(self) -> str:
        return "%s_%s_%s.csv" % (
            self._ddg_repo.structure_id,
            self._ddg_repo.chain_id,
            self._mutationtext
        )

    def _command_result_filenames(self, binary_program_basename: str):
        return (
            binary_program_basename + ".stdout",
            binary_program_basename + ".stderr",
            binary_program_basename + ".exit")

    def _get_previous_exit(self, exitcd_filename: str) -> Tuple[str, int]:
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

