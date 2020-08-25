#!/usr/bin/env python
"""
Calculation-agnostic class to provide directory locations 
and key filenames shared by ddg calcualtions and for 
results retrieval.
The heirarchy is setup by __init__ parameters 
AND a subsequent call to set structure ID and chain
AND finally setting a variant to finalize a calculation directory

Restated:
1) Construct a DDG_repo() object  (calls __init__())
2) Call ddg_repo.set_pdb()/swiss/etc to specify a type of input structure
3) Call ddg_repo.set_variant() which sets final repo calculation directory details
"""
import sys
import os
import re
from typing import Dict, List, Tuple, Union

import logging
import json
import tempfile
import datetime
import configparser
from psb_shared.psb_progress import PsbStatusManager
LOGGER = logging.getLogger(__name__)


class DDG_repo():

    @staticmethod
    def ddg_config_read(ddg_config_filename: str) -> Dict[str, str]:
        LOGGER.info("Attempting to read ddg_config from %s", ddg_config_filename)

        ddg_config = configparser.ConfigParser(allow_no_value=False)

        parsed_file_list = ddg_config.read(ddg_config_filename)
        if not parsed_file_list:
            exit_str = "TERMINATING - Unable to read ddg config file %s" % ddg_config_filename
            LOGGER.critical(exit_str)
            sys.exit(exit_str)

        LOGGER.info("Successfully read ddg_config from " + ddg_config_filename)
        required_ddg_config_items = ['repo_root', 'rosetta_version', 'rosetta_bin_dir']

        assert ddg_config.has_section('ddG'),'The ddg_config file %s lacks a [ddG] section'%ddg_config_filename
        ddg_config_dict = dict(ddg_config.items("ddG"))
        # print('rosetta_bin_dir=',ddg_config_dict['rosetta_bin_dir'])

        for required_item in required_ddg_config_items:
            if not required_item in ddg_config_dict:
                exit_str = "TERMINATING: %s key must be specified in ddg_config file %s" % (
                        required_item, ddg_config_filename)
                LOGGER.critical(exit_str)
                sys.exit(exit_str)

        return ddg_config_dict


    _my_config_dicts = {} # Global dict of read-in ddg configuration files

    def __init__(self, 
            ddg_config_filename:str, 
            calculation_flavor: str):
        """
        Begin process of constructing a ddg repository manager.
        Call a set_pdb/swiss/etc function after contruction.
        Then, specify a variant to complete final initialization

        :param ddg_config_filename Optional if repo_root_dir and rosetta_version specified
         Typically caller finds in ddg repo's config file
        :param calculation_flavor: Ex: "ddg_monomer"  Must be set by caller.
        """

        if ddg_config_filename in DDG_repo._my_config_dicts:
            LOGGER.info("Re-using config from prior read of %s"%ddg_config_filename)
            self._ddg_config_dict = DDG_repo._my_config_dicts[ddg_config_filename]
        else:
            self._ddg_config_dict = DDG_repo.ddg_config_read(ddg_config_filename)
            DDG_repo._my_config_dicts[ddg_config_filename] = self._ddg_config_dict

        self._ddg_config_filename = ddg_config_filename

        self._structure_dir = None  # parent directory for all calculations on a structure.  Contains cleaned PDB and xref
        self._calculation_dir = None  # final directory where all calculucations are conducted
        self._residue_to_clean_xref_filename = None  # Json file to Map original structure residue IDs to the rosetta-ready residue ids
        self._cleaned_structure_pdb_filename = None
        self._structure_config_filename = None
        self._calculation_flavor = calculation_flavor
        self._structure_source = None
        self._residue_to_clean_xref = {}
        self._structure_config = None
        self._structure_id = None
        self._slurm_dir = None
        self._chain_id = None
        self._variant = None
        self._log_filename = None
        self._psb_status_manager = None

        self._ddg_root = os.path.join(self.repo_root_dir, self._calculation_flavor, self.rosetta_version)
        LOGGER.info('DDG calculations will be rooted in %s', self._ddg_root)

    @property
    def repo_root_dir(self):
        return self._ddg_config_dict['repo_root']
    @property
    def rosetta_version(self):
        return self._ddg_config_dict['rosetta_version']
    @property
    def rosetta_bin_dir(self):
        return self._ddg_config_dict['rosetta_bin_dir']
    @property
    def rosetta_database_dir(self):
        return self._ddg_config_dict['rosetta_database_dir']

    def _set_structure_filenames(self):
        assert self._structure_dir, "Must compute _structure_dir before calling this function"
        self._cleaned_structure_pdb_filename = os.path.join(self._structure_dir, "%s_%s.pdb" % (
            self._structure_id, self._calculation_flavor))
        self._residue_to_clean_xref_filename = os.path.join(self._structure_dir, "%s_%s_xref.json" % (
            self._structure_id, self._calculation_flavor))
        self._structure_config_filename = os.path.join(self._structure_dir, "%s_%s_structure.config" % (
            self._structure_id, self._calculation_flavor))
        self._slurm_dir = os.path.join(self._structure_dir,"slurm")

    def set_pdb(self, pdb_id: str, chain_id: str) -> str:
        """
        Part 2 of DDG_repo construction.  Adds directory heirarchy from middle 2 pdb_id positions
        :param pdb_id:
        :param chain_id:
        """

        self._structure_source = 'pdb'
        self._structure_id = pdb_id.lower()
        self._chain_id = chain_id
        self._structure_dir = os.path.join(self._ddg_root, 
            'pdb',
            pdb_id[1:3], 
            pdb_id, 
            self._chain_id)
        self._set_structure_filenames()
        return self._structure_dir

    def set_swiss(self, swiss_id: str, chain_id: str) -> str:
        """
        Part 2 of DDG_repo construction.  Adds directory heirarchy from uniprot ID 2-position segments
        """
        self._structure_source = 'swiss'
        self._structure_id = swiss_id
        self._chain_id = chain_id

        uniprot_id = swiss_id[0:6]
        self._structure_dir = os.path.join(
            self._ddg_root,
            'swiss',
            uniprot_id[0:2],
            uniprot_id[2:4],
            uniprot_id[4:6],
            self._structure_id,
            self._chain_id)

        self._set_structure_filenames()

        return self._structure_dir

    def set_modbase(self, modbase_id: str, chain_id: str) -> str:
        """
        Part 2 of DDG_repo construction.  Adds directory heirarchy from uniprot ID 2-position segments
        """
        self._structure_source = 'modbase'
        self._structure_id = modbase_id
        # Often, ancient modbase models lack a chain ID.  So, assign them 'A' as chain
        # in all these cases
        self._chain_id = chain_id if chain_id else 'A'

        # The old modbase IDs tend to be ENSP0000001234567 format
        # which is a clean ENSP protein ID, then followed by
        # bits of punctuation.  To avoid directory jam, we
        # use last 3 digits as directory name
        ensp_id_regex = re.compile('.*(ENSP00+[0-9]+)(.*)')
        ensp_id_match = ensp_id_regex.match(modbase_id)

        assert ensp_id_match,"The modbase_id of %s seems entirely invalid.  Should be ENSP0000123456 and so forth"%modbase_id

        # modbase_underscore_split = modbase_id.split('_')
        # assert len(modbase_underscore_split) < 3
        # modbase_last_6 = modbase_underscore_split[0][:-6:]

        # modbase_last_6 = modbase_id[-6:]
        # self._structure_dir = os.path.join(
            # self._ddg_root,
            # 'modbase',
            # modbase_last_6[0:3],
            # modbase_last_6[3:6],
            # self._structure_id,
            # self._chain_id)

        ENSP_last_3 = ensp_id_match.group(1)[-3:]
        ENSP_suffix = ensp_id_match.group(2)
        self._structure_dir = os.path.join(
            self._ddg_root,
            'modbase',
            ENSP_last_3,
            ENSP_suffix,
            self._structure_id,
            self._chain_id)

        self._set_structure_filenames()

        return self._structure_dir

    def set_usermodel(self, usermodel: str, chain_id: str) -> str:
        """
        Part 2 of DDG_repo construction.  Adds directory heirarchy from uniprot ID 2-position segments
        """
        self._structure_source = 'usermodel'

        self._structure_id = usermodel.split('.')[0].lower()
        self._chain_id = chain_id

        self._structure_dir = os.path.join(
            self._ddg_root,
            'usermodel',
            self._structure_id,
            self._chain_id)

        self._set_structure_filenames()

        return self._structure_dir

    def set_variant(self, variant: Union[str, List[str]]) -> str:
        """
        Part 3 of DDG_repo construction.
        Still needs work to support "all" flavor calculations
        Setup the final calculation directory by specifying
        :param variant: form A123AS (Ref AA/#/(insert code)/Alt AA or "ALL" or other types of strings
        """
        assert self._structure_dir, "Call set_pdb/set_swiss/set_* before calling this function to set _structure_dir"
        self._variant = variant
        # variant of S123AR (insert code A)
        # results in locating to structure_dir/S123A/R/ for all calculations
        self._calculation_dir = os.path.join(self._structure_dir, self._variant[0:-1], self._variant[-1])

        self._psb_status_manager = PsbStatusManager(self._calculation_dir)
        self._log_filename = os.path.join(self._calculation_dir, "%s_%s_%s.log" % (
            self._calculation_flavor, self._structure_id,self._variant))
        return self._calculation_dir

    @property
    def psb_status_manager(self):
        return self._psb_status_manager

    @property
    def chain_id(self):
        """
        Return the chain_id set in __init__
        """
        return self._chain_id

    @property
    def variant(self):
        """
        Return the variant set in __init__
        """
        return self._variant

    @property
    def log_filename(self):
        """The log file that will be used for all calculations"""
        return self._log_filename

    # set_swiss(self,swiss_id:str,chain_id:str,variant:str)

    @staticmethod
    def _residue_id_tuple_from_str(residue_id_as_raw_str):
        """
        De-serialize a Biopython Residue ID from it's str representation
        :param residue_id_as_raw_str:  resid String like "(' ',123,'A')"

        :return: Biopython format residue ID tuple
        """
        assert residue_id_as_raw_str[0] == '('
        assert residue_id_as_raw_str[-1] == ')'
        # Split the ,-separated elements, and strip padding space from around the quote characters
        three_pieces = list(map(str.strip, residue_id_as_raw_str[1:-1].split(',')))
        assert len(three_pieces) == 3
        return (three_pieces[0][1:-1],
                int(three_pieces[1]),
                three_pieces[2][1:-1])

    @property
    def residue_to_clean_xref(self) -> Dict[Tuple, Tuple]:
        """
        Return the cross-reference dictionary which maps original structure residue IDs
        to calculation-ready (typically 1,2,..N) numbered residue ID tuples
        """
        assert self._residue_to_clean_xref_filename, \
            "One of the set_* functions must first be called to setup a specific structure directory "
        if not self._residue_to_clean_xref:
            # If the json file is not in the repo - this falls through and returns None
            if os.path.exists(self._residue_to_clean_xref_filename):
                with open(self._residue_to_clean_xref_filename) as json_file:
                    xref_raw = json.load(json_file)
                    self._residue_to_clean_xref = {}
                    for structure_residue_id_str in xref_raw:  # For each key in a list of key
                        # Convert the value from string to tuple
                        structure_residue_id = DDG_repo._residue_id_tuple_from_str(structure_residue_id_str)
                        rosetta_residue_id = DDG_repo._residue_id_tuple_from_str(xref_raw[structure_residue_id_str])

                        self._residue_to_clean_xref[structure_residue_id] = rosetta_residue_id
                assert self._residue_to_clean_xref,\
                    ("Unable to load residue_to_clean_xref from %s" %
                     self._residue_to_clean_xref_filename)

        return self._residue_to_clean_xref

    @residue_to_clean_xref.setter
    def residue_to_clean_xref(self, residue_to_clean_xref: Dict[Tuple, Tuple]):
        """
        Save a residue_to_clean_xref dictionary to the filesystem
        :param residue_to_clean_xref:
        """
        assert self._residue_to_clean_xref_filename,\
            "One of the set_* functions must first be called to setup a specific structure directory "
        # Only write if not already on disk, and not set earlier
        if ((not self._residue_to_clean_xref) and
                (not os.path.exists(self._residue_to_clean_xref_filename))):
            xref_raw = {}
            for structure_residue_id in residue_to_clean_xref:
                structure_residue_id_str = str(structure_residue_id)
                rosetta_residue_id_str = str(residue_to_clean_xref[structure_residue_id])
                xref_raw[structure_residue_id_str] = rosetta_residue_id_str

            tempfile_name = None
            # Take care to not create an incomplete file
            with tempfile.NamedTemporaryFile(delete=False, mode='w', dir=self._structure_dir) as json_output_temp_file:
                tempfile_name = json_output_temp_file.name
                json.dump(xref_raw, json_output_temp_file)

            if tempfile_name: # great now try to move the xref into position in the repo
                try:
                    os.rename(tempfile_name, self._residue_to_clean_xref_filename)
                except:
                    # IF another process beat us to the punch, then all is well
                    # IF not, then something is quite wrong with our attempt at rename
                    if not os.path.exists(self._residue_to_clean_xref_filename):
                        LOGGER.exception("Unable to rename %s to %s", tempfile_name,
                                         self._residue_to_clean_xref_filename)
                        sys.exit(1)

            self._residue_to_clean_xref = residue_to_clean_xref

    # Possible over-engineering, but we need to store the structure resolution with the cleaned structure and xref dicts
    @property
    def structure_config(self) -> Dict[str, str]:
        """
        Return the config file options used in all calculations on the referenced structure 
        General idea is to compare these to new calculations, to ensure that key parameters
        are not shifting after the repo is created
        """
        assert self._structure_config_filename, \
            "One of the set_* functions must first be called to setup a specific structure directory "
        if not self._structure_config:
            # If the json file is not in the repo - this falls through and returns None
            if os.path.exists(self._structure_config_filename):
                structure_config_parser = configparser.ConfigParser(comment_prefixes='#')
                structure_config_parser.read(self._structure_config_filename)
                self._structure_config = dict(structure_config_parser.items('pdb_info'))
            else:
                return {}
        return self._structure_config

    @structure_config.setter
    def structure_config(self, structure_config: Dict[str, str]):
        """
        Save critical calculation options as entries in human-readable .config file [pdb_info] group
        """
        assert self._structure_config_filename,\
            """One of the set_* functions must first be called to setup a specific structure directory
            """

        structure_config_parser = configparser.RawConfigParser()

        structure_config_parser.add_section('pdb_info')
        for infokey in structure_config:
            structure_config_parser.set('pdb_info', infokey, structure_config[infokey])

        if (not self._residue_to_clean_xref) and (not os.path.exists(self._structure_config_filename)):
            tempfile_name = None
            # Take care to not create an incomplete file
            with tempfile.NamedTemporaryFile(delete=False, mode='w', dir=self._structure_dir) as structure_configfile:
                tempfile_name = structure_configfile.name
                structure_configfile.write(
                    "# Configuration mined from original PDB file %s\n" % datetime.date)
                structure_config_parser.write(structure_configfile)

            if tempfile_name:  # great now try to move the xref into position in the repo
                try:
                    os.rename(tempfile_name, self._structure_config_filename)
                except OSError:
                    # IF another process beat us to the punch, then all is well
                    # IF not, then something is quite wrong with our attempt at rename
                    if not os.path.exists(self._structure_config_filename):
                        LOGGER.exception(
                            "Unable to rename %s to %s",
                            tempfile_name, self._structure_config_filename)
                        sys.exit(1)

        self._structure_config = structure_config

    @property
    def structure_source(self) -> str:
        """Return pdb/swiss/modbase/usermodel as set in set_* function"""
        return self._structure_source
    @property
    def structure_id(self) -> str:
        """Return the structure_id passed to a set_* function"""
        return self._structure_id



    @property
    def cleaned_structure_filename(self) -> str:
        """Return the name of the structure pdb cleaned for Rosetta input"""
        return self._cleaned_structure_pdb_filename

    def mv_cleaned_structure_in(self, cleaned_structure_filename: str) -> str:
        """Import (simply mv) a cleaned structure filename into the repository"""
        assert self._cleaned_structure_pdb_filename
        try:
            os.rename(cleaned_structure_filename, self._cleaned_structure_pdb_filename)
        except OSError:
            # IF another process beat us to the punch, then all is well
            # IF not, then something is quite wrong with our attempt at rename
            if not os.path.exists(self._cleaned_structure_pdb_filename):
                LOGGER.exception("Unable to rename %s to %s", cleaned_structure_filename,
                                 self._cleaned_structure_pdb_filename)
                sys.exit(1)

        return self._cleaned_structure_pdb_filename



    @property
    def calculation_dir(self) -> str:
        """The calculation directory as set by __init__/set*structure/set_variant"""
        return self._calculation_dir

    def make_calculation_directory_heirarchy(self):
        """
        Call os.makedirs(self.calculation_dir) and set permissions for group access
        This will also create the parent variant directory, and grand-parent structure
        directories.
        """
        
        LOGGER.info("Creating repo variant dir drwxrwx---: %s", self.calculation_dir)
        save_umask = os.umask(0)
        os.makedirs(self.calculation_dir, mode=0o770, exist_ok=True)
        os.umask(save_umask)

    @property
    def structure_dir(self) -> str:
        """The structure directory as set by __init__/set*structure"""
        return self._structure_dir

    def slurm_directory_makedirs(self):
        LOGGER.info("Creating slurm directory drwxrwx---: %s", self.slurm_dir)
        save_umask = os.umask(0)
        os.makedirs(self.slurm_dir, mode=0o770, exist_ok=True)
        os.umask(save_umask)
        return self.slurm_dir



    @property
    def slurm_dir(self) -> str:
        """The calculation directory as set by __init__/set*structure/set_variant"""
        return self._slurm_dir

