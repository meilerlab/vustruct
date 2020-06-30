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
from typing import Dict, List, Tuple, Union

import logging
import json
import tempfile
import datetime
import configparser

LOGGER = logging.getLogger(__name__)


class DDG_repo():
    def __init__(self, repo_root_dir: str, calculation_flavor: str, rosetta_version: str):
        """
        Begin process of constructing a ddg repository manager.
        Call a set_pdb/swiss/etc function after contruction.
        Then, specify a variant to complete final initialization

        :param repo_root_dir:      Ex: /dors/capra_lab/projects/ddg_repo.  
         Typically caller finds in ddg repo's config file
        :param rosetta_version:    Ex: "3.12" Typically caller passes key from config file
        :param calculation_flavor: Ex: "ddg_monomer"  Must be set by caller.
        """

        self._repo_root_dir = repo_root_dir
        self._structure_dir = None  # parent directory for all calculations on a structure.  Contains cleaned PDB and xref
        self._calculation_dir = None  # final directory where all calculucations are conducted
        self._residue_to_clean_xref_filename = None  # Json file to Map original structure residue IDs to the rosetta-ready residue ids
        self._cleaned_structure_pdb_filename = None
        self._structure_config_filename = None
        self._rosetta_version = rosetta_version
        self._calculation_flavor = calculation_flavor
        self._structure_source = None
        self._residue_to_clean_xref = {}
        self._structure_config = None
        self._structure_id = None
        self._chain_id = None
        self._variant = None
        self._log_filename = None

        self._ddg_root = os.path.join(self._repo_root_dir, self._calculation_flavor, self._rosetta_version)
        LOGGER.info('DDG calculations will be rooted in %s', self._ddg_root)

    def _set_structure_filenames(self):
        assert self._structure_dir, "Must compute _structure_dir before calling this function"
        self._cleaned_structure_pdb_filename = os.path.join(self._structure_dir, "%s_%s.pdb" % (
            self._structure_id, self._calculation_flavor))
        self._residue_to_clean_xref_filename = os.path.join(self._structure_dir, "%s_%s_xref.json" % (
            self._structure_id, self._calculation_flavor))
        self._structure_config_filename = os.path.join(self._structure_dir, "%s_%s_structure.config" % (
            self._structure_id, self._calculation_flavor))

    def set_pdb(self, pdb_id: str, chain_id: str) -> str:
        """
        Part 2 of DDG_repo construction.  Adds directory heirarchy from middle 2 pdb_id positions
        :param pdb_id:
        :param chain_id:
        """

        self._structure_id = pdb_id.lower()
        self._chain_id = chain_id
        self._structure_dir = os.path.join(self._ddg_root, pdb_id[1:3], pdb_id, self._chain_id)
        self._set_structure_filenames()
        return self._structure_dir

    def set_swiss(self, swiss_id: str, chain_id: str) -> str:
        """
        Part 2 of DDG_repo construction.  Adds directory heirarchy from uniprot ID 2-position segments
        """
        self._structure_id = swiss_id
        self._chain_id = chain_id

        uniprot_id = swiss_id[0:6]
        self._structure_dir = os.path.join(
            self._repo_root_dir,
            'swiss',
            uniprot_id[0:2],
            uniprot_id[2:4],
            uniprot_id[4:6],
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
        # turns into structure_dir/S123A/R/ for all calculations
        self._calculation_dir = os.path.join(self._structure_dir, self._variant[0:-1], self._variant[-1])
        LOGGER.info("Creating repo variant dir drwxrwx---: %s", self._calculation_dir)
        save_umask = os.umask(0)
        os.makedirs(self._calculation_dir, mode=0o770, exist_ok=True)
        self._log_filename = os.path.join(self._calculation_dir, "%s_%s.log" % (
            self._calculation_flavor, self._structure_id))
        os.umask(save_umask)
        return self._calculation_dir

    @property
    def chain_id(self):
        """
        Return the chain_id set in __init__
        """
        return self._chain_id

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
