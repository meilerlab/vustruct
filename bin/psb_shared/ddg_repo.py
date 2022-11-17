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
import pwd
from base64 import b64encode
import grp
import re
from typing import Dict, List, Tuple, Union
import string

import logging
import json
import tempfile
from datetime import datetime
import configparser
import shutil
from pathlib import Path
from psb_shared.psb_progress import PsbStatusManager
LOGGER = logging.getLogger(__name__)
import socket


class DDG_repo():

    @staticmethod
    def ddg_config_read(ddg_config_filename: str) -> Dict[str, str]:
        LOGGER.debug("Attempting to read ddg_config from %s", ddg_config_filename)

        ddg_config = configparser.ConfigParser(allow_no_value=False)

        parsed_file_list = ddg_config.read(ddg_config_filename)
        if not parsed_file_list:
            exit_str = "TERMINATING - Unable to read ddg config file %s" % ddg_config_filename
            LOGGER.critical(exit_str)
            sys.exit(exit_str)

        LOGGER.info("Successfully read ddg_config from " + ddg_config_filename)

        for ddg_config_section in ['ddG', 'ddG_monomer', 'ddG_cartesian', 'os']:
            assert ddg_config.has_section(ddg_config_section),\
                'The ddg_config file %s lacks a %s section'%(ddg_config_filename,ddg_config_section)

        assert 'repo_root' in ddg_config['ddG'],\
                'The ddG section of %s must contain a repo_root entry'%ddg_config_filename

        # Make sure we have basic entries for monomer and cartesian
        for ddg_type in ['ddG_monomer', 'ddG_cartesian']:
            for required_item in [ 'rosetta_version', 'rosetta_bin_dir']:
                if not required_item in ddg_config[ddg_type]:
                    exit_str = "TERMINATING: %s key must be specified in %s section of ddg_config file %s" % (
                            required_item, ddg_type, ddg_config_filename)
                    LOGGER.critical(exit_str)
                    sys.exit(exit_str)


        return ddg_config


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
        else:
            DDG_repo._my_config_dicts[ddg_config_filename] = \
                DDG_repo.ddg_config_read(ddg_config_filename)
            LOGGER.info("ddG Repo config file read: %s"%ddg_config_filename)

        self._ddg_config_dict =  DDG_repo._my_config_dicts[ddg_config_filename][calculation_flavor]
        self._ddg_config_dict['repo_root'] = DDG_repo._my_config_dicts[ddg_config_filename]['ddG']['repo_root']

        self._ddg_config_filename = ddg_config_filename

        # repo/../structure_id/chain/Annn/  parent directory for all calculations on a structure.
        # Contains cleaned PDB and cross reference file to the 1..N rosetta poses
        self._structure_dir = None
        self._variant_dir = None  # _structure_dir/
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

        os_config_dict = DDG_repo._my_config_dicts[ddg_config_filename]['os']
        if 'makedir_mode' in os_config_dict:
            self._os_makedir_mode = int(os_config_dict['makedir_mode'],8)
        else:
            self._os_makedir_mode = 0o770 # Create directories with rwxrwx---

        if 'makedir_group' in os_config_dict:
            self._gr_gid = grp.getgrnam(os_config_dict['makedir_group']).gr_gid
        else:
            self._gr_gid =  os.getgid()
        if 'file_create_mode' in os_config_dict:
            self._file_create_mode = int(os_config_dict['file_create_mode'],8)
        else:
            self._file_create_mode = 0o770 # Create directories with rwxrwx---

        if 'file_create_group' in os_config_dict:
            self._file_create_group = grp.getgrnam(os_config_dict['file_create_group']).gr_gid
        else:
            self._file_create_group =  os.getgid()

        if 'run_umask' in os_config_dict:
            self._run_umask = int(os_config_dict['run_umask'],8)
        else:
            self._run_umask=  0o006


        self._ddg_root = os.path.join(self.repo_root_dir, self._calculation_flavor, self.rosetta_version)
        LOGGER.info('DDG calculations will be rooted in %s', self._ddg_root)

    def __repr__(self):
        return "DDG_repo: %s   repo_root=%s bin_dir=%s"%(
            self._ddg_config_filename,self._ddg_root,self.rosetta_bin_dir
        )

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
    @property
    def rosetta_ld_library_path(self):
        if 'rosetta_ld_library_path' in self._ddg_config_dict:
            return self._ddg_config_dict['rosetta_ld_library_path']
        else:
            return None

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

    def set_alphafold(self, alphafold_id: str, chain_id: str) -> str:
        """
        Part 2 of DDG_repo construction.  Adds directory heirarchy from uniprot ID 2-position segments
        """
        self._structure_source = 'alphafold'
        self._structure_id = alphafold_id
        self._chain_id = chain_id

        # Format is AF-UNPUNP...
        uniprot_id = alphafold_id[3:9]
        self._structure_dir = os.path.join(
            self._ddg_root,
            'alphafold',
            uniprot_id[0:2],
            uniprot_id[2:4],
            uniprot_id[4:6],
            self._structure_id,
            self._chain_id)

        self._set_structure_filenames()

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
        ensp_id_regex = re.compile('[A-Z_]+([0-9]+)(.*)')
        ensp_id_match = ensp_id_regex.match(modbase_id)

        assert ensp_id_match,"The modbase_id of %s seems entirely invalid.  Should be ENSP0000123456 and so forth"%modbase_id

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
        self._variant_dir = os.path.join(self._structure_dir, self._variant[0:-1], self._variant[-1])

        # Set the psb status manager so that it uses our own member functions
        # to perform the os file descriptor creation
        self._psb_status_manager = PsbStatusManager(self._variant_dir,self)
        self._log_filename = os.path.join(self._variant_dir, "%s_%s_%s.log" % (
            self._calculation_flavor, self._structure_id,self._variant))
        return self._variant_dir

    @property
    def psb_status_manager(self):
        return self._psb_status_manager

    @property
    def calculation_flavor(self):
        """

        @return: 'ddg_cartesian' or 'ddg_monomer'  Handy property
        """
        return self._calculation_flavor

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
        # assert self._calculation_flavor != "ddg_cartesian", \
        #    "With ddg_cartesian, there is no cross-referencing of original residue numbers to 1..N"
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
                LOGGER.info("Cross-reference dictionary to cleaned structure parsed from %s",
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

            # Take care to not create an incomplete file
            tempfile_name,tempfile_fd = self.mkstemp(read_or_write='w', dir=self._structure_dir)
            with os.fdopen(tempfile_fd,mode='w') as json_output_temp_file:
                json.dump(xref_raw, json_output_temp_file)

            # self.os_file_chmod(tempfile_name)
            # self.set_group(tempfile_name)

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
            LOGGER.info("Cross-reference dictionary to cleaned structure saved to %s",
                        self._residue_to_clean_xref_filename)

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
                LOGGER.info("Loading structure_config from %s", self._structure_config_filename)
                structure_config_parser = configparser.ConfigParser(comment_prefixes='#')
                structure_config_parser.read(self._structure_config_filename)
                self._structure_config = dict(structure_config_parser.items('pdb_info'))
            else:
                LOGGER.info("No structure config file: %s", self._structure_config_filename)
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
            # Take care to not create an incomplete file

            tempfile_name,tempfile_fd = self.mkstemp(read_or_write='w', dir=self._structure_dir)
            with os.fdopen(tempfile_fd, mode='w') as structure_configfile:
                structure_configfile.write(
                    "# Configuration mined from original PDB file %s\n" % datetime.now)
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
            shutil.move(cleaned_structure_filename, self._cleaned_structure_pdb_filename)
        except OSError:
            # IF another process beat us to the punch, then all is well
            # IF not, then something is quite wrong with our attempt at rename
            if not os.path.exists(self._cleaned_structure_pdb_filename):
                LOGGER.exception("Unable to rename %s to %s", cleaned_structure_filename,
                                 self._cleaned_structure_pdb_filename)
                sys.exit(1)

        self.os_file_chmod(self._cleaned_structure_pdb_filename)
        self.set_group(self._cleaned_structure_pdb_filename)
        return self._cleaned_structure_pdb_filename



    @property
    def variant_dir(self) -> str:
        """The calculation directory as set by __init__/set*structure/set_variant"""
        return self._variant_dir

    def make_variant_directory_heirarchy(self):
        """
        Call os.makedirs(self.variant_dir) and set permissions for group access
        This will also create the parent variant directory, and grand-parent structure
        directories.
        """
        
        LOGGER.info("Creating repo variant dir: %s", self.variant_dir)
        self.makedirs(self.variant_dir, exist_ok=True)

    @property
    def structure_dir(self) -> str:
        """The structure directory as set by __init__/set*structure"""
        return self._structure_dir

    def slurm_directory_makedirs(self):
        LOGGER.info("Creating slurm directory: %s", self.slurm_dir)
        self.makedirs(self.slurm_dir, mode=0o770, exist_ok=True)
        return self.slurm_dir



    @property
    def slurm_dir(self) -> str:
        """The calculation directory as set by __init__/set*structure/set_variant"""
        return self._slurm_dir

    def set_group(self,path : str,logging_active = True) -> None:
        """
        Recursively set the group of all the components of a path
        so that group ownership matches the configuration of the ddg_repo
        """

        save_full_path = path

        # Track chgrp fails for debugging.  Normally these fails are not a big problem
        chown_failures = []
        ddg_root_abspath = os.path.abspath(self._ddg_root)
        chown_attempted = 0
        while path and os.path.abspath(path) != ddg_root_abspath:
            chown_needed = True

            # Often the file group is already fine.  Don't repeat work, warnings, etc if so
            try:
                stat_info = os.stat(path)
                if stat_info.st_gid == self._gr_gid:
                    chown_needed = False
            except FileNotFoundError:
                LOGGER.warning('Attempting set_group() on missing file %s',path)

            if chown_needed:
                chown_attempted += 1
                try:
                    os.chown(path,-1,self._gr_gid)
                except PermissionError: # This is fine, because it just means some other user set the group
                    if len(chown_failures) < 2:
                        chown_failures.append(path)
                    else: # We only track the first failure, and the last, then print them backwards
                        chown_failures[1] = path
                    if logging_active:
                        LOGGER.debug("Failed to os.chown(%s,-1,%d)"%(path,self._gr_gid))

            head, tail = os.path.split(path)
            if tail:
                path= head
            if path and path[-1] == '/':
                path = path[0:-1]

        if chown_attempted > 0 and logging_active:
            if len(chown_failures) == 0:
                LOGGER.info("Group %s=%d set successfully for path %s", \
                    grp.getgrgid(self._gr_gid)[0],self._gr_gid, save_full_path)
            else:
                LOGGER.info("Likely harmless permissions errors setting group %d for %s", \
                    self._gr_gid, str(chown_failures))

    def makedirs(self,name : str, exist_ok = True ) -> None:
        """
        name: fullpath of all directories to be made
        """

        try:
            # Without overriding the prevailing umask, os.makedirs will fail
            old_umask = os.umask(0)

            os.makedirs(name,self._os_makedir_mode,exist_ok)
        finally:
            os.umask(old_umask)

        self.set_group(name)

    def os_open(self, filename: str, read_or_write: str, logging_active=True) -> int:
        """
        Return an integer filedescriptor to a file opened with the permissions
        and group ownership known to ddg_repo through ddg_repo's initialization

        It is imperative to turn off logging _if_ this function is being called
        form inside our custom logger handler 
        """
        assert read_or_write in ['r','a','w']

        file_create_mode =(self._file_create_mode if (read_or_write in ['w','a']) else 0)
        # file_create_mode = 0o660

        old_umask = os.umask(0)

        os_flags = os.O_RDONLY
        if read_or_write == 'a':
            os_flags = os.O_CREAT | os.O_WRONLY | os.O_APPEND
        elif read_or_write == 'w':
            os_flags = os.O_WRONLY | os.O_CREAT | os.O_TRUNC

        fd = os.open(filename,
                     os_flags,
                       file_create_mode, # file_create_mode,
                       )
        if fd < 0:
            message = "Unable to ddg_repo.os.open(%s,%s)"%(filename,read_or_write)
            if logging_active:
                LOGGER.critical(message)
            sys.exit(message)

        if logging_active:
            LOGGER.debug("Success: ddg_repo.os.open(%s,%s,'%s') returning %d",
                    os.path.abspath(filename),
                    oct(file_create_mode),
                    read_or_write,
                    fd)
        self.set_group(filename,logging_active=logging_active)
        os.umask(old_umask)
        return fd

    def touch(self, filename: str):
        old_umask = os.umask(0)
        Path(filename).touch(mode=self._file_create_mode,exist_ok= True)
        os.umask(old_umask)
        self.set_group(filename)

    def os_file_chmod(self,filename: str) -> None:
        old_umask = os.umask(0)
        os.chmod(filename, self._file_create_mode)
        os.umask(old_umask)

    def os_replace(self, source, destination) -> None:
        """
        Call os.replace and add set ddg_repo's group ownership
        """
        os.replace(source,destination)
        self.set_group(destination)

    def os_run_umask_start(self):
        assert not hasattr(self,"_save_mask"),"You did not call os_run_umask_end previously"
        self._save_umask = os.umask(self._run_umask)

    def os_run_umask_end(self):
        assert hasattr(self,"_save_umask"),"You did not call os_run_umask_start"
        os.umask(self._save_umask)
        del self._save_umask

    def _temp_name(self,suffix:str = '', prefix:str = 'tmp_', dir:str = '.') -> str:
        # We just seal the deal with some entropy in the final part of the filename
        base64_encoding_of_8_random_bytes = b64encode(os.urandom(8)).decode('utf-8')
        # BUT - we cannot allow / and other punctuation into the tmp directory name
        random_bytes_without_punctuation = base64_encoding_of_8_random_bytes.translate(
            str.maketrans('','',string.punctuation))

        hostname_without_punctuation = socket.gethostname().translate(
            str.maketrans('','',string.punctuation))


        temp_name = os.path.join(dir,"%s_%s_%s_%d_%s_%s%s"%(
            prefix,
            pwd.getpwuid(os.getuid()).pw_name,
            hostname_without_punctuation,
            os.getpid(),
            datetime.now().strftime('%Y%m%d%H%M%S_%f'), # Timestring down to milliseconds
            random_bytes_without_punctuation,
            suffix
            ))
        return temp_name

    def mkdtemp(self,suffix:str = '', prefix:str = 'tmp_', dir:str = '.') -> str:
        """
        Create a unique directory name composed of 
        os.path.join(dir,prefix+timestamp_milliseconds+randombytes+suffix)
        """

        temp_dirname = self._temp_name(suffix,prefix,dir)
        self.makedirs(temp_dirname,exist_ok = False)
        return temp_dirname

    def mkstemp(self, read_or_write: str, dir='.', logging_active=True) -> (str,int):
        """
        Make a temporary filename, and it along with open int file descriptor
        """
        temp_filename = self._temp_name('.temp','tmpfile',dir)
        return temp_filename,self.os_open(temp_filename,read_or_write,logging_active)


from logging.handlers import RotatingFileHandler
class DDG_repo_RotatingFileHandler(RotatingFileHandler):
    """
    Replace the _builtin_open of the low level FileHandler
    So that log files are opened with proper ddg-repo compatiable
    permissions
    """
    _ddg_repo = None
    def __init__(self,ddg_repo: DDG_repo,backupCount):
        DDG_repo_RotatingFileHandler._ddg_repo = ddg_repo
        super().__init__(ddg_repo.log_filename, backupCount=backupCount)
        self._builtin_open = DDG_repo_RotatingFileHandler.open_func

    def _open(self):
        """
        Open the current base file with the (original) mode and encoding.
        Return the resulting stream.
        """
        return os.fdopen(DDG_repo_RotatingFileHandler._ddg_repo.os_open(self.baseFilename, self.mode, logging_active=False),self.mode, encoding=self.encoding)
        # return DDG_repo_RotatingFileHandler.open_func(self.baseFilename, self.mode, encoding=self.encoding)


    @staticmethod
    def open_func(filename: str, mode: str, encoding, errors='strict'):
        return os.fdopen(DDG_repo_RotatingFileHandler._ddg_repo.os_open(filename,mode, logging_active=False),mode=mode,encoding=encoding)
