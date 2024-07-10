#!/usr/bin/env python
"""
Define PSB_Status class which is used by ddg - and hopefully later integrated into pathprox3
to process the outputdirectory/self._status_dir files that record updates that monitor applications
can see
"""
import grp
import logging
import os
import stat
import sys
import inspect
import logging
from pathlib import Path
from typing import Tuple

LOGGER = logging.getLogger(__name__)

# Provide a simple file opener in case it is not overridden
class _DefaultOsInterface(object):
    def os_open(self, filename: str, read_or_write: str) -> int:
        """
        Return an integer filedescriptor to a file opened with the permissions
        and group ownership known to ddg_repo through ddg_repo's initialization
        """
        assert read_or_write == 'r' or read_or_write == 'w'

        # file_create_mode =(self._file_create_mode if (read_or_write == 'w') else 0)
        file_create_mode = 0o660

        old_umask = os.umask(0)

        fd = os.open(filename,
                     (os.O_CREAT | os.O_WRONLY | os.O_TRUNC) if (read_or_write == 'w') else os.O_RDONLY,
                       )

        if fd < 0:
            message = "Unable to os.open(%s,%s)"%(filename,read_or_write)
            LOGGER.critical(message)
            sys.exit(message)
        LOGGER.info("Success: os.open(%s,%s,'%s') returning %d",
                    os.path.abspath(filename),
                    oct(file_create_mode),
                    read_or_write,
                    fd)
        # self.set_group(filename)
        os.umask(old_umask)
        return fd

    @property
    def _os_makedir_mode(self) -> int:
        return 0o660


    def makedirs(self,name : str, exist_ok = True ) -> None:
        """
        name: fullpath of all directories to be made
        """
        os.makedirs(name,mode=self._os_makedir_mode,exist_ok=exist_ok)

    def touch(self,filename: str):
        Path(filename).touch(exist_ok=True,mode=0o770)


class PsbStatusManager(object):
    def __init__(self,  status_dir_parent: str, os_interface = _DefaultOsInterface()) -> None:
        self._status_dir = os.path.abspath(os.path.join(status_dir_parent,"status"))
        self._complete_filename = os.path.join(self._status_dir, "complete")
        self._failed_filename = os.path.join(self._status_dir, "FAILED")
        self._progress_filename = os.path.join(self._status_dir, "progress")
        self._info_filename = os.path.join(self._status_dir, "info")
        self._os_interface = os_interface

    def clear_status_dir(self) -> None:
        if not os.path.exists(self._status_dir):
            # import pdb; pdb.set_trace()
            try:
                self._os_interface.makedirs(self._status_dir)
            except:
                pass
        # set_capra_group_sticky(self._status_dir)

        for the_file in os.listdir(self._status_dir):
            file_path = os.path.join(self._status_dir, the_file)
            LOGGER.info('Deleting old file %s' % file_path)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
            except Exception as e:
                LOGGER.exception("Unable to delete file in status directory")
                sys.exit(1)

    @property
    def status_dir(self):
        return self._status_dir


    def mark_complete(self) -> None:
        LOGGER.info("Creating %s file to mark process successful end"%self._complete_filename)
        self._os_interface.touch(self._complete_filename)

    @property
    def complete_file_present(self) -> float:
        if os.path.exists(self._complete_filename):
            return os.path.getmtime(self._complete_filename)
        return False

    def mark_failed(self) -> None:
        LOGGER.info("Creating %s file to mark process final failure"%self._failed_filename)
        self._os_interface.touch(self._failed_filename)

    @property
    def failed_file_present(self) -> float:
        if os.path.exists(self._failed_filename):
            return os.path.getmtime(self._failed_filename)
        return False


    def read_info_progress(self) -> Tuple[str,str]:
        """
        Read any status info and return any info that has been written as the program progresses
        
        :return Tuple(info,progress) both str
        """
        info_str = None
        progress_str = None

        if os.path.exists(self._info_filename):
            try:
                with open(self._info_filename,'r') as fp:
                    info_str = fp.read()
            except OSError as err:
                pass

        if os.path.exists(self._progress_filename):
            try:
                with open(self._progress_filename,'r') as fp:
                    progress_str = fp.read()
            except OSError as err:
                pass

        return (info_str,progress_str)

    def write_info(self, info):
        new_progress_filename = os.path.join(self._status_dir, "progress.new")
        with os.fdopen(self._os_interface.os_open(new_progress_filename, 'w'),'w') as f:
            f.write("%s: %s\n" % (__file__, inspect.currentframe().f_back.f_lineno))
        os.replace(new_progress_filename,  self._progress_filename)

        new_info_filename = os.path.join(self._status_dir,"info.new")
        with os.fdopen(self._os_interface.os_open(new_info_filename,'w'),'w') as f:
            f.write(info + '\n')
        os.replace(new_info_filename,self._info_filename)
        LOGGER.info("%s now contains: %s"%(self._info_filename,info))

    def _write_info_progress(self,info, progress):
        with os.fdopen(self._os_interface.os_open('%s/info.new' % self._status_dir,'w'), 'w') as f:
            f.write(info + '\n')
        os.replace('%s/info.new' % self._status_dir, '%s/info' % self._status_dir)
        with os.fdopen(self._os_interface.os_open('%s/progress.new' % self._status_dir, 'w'),'w') as f:
            f.write(progress)
        os.rename('%s/progress.new' % self._status_dir, '%s/progress' % self._status_dir)

    def sys_exit_failure(self,info_str):
        self._write_info_progress(info_str, "%s: %s\n" % (__file__, inspect.currentframe().f_back.f_lineno))
        # Mark this job as failed
        self.mark_failed()
        LOGGER.critical("Terminating in sys_exit_failure(): %s", info_str)
        sys.exit(info_str)
