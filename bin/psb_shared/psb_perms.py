#!/usr/bin/env python
import os
import grp
import stat
from pathlib import Path

import logging
LOGGER = logging.getLogger(__name__)



class PsbPermissions(object):

    def __init__(self, config_dict):
        self.chgrp_to = config_dict['chgrp_to']
        self.gr_gid = grp.getgrnam(self.chgrp_to).gr_gid

    def os_open(self, filename: str, read_or_write: str) -> int:
        """
        Return an integer filedescriptor to a file opened with the permissions
        and group ownership known to ddg_repo through ddg_repo's initialization
        """
        assert read_or_write == 'r' or read_or_write == 'w'

        # file_create_mode =(self._file_create_mode if (read_or_write == 'w') else 0)
        file_create_mode = 0o660

        # old_umask = os.umask(0)

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
        # os.umask(old_umask)
        return fd

    def touch(self,filename: str):
        Path(filename).touch(exist_ok=True,mode=0o770)


    def set_dir_group_and_sticky_bit(self, dirname):
        """Attempt to mark a directory as being owned by the configured group.
        Then, attempt to mark the directory's sticky bit so that this group
        id propagates to any files created inside this directory"""
        return

        try:
            os.chown(dirname, -1, self.gr_gid)
        except OSError:
            pass

        # Setting the sticky bit on directories also fantastic
        # Simultaneously, make sure users and group have full
        # Read/Write/Execute on the directory
        try:
            os.chmod(dirname,
                     stat.S_IRWXU | stat.S_IRWXG | stat.S_ISGID)
        except OSError:
            pass

    def makedirs(self, dest_path):
        os.makedirs(dest_path, mode=0o777, exist_ok=True)

        self.set_dir_group_and_sticky_bit(dest_path)
