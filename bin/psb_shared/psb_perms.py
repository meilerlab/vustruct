#!/usr/bin/env python
import os
import grp
import stat


class PsbPermissions(object):

    def __init__(self, config_dict):
        self.chgrp_to = config_dict['chgrp_to']
        self.gr_gid = grp.getgrnam(self.chgrp_to).gr_gid

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
        os.makedirs(dest_path, exist_ok=True)

        self.set_dir_group_and_sticky_bit(dest_path)
