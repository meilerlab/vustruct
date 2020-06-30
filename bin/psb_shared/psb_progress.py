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

LOGGER = logging.getLogger(__name__)

capra_group = grp.getgrnam('capra_lab').gr_gid

def set_capra_group_sticky(dirname):
    try:
        os.chown(dirname, -1, capra_group)
    except:
        pass

    # Setting the sticky bit on directories also fantastic
    try:
        os.chmod(dirname, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IXGRP | stat.S_ISGID);
    except:
        pass


class PsbStatusManager(object):
    def __init__(self,  status_dir_parent: str) -> None:
        self._status_dir =os.path.join(status_dir_parent,"status")

    def clear_status_dir(self) -> None:
        if not os.path.exists(self._status_dir):
            try:
                os.makedirs(self._status_dir)
            except:
                pass
        set_capra_group_sticky(self._status_dir)

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

    @property
    def complete_file_present(self) -> float:
        complete_filename = os.path.join(self._status_dir, "complete")
        if os.path.exists(complete_filename):
            return os.path.getmtime(complete_filename)
        return False


    def write_info(self, info):
        new_progress_filename = os.path.join(self._status_dir, "progress.new")
        with open(new_progress_filename, 'w') as f:
            f.write("%s: %s\n" % (__file__, inspect.currentframe().f_back.f_lineno))
        os.replace(new_progress_filename,  os.path.join(self._status_dir,"progress"))

        new_info_filename = os.path.join(self._status_dir,"info.new")
        with open(new_info_filename,'w') as f:
            f.write(info + '\n')
        final_info_filename = os.path.join(self._status_dir,"info")
        os.replace(new_info_filename,final_info_filename)
        LOGGER.info("%s now contains: %s"%(final_info_filename,info))

    def _write_info_progress(self,info, progress):
        with open('%s/info.new' % self._status_dir, 'w') as f:
            f.write(info + '\n')
        os.replace('%s/info.new' % self._status_dir, '%s/info' % self._status_dir)
        with open('%s/progress.new' % self._status_dir, 'w') as f:
            f.write(progress)
        os.rename('%s/progress.new' % self._status_dir, '%s/progress' % self._status_dir)

    def sys_exit_failure(self,info_str):
        self._write_info_progress(info_str, "%s: %s\n" % (__file__, inspect.currentframe().f_back.f_lineno))
        # Mark this job as failed
        open("%s/FAILED" % self._status_dir, 'w').close()
        LOGGER.critical("Terminating in sys_exit_failure(): %s", info_str)
        sys.exit(info_str)
