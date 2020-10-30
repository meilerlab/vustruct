import grp
import os

## Function Definitions ##
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

def makedirs_capra_lab(DEST_PATH, module_name):
    if not os.path.exists(DEST_PATH):
        os.makedirs(DEST_PATH)
    try:
        assert os.path.exists(DEST_PATH)
    except:
       sys_exit_failure("Fatal: Function %s failed to create destination path %s" % (module_name,DEST_PATH))
    set_capra_group_sticky(DEST_PATH)
