import grp
import os

## Function Definitions ##
# DO NOT CHECK THIS IN BELOW IT IS CRAP
try:
    capra_group = grp.getgrnam('capra_lab').gr_gid
except:
    capra_group = os.getegid()
# DO NOT CHECK THIS IN ABOVE IT IS CRAP

def set_capra_group_sticky(dirname):
  return
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
