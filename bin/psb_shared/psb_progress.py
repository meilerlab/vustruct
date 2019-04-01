#!/usr/bin/env python
def sys_exit_failure(info):
  __write_info_progress(info, "%s: %s\n"%(__file__, inspect.currentframe().f_back.f_lineno))
  # Mark this job as failed
  open("%s/FAILED"%statusdir, 'w').close()
  LOGGER.critical("Terminating in sys_exit_failure(): %s", info)

  sys.exit(info)

def __write_info_progress(info, progress):
  with open('%s/info.new'%statusdir, 'w') as f:
    f.write(info + '\n')
  os.rename('%s/info.new'%statusdir, '%s/info'%statusdir)
  with open('%s/progress.new'%statusdir, 'w') as f:
    f.write(progress)
  os.rename('%s/progress.new'%statusdir, '%s/progress'%statusdir)




