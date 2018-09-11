# Chimera executable
# This simple script loads a .pdb into chimera and then writes it back
# In doing so, chimera will add HELIX and SHEET notations to the .pdb
# which is great for visualization with ngl
if __name__ == '__main__':
  from chimera import runCommand as rc
  from chimera import replyobj
  import sys,os,stat
  args = sys.argv

  rc("open %s"%args[1])
  rc("write format pdb 0 %s"%args[2])

  rc("stop now")
