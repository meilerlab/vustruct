# Chimera executable
if __name__ == '__main__':
  from chimera import runCommand as rc
  from chimera import replyobj
  import sys,os,stat
  args = sys.argv
  label = args[1]
  sid   = args[2]
  bio   = args[3]
  pdbf     = "%s_%s_%s.cif"%(label,sid,bio)      # pdb file
  pdbrf    = "%s_%s_%s_renum.cif"%(label,sid,bio)# pdb file (renumbered)
  nvattrf  = "%s_neutral.attr"%label             # neutral variants
  nvrattrf = "%s_renum_neutral.attr"%label       # neutral variants (renumbered)
  pvattrf  = "%s_pathogenic.attr"%label          # pathogenic variants
  pvrattrf = "%s_renum_pathogenic.attr"%label    # pathogenic variants (renumbered)
  avattrf  = "%s_variants.attr"%label            # all variants
  avrattrf = "%s_renum_variants.attr"%label      # all variants (renumbered)
  qattrf   = "%s_quantitative.attr"%label        # quantitative traits
  qrattrf  = "%s_renum_quantitative.attr"%label  # quantitative traits
  ncattrf  = "%s_neutcon.attr"%label             # neutral constraint
  ncrattrf = "%s_renum_neutcon.attr"%label       # neutral constraint (renumbered)
  pcattrf  = "%s_pathcon.attr"%label             # pathogenic constraint
  pcrattrf = "%s_renum_pathcon.attr"%label       # pathogenic constraint (renumbered)
  ppattrf  = "%s_pathprox.attr"%label            # path prox
  pprattrf = "%s_renum_pathprox.attr"%label      # path prox (renumbered)
  qpattrf  = "%s_qtprox.attr"%label              # quantitative path prox
  qprattrf = "%s_renum_qtprox.attr"%label        # quantitative path prox (renumbered)
  out      = "%s"%label                          # output label
  outr     = "%s_renum"%label                    # output label (renumbered)

  print "\nRunning in directory: %s"%os.getcwd()

  def visualize(pdbf,nvattrf,pvattrf,avattrf,qattrf,ncattrf,pcattrf,ppattrf,qpattrf,out,png=True):
    # Visualize with Chimera
    # Assign the annotation to relevant residues
    rc("open %s"%pdbf)
    # Initialize the settings
    rc("background solid white")
    rc("ribbackbone")
    rc("~disp")
    rc("ribspline card smooth strand")
    rc("represent bs")
    rc("setattr m autochain 0")
    rc("setattr m ballScale .7")
    rc("color grey,r")
    ##UPDATE: Removed until Chimera improves "orient" / "fit" functionality
    # Align to neutral variants so that ClinVar and COSMIC match
    # If neutral variants unavailable, use the default orientation
    # if os.path.exists(nvattrf):
      # rc("defattr %s raiseTool false"%nvattrf)
      # rc("define plane name p1 :/neutral")
    # rc("align p1")
    # rc("~define")
    # rc("center")
    # rc("window")
    # Display structure only (include VUS if any specified)
    if png:
      rc("defattr %s raiseTool false"%avattrf)
      rc("disp :/pathogenicity<0 & @ca")
      rc("color black,a :/pathogenicity<0 & @ca")
      rc("transparency 70,r")
      rc("copy file %s_structure.png width 5 height 4 units inches dpi 300"%out)
    # Reduce ribbon transparency for variant plots
    rc("transparency 70,r")
    # Display both pathogenic and neutral
    if os.path.exists(avattrf) and png:
      rc("defattr %s raiseTool false"%avattrf)
      rc("disp :/pathogenicity>-2 & @ca")
      rc("rangecolor pathogenicity,a 1 red 0 blue")
      rc("disp :/pathogenicity<0")
      rc("color black,a :/pathogenicity<0")
      rc("copy file %s_variants.png width 3 height 3 units inches dpi 300"%out)
      rc("~disp")
    # Display neutral only
    if os.path.exists(nvattrf) and png:
      rc("defattr %s raiseTool false"%nvattrf)
      rc("disp :/neutral & @ca")
      rc("color blue,a")
      rc("disp :/pathogenicity<0")
      rc("color black,a :/pathogenicity<0")
      rc("copy file %s_neutral.png width 3 height 3 units inches dpi 300"%out)
      rc("~disp")
    # Display pathogenic only
    if os.path.exists(pvattrf) and png:
      rc("defattr %s raiseTool false"%pvattrf)
      rc("disp :/pathogenic & @ca")
      rc("color red,a")
      rc("disp :/pathogenicity<0")
      rc("color black,a :/pathogenicity<0")
      rc("copy file %s_pathogenic.png width 3 height 3 units inches dpi 300"%out)
      rc("~disp")
    rc("transparency 0,r") # confirm opacity
    # Display quantitative trait
    if os.path.exists(qattrf) and png:
      rc("defattr %s raiseTool false"%qattrf)
      rc("rangecolor quantitative,a min white max red")
      rc("disp :/quantitative & @ca")
      rc("disp :/pathogenicity<0 & @ca")
      rc("copy file %s_quantitative.png width 3 height 3 units inches dpi 300"%out)
    # Display neutral constraint
    if os.path.exists(ncattrf) and png:
      rc("defattr %s raiseTool false"%ncattrf)
      rc("rangecolor neutcon,r min red max white")
      rc("disp :/pathogenicity<0")
      rc("color black,a :/pathogenicity<0")
      rc("copy file %s_neutcon.png width 3 height 3 units inches dpi 300"%out)
    # Display pathogenic constraint
    if os.path.exists(pcattrf) and png:
      rc("defattr %s raiseTool false"%pcattrf)
      rc("rangecolor pathcon,r max red min white")
      rc("disp :/pathogenicity<0")
      rc("color black,a :/pathogenicity<0")
      rc("copy file %s_pathcon.png width 3 height 3 units inches dpi 300"%out)
    # Display PathProx
    if os.path.exists(ppattrf) and png:
      rc("defattr %s raiseTool false"%ppattrf)
      rc("rangecolor pathprox,r max red min blue 0 white")
      rc("disp :/pathogenicity<0")
      rc("color black,a :/pathogenicity<0")
      rc("copy file %s_pathprox.png width 3 height 3 units inches dpi 300"%out)
    # Display qtProx
    if os.path.exists(qpattrf) and png:
      rc("defattr %s raiseTool false"%qpattrf)
      rc("rangecolor qtprox,r max red min white")
      rc("rangecolor qtprox,a max red min white")
      rc("disp :/pathogenicity<0")
      rc("copy file %s_qtprox.png width 3 height 3 units inches dpi 300"%out)
    # Display all variants
    rc("disp :/pathogenicity>-9 & @ca")
    # Do not alter atom colors if displaying quantitative trait
    if not os.path.exists(qattrf):
      rc("rangecolor pathogenicity,a 1 red 0 blue")
    rc("disp :/pathogenicity<0")
    rc("color black,a :/pathogenicity<0")
    rc("save ./%s.py"%out)
    rc("close all")
    # Fix chimera permissions
    os.system("chmod +rx ./%s.py"%out)

  # Visualize with original PDB numbering
  print "\nChimera visualization using original PDB numbering..."
  visualize(pdbf,nvattrf,pvattrf,avattrf,qattrf,ncattrf,pcattrf,ppattrf,qpattrf,out)

  # Visualize with renumbered PDB
  print "\nChimera visualization using renumbered PDB..."
  visualize(pdbrf,nvrattrf,pvrattrf,avrattrf,qrattrf,ncrattrf,pcrattrf,pprattrf,qprattrf,outr,png=False)

  rc("stop now")

