### Steps to create a Rosetta3.7_ddGMonomer.sif singularity container.
1. Install singularity / apptainer  (I have had success building the container on the Windows Subsystem for Linux, and Linux itself.  root is not required for this container)

2. Verify your licensing at
  https://github.com/RosettaCommons/rosetta/blob/main/LICENSE.md
  or
  https://els2.comotion.uw.edu/product/rosetta

3. Download the 8.8GB Rosetta3.7 bundle rosetta_bin_linux_3.7_bundle.tgz 
   which includes all static compiled binaries from https://downloads.rosettacommons.org/downloads/academic/
   
4. Modify the xfer_rosetta_3.7_static.bash script to point to your downloaded tar archive, and a local
   directory to which the needed binaries and databases will be extracted
   Run the script to perform the tar file extraction.
   
6. Modify the %files statement of Rosetta3.7_ddGMonomer.def to point to the local directory to which the needed binaries were extracted above

7. Run build_rosetta3.7_singularlity.bash to create the singularity container Rosetta3.7_ddGMonomer.sif
   (This step requires a dockerhub login to download a Rocky9 base image)
