### Steps to create a Rosetta3.13_ddGCartesian.sif singularity container.
1. Install singularity / apptainer

2. Verify your licensing at
  https://github.com/RosettaCommons/rosetta/blob/main/LICENSE.md
  or
  https://els2.comotion.uw.edu/product/rosetta

3. Download the 18GB Rosetta3.13 bundle rosetta_bin_linux_3.13_bundle.tgz 
   which includes all static compiled binaries from https://downloads.rosettacommons.org/downloads/academic/
   
4. Modify the xfer_rosetta_3.13_static.bash script to point to your downloaded tar archive, and a local
   directory to which the needed binaries and databases will be extracted
   Run the script to perform the tar file extraction.
   
6. Modify the %files statement of Rosetta3.13_ddGCartesian.def to point to the local directory to which the needed binaries were extracted above

7. Run build_rosetta3.13_singularlity.bash to create the singularity container Rosetta3.13_ddGCartesian.sif
   (This step requires a dockerhub login to download a Rocky9 base image)
