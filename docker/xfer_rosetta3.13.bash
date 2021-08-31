#!/bin/bash
# Script to transfer ONLY needed files from a full Rosetta 3.13 installation to a docker
# image assembly area to accomplish:
# ddg_monomer calculations.  Rosetta 3.13 required
#
# Rosetta 3.13 must be built elsewhere on the same OS environment as will be
# placed in the Docker container, presently Centos 7
# Point ROSETTA_MAIN_SOURCE to the built file area
ROSETTA_MAIN_SOURCE=/TB4/mothcw/rosetta3.13/rosetta_src_2021.16.61629_bundle/main
# POINT ROSETTA_MAIN_DEST to the area where the Docker phase9 image will look for it
ROSETTA_MAIN_DEST=~/psbadmin/rosetta3.13/main

cmd="rm -rfv $ROSETTA_MAIN_DEST"
echo $cmd
eval $cmd

# Create a new area to house the transferred final binary set, libraries, symlinks
mkdir -vp $ROSETTA_MAIN_DEST/source/bin

# Copy over weights, and score function files
declare -a full_directories=(
  'database/citations'
  'database/chemical'
  'database/scoring/weights'
  'database/scoring/score_functions/PairEPotential'
  'database/scoring/score_functions/InterchainPotential'
  'database/scoring/score_functions/EnvPairPotential'
  'database/scoring/score_functions/bondlength_bondangle'
  'database/scoring/score_functions/hbonds'
  'database/scoring/score_functions/disulfides'
  'database/scoring/score_functions/omega'
  'database/scoring/score_functions/rama'
  'database/scoring/score_functions/P_AA_pp'
  'database/input_output'
  'database/sampling/relax_scripts'
  'database/rotamer/shapovalov'
)

cd $ROSETTA_MAIN_DEST
for directory in ${full_directories[@]}
do
    cmd="mkdir --verbose --parents $(dirname $ROSETTA_MAIN_DEST/$directory); cp --preserve --recursive --verbose $ROSETTA_MAIN_SOURCE/$directory/ $ROSETTA_MAIN_DEST/$directory/"
    echo $cmd
    echo
    eval $cmd
done


# Grab only the needed Rosetta binaries for ddg cartesian
mkdir -p $ROSETTA_MAIN_DEST/source/bin

declare -a rosetta_binaries=(
'source/bin/rosetta_scripts.default.linuxgccrelease'
'source/bin/cartesian_ddg.default.linuxgccrelease'
)

for rosetta_binary in ${rosetta_binaries[@]}
do
    cmd="mkdir --verbose --parents $(dirname $rosetta_binary); cp --preserve --recursive --verbose $ROSETTA_MAIN_SOURCE/$rosetta_binary $ROSETTA_MAIN_DEST/$rosetta_binary"
    echo $cmd
    echo
    eval $cmd
    needed_libraries=$(ldd $ROSETTA_MAIN_SOURCE/$rosetta_binary | grep -oE "source/build/.*.so")
    echo "Needed libraries $needed_libraries"
    for needed_library in $needed_libraries
    do
        mkdir --verbose --parents $(dirname $ROSETTA_MAIN_DEST/$needed_library)
        cmd="cp --no-clobber --preserve --verbose --recursive $ROSETTA_MAIN_SOURCE/$needed_library $ROSETTA_MAIN_DEST/$needed_library"
        echo $cmd
        eval $cmd

    done
done
