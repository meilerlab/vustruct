#!/bin/bash
# Script to transfer ONLY needed files from a full Rosetta 3.7 installation to a docker
# image assembly area to accomplish:
# ddg_monomer calculations.  Rosetta 3.7 required
#
# Rosetta 3.7 must be built elsewhere on the same OS environment as will be
# placed in the Docker container, presently Centos 7
# Point ROSETTA_MAIN_SOURCE to the built file area
ROSETTA_MAIN_SOURCE=/TB4/mothcw/rosetta3.7/rosetta_src_2016.32.58837_bundle/main
# POINT ROSETTA_MAIN_DEST to the area where the Docker phase9 image will look for it
ROSETTA_MAIN_DEST=~/psbadmin/rosetta3.7/main

cmd="rm -rfv $ROSETTA_MAIN_DEST"
echo $cmd
eval $cmd

# Create a new area to house the transferred final binary set, libraries, symlinks
mkdir -vp $ROSETTA_MAIN_DEST/source/bin

# Copy over needed weights  and scoring files
mkdir -vp $ROSETTA_MAIN_DEST/database/scoring/weights
declare -a needed_weights=(
  'talaris2014.wts'
  'soft_rep_design.wts'
)

for weight_file in ${needed_weights[@]}
do
   cmd="cp --verbose --preserve $ROSETTA_MAIN_SOURCE/database/scoring/weights/$weight_file $ROSETTA_MAIN_DEST/database/scoring/weights"
   echo $cmd
   eval $cmd
done

declare -a full_directories=(
  'database/chemical'
  'database/scoring/score_functions'
  'database/input_output'
  'database/rotamer/ExtendedOpt1-5'
)

cd $ROSETTA_MAIN_DEST
for directory in ${full_directories[@]}
do
    cmd="mkdir --verbose --parents $(dirname $ROSETTA_MAIN_DEST/$directory); cp --preserve --recursive --verbose $ROSETTA_MAIN_SOURCE/$directory/ $ROSETTA_MAIN_DEST/$directory/"
    echo $cmd
    echo
    eval $cmd
done


# Grab only the needed Rosetta binaries for ddg monomer
mkdir -p $ROSETTA_MAIN_DEST/source/build/src/release/linux/3.10/64/x86/gcc/4.8/default

declare -a rosetta_binaries=(
'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/minimize_with_cst.default.linuxgccrelease'
'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/ddg_monomer.default.linuxgccrelease'
'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/per_residue_energies.default.linuxgccrelease'
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
