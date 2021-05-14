#!/bin/bash
# Script to transfer needed files from a full Rosetta installation to accomplish:
# ddg_monomer calculations.  Rosetta 3.7 required
ROSETTA_MAIN_SOURCE=/TB4/mothcw/rosetta3.7/rosetta_src_2016.32.58837_bundle/main
# ROSETTA_MAIN_SOURCE=/TB4/mothcw/rosetta3.7.was/rosetta_src_2016.32.58837_bundle/main
ROSETTA_MAIN_DEST=./rosetta3.7/main

cmd=rm -rf $ROSETTA_MAIN_DEST
echo $cmd
eval $cmd

# Create a new area to house the transferred final binary set, libraries, symlinks
mkdir -vp $ROSETTA_MAIN_DEST/bin

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
  'database/scoring/score_functions/PairEPotential'
  'database/scoring/score_functions/hbonds'
  'database/scoring/score_functions/disulfides'
  'database/scoring/score_functions/rama'
  'database/scoring/score_functions/P_AA_pp'
  'database/input_output'
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

#declare -a rosetta_internal_libs=(
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libdevel.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libprotocols.7.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libprotocols.6.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libprotocols_e.5.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libprotocols_d.5.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libprotocols_c.5.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libprotocols_b.5.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libprotocols_a.5.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libprotocols_h.4.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libprotocols_g.4.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libprotocols_f.4.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libprotocols_e.4.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libprotocols_d.4.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libprotocols_c.4.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libprotocols_b.4.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libprotocols_a.4.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libprotocols.3.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libprotocols_b.2.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libprotocols_a.2.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libprotocols.1.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libcore.5.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libcore.4.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libcore.3.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libcore.2.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libcore.1.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libbasic.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libnumeric.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libutility.so'
#'source/build/src/release/linux/3.10/64/x86/gcc/4.8/default/libObjexxFCL.so'
#)
#
#declare -a rosetta_external_libs=(
#'source/build/external/release/linux/3.10/64/x86/gcc/4.8/default/libcppdb.so'
#'source/build/external/release/linux/3.10/64/x86/gcc/4.8/default/libsqlite3.so'
#'source/build/external/release/linux/3.10/64/x86/gcc/4.8/default/libcifparse.so'
#'source/build/external/release/linux/3.10/64/x86/gcc/4.8/default/libxml2.so'
#)
#
#ROSETTA_EXTERNAL_LIB_DIR=$ROSETTA_MAIN_DEST/source/build/external/release/linux/3.10/64/x86/gcc/4.8/default
#mkdir -pv $ROSETTA_EXTERNAL_LIB_DIR
#for library in ${rosetta_external_libs[@]}
#do
#    echo $ROSETTA_MAIN_DEST/$directory
#    mkdir --verbose --parents $(dirname $ROSETTA_MAIN_DEST/$needed_library)
#    cmd="cp --no-clobber --preserve --verbose --recursive $ROSETTA_MAIN_SOURCE/$needed_library $ROSETTA_MAIN_DEST/$needed_library"
#    echo $cmd
#    eval $cmd
#     
#done
#
