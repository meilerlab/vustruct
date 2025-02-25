#!/bin/bash
# Script to transfer ONLY needed files from a full Rosetta 3.13 installation
# from the static binary download provided by rosetta commons.
#
# See xfer_rosetta3.13.bash for extraction of binaries and dynamic linked libraries
# that you have built on Rocky linux 8 from the source download from rosetta commons

# Following "set -e", any failure of a commend exits the script
set -e

# ROSETTA_TAR_BUNDLE=/rhenium/ssd1/rosetta_bin_linux_3.13_bundle.tgz
ROSETTA_TAR_BUNDLE=../Downloads/rosetta_bin_linux_3.13_bundle.tgz
# POINT ROSETTA_MAIN_DEST to the area where container will look for it
# During extraction the main/ directory will be expanded under this one
ROSETTA_MAIN_DEST=../rosetta3.13_ddGCartesian

# These below should not need to be changed
ROSETTA_BUNDLE_ROOT=rosetta_bin_linux_2021.16.61629_bundle
ROSETTA_BUNDLE_DATABASES=main/database
ROSETTA_BUNDLE_BINARIES=main/source/build/src/release/linux/3.10/64/x86/gcc/4.8/static
ROSETTA_BUNDLE_SYMLINKS=main/source/bin
ROSETTA_STATIC_BINARY_SUFFIX=.static.linuxgccrelease


echo This script will un-tar needed ddGCartesian files from $ROSETTA_TAR_BUNDLE to $ROSETTA_MAIN_DEST

cmd="rm -rfv $ROSETTA_MAIN_DEST"
echo $cmd
eval $cmd

echo Create a new area to house the transferred final binary set, libraries, symlinks
mkdir --verbose --parents $ROSETTA_MAIN_DEST/source/bin

ddG_filelist_file=/tmp/rosetta.filelist
rm -f $ddG_filelist_file
echo We create a list of needed files in $ROSETTA_FILES
echo Copy over weights, and score function database files
declare -a databases=(
  'citations'
  'chemical'
  'scoring/weights'
  'scoring/score_functions'
  'input_output'
  'sampling'
  'rotamer/shapovalov'
)

cd $ROSETTA_MAIN_DEST
for database in ${databases[@]}
do
    echo $ROSETTA_BUNDLE_ROOT/$ROSETTA_BUNDLE_DATABASES/$database | tee -a $ddG_filelist_file
done

echo Grab only the needed Rosetta binaries for ddg cartesian

declare -a rosetta_binaries=(
'rosetta_scripts'
'cartesian_ddg'
)

for binary in ${rosetta_binaries[@]}
do
    echo $ROSETTA_BUNDLE_ROOT/$ROSETTA_BUNDLE_BINARIES/$binary$ROSETTA_STATIC_BINARY_SUFFIX | tee -a $ddG_filelist_file
    echo $ROSETTA_BUNDLE_ROOT/$ROSETTA_BUNDLE_SYMLINKS/$binary$ROSETTA_STATIC_BINARY_SUFFIX | tee -a $ddG_filelist_file
done

echo Unpacking the needed ddG files from the downloaded rosetta bundle.
echo The files will be unpacked to $ROSETTA_MAIN_DEST
echo The leading $ROSETTA_BUNDLE_ROOT will be stripped
cmd="tar xvhf $ROSETTA_TAR_BUNDLE -C $ROSETTA_MAIN_DEST --files-from $ddG_filelist_file --strip-components=1"
echo $cmd
eval $cmd

echo All script commands completed without error
echo All needed ddGCartesian files have been un-tar\'d from $ROSETTA_TAR_BUNDLE to $ROSETTA_MAIN_DEST
echo "$ROSETTA_MAIN_DEST/source contains $(find $ROSETTA_MAIN_DEST/source -type f | wc -l) files"
echo "$ROSETTA_MAIN_DEST/database contains $(find $ROSETTA_MAIN_DEST/database -type f | wc -l) files"
echo Containers should copy these files to /opt and prefix the PATH and LD_LIBRARY_PATH with:
echo PATH=/opt/rosetta3.13_ddGCartesian/main/source/bin:\$PATH
