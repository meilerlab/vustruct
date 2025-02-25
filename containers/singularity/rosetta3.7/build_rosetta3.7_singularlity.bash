#!/usr/bin/bash
echo Running Singularity build to create rosetta3.7_ddGMonomer.sif container image file
set -e
cmd='singularity build --fakeroot --docker-login Rosetta3.7_ddGMonomer.sif Rosetta3.7_ddGMonomer.def'
echo "Executing: $cmd"
eval $cmd
