#!/usr/bin/bash
echo Running Singularity build to create rosetta3.13_ddGCartesian.sif container image file
set -e
cmd='singularity build --fakeroot --docker-login Rosetta3.13_ddGCartesian.sif Rosetta3.13_ddGCartesian.def'
echo "Executing: $cmd"
eval $cmd
