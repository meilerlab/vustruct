#!/bin/bash
# Script to create "final" pipeline docker image file chrismoth/image_phase9 
# Must be run as root
# Requires that Rosetaa 3.13 and 3.7 files have been copied in to the bulid area
# Using hte xfer_rosetta*.bash sciprts
# Also requires that chimera_headless has been copied in
rm -f image_phase9.simg
docker build --rm -f docker/Dockerfile.0.centos_anaconda -t chrismoth/image_phase0 . 2>&1 | tee image_phase0.stderr_out  && \
docker build --rm -f docker/Dockerfile.1.perl            -t chrismoth/image_phase1 . 2>&1 | tee image_phase1.stderr_out  && \
docker build --rm -f docker/Dockerfile.2.ensembl_bioperl -t chrismoth/image_phase2 . 2>&1 | tee image_phase2.stderr_out  && \
docker build --rm -f docker/Dockerfile.3.chimera         -t chrismoth/image_phase3 . 2>&1 | tee image_phase3.stderr_out  && \
docker build --rm -f docker/Dockerfile.9.pip_rosetta     -t chrismoth/image_phase9 . 2>&1 | tee image_phase9.stderr_out 
