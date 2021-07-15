#!/bin/bash
rm -f image_phase9.simg
docker build --rm -f Dockerfile.phase0 -t chrismoth/image_phase0 . 2>&1 | tee image_phase0.stderr_out  && \
docker build --rm -f Dockerfile.phase1 -t chrismoth/image_phase1 . 2>&1 | tee image_phase1.stderr_out  && \
docker build --rm -f Dockerfile.phase2 -t chrismoth/image_phase2 . 2>&1 | tee image_phase2.stderr_out  && \
docker build --rm -f Dockerfile.phase3 -t chrismoth/image_phase3 . 2>&1 | tee image_phase3.stderr_out  && \
docker build --rm -f Dockerfile.phase9 -t chrismoth/image_phase9 . 2>&1 | tee image_phase9.stderr_out 
