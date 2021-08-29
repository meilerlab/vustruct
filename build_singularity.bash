#!/bin/bash 
rm -f image_phase9.simg
singularity build --docker-login image_phase9.simg Singularity
# had pull before
# docker://chrismoth/image_phase9:latest
