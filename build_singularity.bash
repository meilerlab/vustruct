#!/bin/bash 
echo ======================================================
echo Reminder: Did you pugh image_phase9 up to he cloud?
echo ======================================================
rm -f image_phase9.simg
singularity build --docker-login image_phase9.simg Singularity
# had pull before
# docker://chrismoth/image_phase9:latest
