#!/bin/bash 
echo ======================================================
echo Reminder: Did you push musite_deep up to he cloud?
echo ======================================================
rm -f MusiteDeep.simg
singularity build --docker-login musitedeep.simg Singularity
# had pull before
# docker://chrismoth/image_phase9:latest
