#!/bin/bash 
echo ======================================================
echo Reminder: Did you push chrismoth/diep up to he cloud?
echo ======================================================
rm -f DIEP.simg
singularity build --docker-login DIEP.simg Singularity.DIEP
