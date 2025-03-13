#!/bin/bash
# Taken directly from http://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html#singularity
singularity pull --name vep.sif docker://ensemblorg/ensembl-vep
# You must then run the container to install compatible cache files from same vesion using a command like
# singularity exec vep.sif INSTALL.pl -c $HOME/somewhere/vep_data -a cf -s homo_sapiens -y GRCh38
