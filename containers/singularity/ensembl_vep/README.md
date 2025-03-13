# Building the ENSEMBL Variant Effect Predictor (vep.sif)
The apptainer container is built from the .bash file in this directory, which simply copies the advice
given from ENSEMBL at this URL:
https://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html#singularity

Note that matching vepcache data also must be downloaded separately from ENSEMBL per
the above documentation and as documented in teh .bash file.

At run time, the cache data directory must be "bound" (apptainer provides numerous options)
to run.

The GRCh38 release version must match between the vep.sif buuld and the downloaded vep cache.
