All backend VUStruct pipeline functions are run out of singularity containers.

Broadly speaking, we prefer to directly build singularity containers.  However, in 
some cases, it is still easier to first build a docker image, push it up 
to the vustruct area at hub.docker.com, and then integrate that into a singularity
container build.

For each container, there is a singularity/* subdirectory with the image definition file.

The VUStruct containers include:

### vustruct.simg
The root container in our pipeline, vustruct contains all the python codes in this repo

### ensembl-perlapi
Interconverts ENSEMBL transcript IDs, genomic coordinates, and Amino Acid sequences.
Built from the docker image at: [Dockerhub Repo](https://hub.docker.com/repositories/vustruct)

### ddG3.13
Hosts the Rosetta ddG Cartesian code in a version which gives results closest
to Frenz et al https://pmc.ncbi.nlm.nih.gov/articles/PMC7579412/
You need an (academic) license from Rosetta to build this container

### ddG3.7
Hosts the Rosetta ddG Cartesian code in a version which gives results closest
to Frenz et al https://pmc.ncbi.nlm.nih.gov/articles/PMC7579412/
