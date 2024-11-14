All backend VUStruct pipeline functions are run out of singularity containers, which are built from intermedia Docker containers which we push to dockerhub.
For some functions, containers call containers.

For each container, there is a separate subdirectory here which contains a Docker

[Dockerhub Repo](https://hub.docker.com/repositories/vustruct)

TH VUStruct containers are:

### vustruct.simg
The root container in our pipeline, vustruct contains all the python codes in this repo

### ensembl-perlapi
Interconverts ENSEMBL transcript IDs, genomic coordinates, and Amino Acid sequences.

### 


| Container     | Purpose          | Cool  |
| ------------- | ------------- | -----:|
| vustruct.simg | Python Code Base | $1600 |
| col 2 is      | centered         |   $12 |
| zebra stripes | are neat         |    $1 |
