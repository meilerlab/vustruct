# VUStruct Pipeline Installation Guide for a SLURM HPC

## Introduction

This guide provides step-by-step instructions for installing and configuring the VUStruct pipeline on a High-Performance
Computing (HPC) cluster running SLURM workload manager. The VUStruct pipeline is designed for variant annotation and
structural analysis.

**Prerequisites:**

- Access to a SLURM-based HPC cluster
- User account with appropriate permissions
- A Mariadb or MYSQL database that can be accessed by pipeline code\
running on any compute node.
- Basic knowledge of the Linux command line
- Sufficient storage space for reference genomes and databases

## HPC Filesystem Pre-planning

Before proceeding to installation, it is helpful to sketch out a basic filesystem plan.  
In a production environment there are often constraints on available storage, limits,
or data retention time. Wide variations in storage cost and I/O performance -
can be dependent on where data files are located.

For clarity and generality of documentation, this document greedily locates key sets of VUStruct
files directly off the user $HOME directory.  

| directory for documentation | function                              |
|:----------------------------|:--------------------------------------|
| ~/vustruct                  | VUStruct software                     |
| ~/VUS                       | $VUS, Each case run is a subdirectory |
| ~/data                      | downloaded data files                 |
| ~/containers                | apptainer images                      |
| ~/ddg_repo                  | Stored prior ddG calculation results  |

In this case, we accept the thought experiment that all storage
below $HOME is unlimited,
and that the storage is performant and persistent without expiration.

#### EMPHASIS:
**In practice, you will most likely use a variety 
of paths to house these functions**

( We have not yet discussed the SQL database.  Those files be local to the 
MariaDB instance.  You should have 2TB of disk storage, but the filesystem
under the SQL layer isnot seen by any pipeline code of course. )

## Loading VUStruct Software

Populating your ~/vustruct directory requires two git clones.

First you clone this github repo to create ~/vustruct:

```
cd ~
git clone https://github.com/meilerlab/vustruct.git

cd vustruct
git status 
```

Second, clone the supporting pdbmap supoprting libraru beneath vustruct.  
Continuing from the vustruct directory above:

```
cd ~/vustruct   # You should already be here from above.
git clone https://github.com/CapraLab/pdbmap/
cd pdbmap
git status
```

Developer note: Although pdbmap resides "under" vustruct, it does not 
"pollute" the vustruct git repo due to a "pdbmap/" entry in the 
~/vustruct/.gitignore file.

## Downloading ~/data with the download scripts
The VUStruct pipeline makes extensive use of pre-downloaded external data.
As a general design constraint, we demand that VUStruct can run patient
cases when external sources are off-line.  

For each type of data, there is a mature .bash script 
with documentation in the top lines, which you should quickly read
for warnings and tips, before launching the script.

| data source | function                                                  | SQL integration required? |
|:------------|:----------------------------------------------------------|:-------------------------|
| alphafold   | Alphafold 2 structures                                    |                          |
| clinvar     | Clinvar (includes pathogenic) variants for pathprox calcs | X                        |
| ensembl     | vep cache file and huge GRCh38 SQL files                  | no and X                 |
| gnomad      | Gnomad variants (population diversity) for pathprox       | X                        |                        
| interpro    | Interpro protein domains                                  |                          |
| modbase     | Modbase protein homology models                           |                          |
| pdb         | Structured mirror of all experimental depositions         |                          |
| sifts       | Alignments from Uniprot to PDB chains                     | X                        |
| swissmodel  | Structured mirror of all experimental depositions         | |                         |
| uniprot     | Extracts the critical IDmapping file.  | X                        |


When downloading an individual data type, it is essential that _all_ the data is downloaded,
that .yaml version files are created, and that occasional complex  post processing processes are copiously logged.

In some cases, after download, data must be integrated in 
the pipeline's SQL database.

To get started with the download scripts for first time, do

```
cd ~
git clone https://github.com/CapraLab/pipeline_download_scripts.git
```
A challenge with the download scripts is that the data sources are continually changing their
approaches for publishing data.   Please get in touch if you encounter problems.



** STOP !! **  I will resume working on this documentation Thursday

## Integrating downloaded data into the SQL database
Over recent years, I have strived to bring the pipeline closer to the supplied source data,
and less dependent on pre-processing and saving of processed data to SQL.  For example,
the gnomad files are now read directly by the pathoprox program using tabix.  Previously,
lookups of population frequencies in gnomad files required SQL queries.

The 

### Current 

#### Ensembl Genome
-- Describe here how to use the public servers as alternative.

#### Clinvar cross-reference to genome

## Download of the Ensembl Variant Effect Predictor (VEP) and supporting files
Note that this is run frmo a container



With the exception of Your plan will encompass the varying types of storage that are often available on an HPC cluster. Consider the following:

- **Reference Genomes:** Store reference genomes in a dedicated directory to ensure easy access and management.
- **Databases:** Organize databases in a separate directory to maintain a clean and organized filesystem.
- **Pipeline Data:** Create a directory for pipeline-specific data, such as intermediate results and final outputs.

### VUStruct backend software
If convenient, 
It can be convenient to clone the 


Consider the following:

- **Reference Genomes:** Store reference genomes in a dedicated directory to ensure easy access and management.
- **Databases:** Organize databases in a separate directory to maintain a clean and organized filesystem.
- **Pipeline Data:** Create a directory for pipeline-specific data, such as intermediate results and final outputs.

## HPC Environment Preparation via .bashrc
We have found the default scientific software stacks to be excellent points
of departure for installation of the customized python 3.11 stack.  On Vanderbilt's
ACCRE cluster, make the following changes to the default .bashrc:

1.  So that commands to VUStruct can be sent via sshkeys, it is\
imperative that the .bashrc not output text when a shell is non-interactive.
You can test this with "ssh -q your_id@yourcluster_login_node"
```
if [[ $- = *i* ]]
then
# then it is an interactive login and you can output text freely
# for the human user to read
fi
```
2. Setup your scientific software stack.  This will vary wildly from
cluster to cluster.  We use at ACCRE:

```
# Load ACCRE software stack
setup_accre_software_stack
export FLEXIBLAS=OPENBLAS
module load scipy-stack/2024b
```



```bash
module load python/3.11
```

sufficient for our needs.\
However, if you need to customize your environment, you can modify your `.bashrc` file to include any necessary module loads or environment variables.

### Python Environment Setup

1. Python 3.11 is used by VUSTruct to run the interactive command line \
interface programs, the pathprox3.py code, and the container launchers.\
While many imaginable routes to  exist for You will create a python3 virtual environment using `conda` or `venv`:
2. 


Our goal remains containerization of this "outside the container code" but the specter of 


   python3 -m venv vustruct
   source vustruct/bin/activate

### .bashrc
Paths into the above.