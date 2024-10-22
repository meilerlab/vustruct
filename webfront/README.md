# Web front end to vustruct
This directory contains two python scripts which
implement the bulk of a web front-end to the
vustruct pipeline.  I also check in the tail
of the wordpress functions.php which POSTS
to the vustruct_* flask routines

## <ins>vustruct_flask</ins> 

A flask program which 
- runs on an HPC compute node and
- implements RESTAPI calls which:
1. return unique job identifiers to wordpress .php code and
2. launch the needed psb_*.py command line programs on the compute cluster
3. record and report  emerging website availability
4. inform vustruct_webupdate of running job IDs which are
needing a web refresh.

Because vustruct_flask will launch prep tasks and
cluster tasks, it is executed outside the container.

## <ins>vustruct_webupdate</ins>

A second flask program which 
- runs local alongside the wordpress server and
- implements RESTAPI calls which:
1. record job identifiers which are in process
2. rebuild local websites in file system from emerging vustruct .tar.gz files

## <ins>functions.php.tail</ins>

Hook forms to validate gravity forms user interface entries
and launch the pipeline as final part of SUBMIT processing

