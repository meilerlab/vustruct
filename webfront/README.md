# Web front end to vustruct
This directory contains two python scripts which
implement the bulk of a web front-end to the
vustruct pipeline 

## <ins>vustruct_flask</ins> 

Via flask, implements RESTAPI calls which 
1. return unique job identifiers and
2. launch the needed psb_*.py command line programs on the compute cluster
3. record and report  emerging website availability
4. inform vustruct_webupdate of running job IDs which are
needing a web refresh.

Because vustruct_flask will launch prep tasks and
cluster tasks, it is executed outside the container.

## <ins>vustruct_flask</ins> 
