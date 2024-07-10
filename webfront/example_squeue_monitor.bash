#!/bin/bash
most_recent=`ls -1t $VUS | grep external_user | head -1`
len=${#most_recent}
uuid=${most_recent:$(($len - 36)):36}
echo UUID=$uuid

curl -X POST https://api.vgi01.accre.vanderbilt.edu/squeue_monitor \
-d '{"case_uuid": "'$uuid'"}' \
-H 'Authorization: password' -H "Content-Type: application/json"

