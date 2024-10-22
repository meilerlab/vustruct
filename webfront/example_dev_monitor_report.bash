#!/bin/bash
most_recent=`ls -1t $VUS | grep external_user | head -1`
len=${#most_recent}
uuid=${most_recent:$(($len - 36)):36}
echo UUID=$uuid

relaunch=', "relaunch": "True"'
# relaunch=
# unset relaunch
echo relaunch=$relaunch

curl -X POST https://api.vgi01.accre.vanderbilt.edu/dev_monitor_report_loop \
-d '{"case_uuid": "'$uuid'", "case_id": "SUNDAY2PM"'"$relaunch"'}' \
-H 'Authorization: password' -H "Content-Type: application/json"

