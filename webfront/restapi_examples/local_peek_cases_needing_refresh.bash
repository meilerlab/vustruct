#!/usr/bin/bash
#Demonstrate a call from the user-facing webserver, back to vustruct_flask,
#to ask "What runnig cases have updated web reports to publish to website?"
set -ex
curl -X GET http://localhost:3080/peek_cases_needing_refresh -H 'Authorization: ' -H 'Content-Type: application/json'

