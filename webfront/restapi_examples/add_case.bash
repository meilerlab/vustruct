#!/usr/bin/bash
curl -X POST http://localhost:3080/add_case \
    -H "Authorization: password" \
    -H 'Content-Type: application/json' \
    -d '{"data_format": "Vanderbilt UDN Case Spreadsheet", 
        "case_uuid": "b60003f4-08bb-4821-83c5-bf6ddaa7e5fa", 
        "case_id": "SAMPLEJUNE2024_3PM"}'
