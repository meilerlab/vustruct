#!/bin/bash
cd bin
echo currently in `pwd`
python -m parse_udn_report.py
cd ../pdbmap
echo currently in `pwd`
python -m pdbmap
if [ $? -ne 0 ]; then
  echo `basename "$0"` halting due to errors
  exit 1 
fi

