#!/usr/bin/env bash
# A sample of how to run vustruct_flask from bash shell
export FLASK_APP=./vustruct_flask
export UDN=/dors/capra_lab/users/mothcw/UDNtests
export PATH=/dors/capra_lab/users/mothcw/psbadmin/bin:$PATH
export PYTHONPATH=/dors/capra_lab/users/mothcw/psbadmin/pdbmap:/dors/capra_lab/users/mothcw/psbadmin/bin

export FLASK_ENV=development
cmd='flask run --no-reload -p 3080 --host="0.0.0.0"'
# cmd='flask run --no-reload -p 3080 --host="127.0.0.1"'
echo "Launching: $cmd"
eval $cmd
