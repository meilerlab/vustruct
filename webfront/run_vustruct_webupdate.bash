#!/usr/bin/env bash
# A sample of how to run vustruct_flask from bash shell
set -x
export FLASK_APP=./vustruct_webupdate
export UDN=/dors/capra_lab/users/mothcw/VUStruct
export PYTHONPATH=/dors/capra_lab/users/mothcw/psbadmin/pdbmap:/dors/capra_lab/users/mothcw/psbadmin/bin

# export FLASK_ENV=development
export FLASK_ENV=development
export FLASK_DEBUG=True
cmd='flask --debug run --no-reload -p 3000'
echo "Launching: $cmd"
eval $cmd
