#!/usr/bin/env bash
# A sample of how to run vustruct_flask from bash shell
export FLASK_APP=./vustruct_webupdate
export UDN=/dors/capra_lab/users/mothcw/VUStruct
export PYTHONPATH=~/psbadmin/pdbmap:~/psbadmin/bin

export FLASK_ENV=development
cmd='flask run --no-reload -p 3000'
echo "Launching: $cmd"
eval $cmd
