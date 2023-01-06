#!/usr/bin/env bash
# A sample of how to run vustruct_flask from bash shell
export FLASK_APP=./vustruct_flask
export UDN=/dors/capra_lab/users/mothcw/UDNtests
export PYTHONPATH=~/psbadmin/pdbmap:~/psbadmin/bin

# export FLASK_ENV=development
flask run
