#!/usr/bin/env bash
# Let the world see what we are doing
set -x
# A sample of how to run vustruct_flask from bash shell
export FLASK_APP=./vustruct_flask
# export FLASK_APP="vustruct_flask:create_app('--config=$UDN/config/global.config --userconfig=$UDN/config/mothcw.config')"

export PATH=/dors/capra_lab/users/mothcw/psbadmin/bin:$PATH
export PYTHONPATH=/dors/capra_lab/users/mothcw/psbadmin/pdbmap:/dors/capra_lab/users/mothcw/psbadmin/bin
export FLASK_RUN_PORT=3080
export FLASK_RUN_HOST="0.0.0.0"

# deprecated
export FLASK_ENV=development
export FLASK_DEBUG=True

# Arguments read by the flast program
UDN=/dors/capra_lab/users/mothcw/VUStruct
export VUSTRUCT_CONFIG=$UDN/config/global.config
export VUSTRUCT_USERCONFIG=$UDN/config/mothcw.config

flask --debug run
exit 0

while [ 1 ]
do
if netstat -ltnp | grep -q $FLASK_RUN_PORT; then
    echo "flask running on port $FLASK_RUN_PORT"
    echo "killall flask"
    killall flask
    sleep 5
fi

echo restarting vustruct_flask | mail -s 'Flask restart' chris.moth@vanderbilt.edu
flask --debug run 
echo Halted vustruct_flask | mail -s 'Urgen Flask halt' chris.moth@vanderbilt.edu
sleep 60
done

