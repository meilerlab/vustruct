#!/usr/bin/env bash
#
# run_vustruct_flask.bash launches the pipeline process manager which runs on a "gateway node"
# of our slurm cluster.
# The present code has far too many hard coded references to directories
# which are no longer used, especially PERL modules which are running in containers anyway
# 
# Let the world see what we are doing
set -x
set -e
mkdir -p log
echo echo $(date) "Starting run_vustruct.bash " >> log/run_vustruct_flask.log

###################################################################3
# Setup up the environment as if we logged in
###################################################################3
export psbbin=/dors/capra_lab/users/mothcw/psbadmin/bin

PATH=/dors/capra_lab/users/mothcw/psbadmin/bin:/dors/capra_lab/users/mothcw/psbadmin/pdbmap:/dors/capra_lab/users/mothcw/psbadmin/pathprox:$PATH
PATH=/home/mothcw/anaconda3/bin:$PATH
export PYTHONPATH=/dors/capra_lab/users/mothcw/psbadmin/bin:/dors/capra_lab/users/mothcw/psbadmin/pdbmap:$PYTHONPATH

##############################################################################
# POINT TO OUR LOCAL PERL INSTALLATION
#
PERL_LOCAL_ROOT="/dors/capra_lab/users/mothcw/localperl"
PATH=$PERL_LOCAL_ROOT/bin:$PATH
export PERL5LIB=$PERL_LOCAL_ROOT/lib

# TWO Arcane flags which instruct cpanm where to install PERL modules
# Since our PERL built robustly for pipeline, doubt this is super critical
export PERL_MB_OPT="--install_base $PERL_LOCAL_ROOT"
export PERL_MM_OPT="INSTALL_BASE=$PERL_LOCAL_ROOT"
##############################################################################

PERL_ENSEMBL_ROOT=/dors/capra_lab/users/mothcw/ensembl
# newPERL5LIB=$newPERL5LIB:$PSBADMIN_PERL_ROOT/localperl/lib/perl5/x86_64-linux/
newPERL5LIB=$PERL_ENSEMBL_ROOT/ensembl/modules
newPERL5LIB=$newPERL5LIB:$PERL_ENSEMBL_ROOT/ensembl-variation/modules
newPERL5LIB=$newPERL5LIB:$PERL_ENSEMBL_ROOT/ensembl-compara/modules
newPERL5LIB=$newPERL5LIB:$PERL_ENSEMBL_ROOT/ensembl-funcgen/modules
newPERL5LIB=$newPERL5LIB:$PERL_ENSEMBL_ROOT/ensembl-tools/modules
newPERL5LIB=$newPERL5LIB:$PERL_ENSEMBL_ROOT/ensembl-io/modules
newPERL5LIB=$newPERL5LIB:$PERL_ENSEMBL_ROOT/bioperl-live/
export PERL5LIB=$newPERL5LIB:$PERL5LIB

# export ENSEMBL_REGISTRY=/dors/capra_lab/users/mothcw/psbadmin/pdbmap/EnsEMBL/ensembl_registry_GRCh38.conf
export ENSEMBL_REGISTRY=/dors/capra_lab/users/mothcw/VUStruct/config/ensembl_registry_GRCh38.conf
unset OPT

##############################################################################
# The next line updates PATH for the Google Cloud SDK.
# if [ -f '/dors/capra_lab/data/gnomad/google-cloud-sdk/path.bash.inc' ]; then . '/dors/capra_lab/data/gnomad/google-cloud-sdk/path.bash.inc'; fi

# The next line enables shell command completion for gcloud.
# if [ -f '/dors/capra_lab/data/gnomad/google-cloud-sdk/completion.bash.inc' ]; then . '/dors/capra_lab/data/gnomad/google-cloud-sdk/completion.bash.inc'; fi
##############################################################################

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/mothcw/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/mothcw/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/home/mothcw/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/mothcw/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<


###################################################################3
# END of Environment setup
###################################################################3


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
echo Halted vustruct_flask | mail -s 'Urgent Flask halt' chris.moth@vanderbilt.edu
sleep 60
done

