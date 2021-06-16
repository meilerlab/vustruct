echo "Preparing PDB Pipeline PATH and environment settings"
echo "This file should only be sourced from a bash shell"
echo
genome='GRCh38'
if [[ $0 != $BASH_SOURCE ]]; then
PIPELINE_ROOT=$(readlink -f `dirname $BASH_SOURCE`)
echo "Script in $PIPELINE_ROOT was correctly sourced"
#
# Sometimes accre changes the first path component from /dors to /dors1
# So we compare the _next_/trailing path components to determine
# if we are isProduction mode.....
if [[ ${PIPELINE_ROOT:(-25)} = '/capra_lab/users/psbadmin' ]]; then
isProduction=true
else
isProduction=false
false
fi
if $isProduction ; then
echo "PRODUCTION INITIALIZATION"
else
echo "** DEVELOPMENT (NOT PRODUCTION) ENVIRONMENT INITIALIZING **"
fi
else
echo "Script cannot be run as you invoked it.  You must instead 'source $(readlink -f $0)'"
exit
fi

# When accessing genomic coordinates, we open the database specified in this configuration file
# This .conf file must be updated when new versions of the genome, or the Ensembl PERL API are loaded
if [[ $genome == 'GRCh37' ]]; then
echo " **********   OBSOLETE GRCh37 DETECTED ********** "
if [[ $- = *i* ]]
then 
# Interactive shell
while true; do
    read -p "Are you SURE you want to continue Y/N?" yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done
fi # End interactive shell
export ENSEMBL_REGISTRY=$PIPELINE_ROOT/pdbmap/EnsEMBL/ensembl_registry.conf
else
export ENSEMBL_REGISTRY=$PIPELINE_ROOT/pdbmap/EnsEMBL/ensembl_registry_GRCh38.conf
fi

echo $genome selected with ENSEMBL_REGISTRY=$ENSEMBL_REGISTRY

export LC_ALL="en_US.UTF-8"
# PERL binaries, and Ensembl are NOT copied to
# the development environment.  Not worth trouble
PSBADMIN_PERL_ROOT=/dors/capra_lab/users/psbadmin
newPATH=''
# Directories that must be added to your path
newPATH=$PIPELINE_ROOT/anaconda3/bin
newPATH=$newPATH:$PSBADMIN_PERL_ROOT/localperl/bin
newPATH=$newPATH:$PIPELINE_ROOT/bin
newPATH=$newPATH:$PIPELINE_ROOT/pdbmap
newPATH=$newPATH:$PIPELINE_ROOT/pathprox
newPATH=$newPATH:$PIPELINE_ROOT/htslib/bin
newPATH=$newPATH:/dors/capra_lab/bin/vcftools/bin
newPATH=$newPATH:/dors/capra_lab/bin/vcftools/perl
newPATH=$newPATH:/dors/capra_lab/users/psbadmin/psb_prep.bash
newPATH=$newPATH:/dors/capra_lab/bin/wkhtmltox/bin
echo -ne "Prepending to PATH:\n\t"
echo $newPATH | sed 's/:/\n\t/g'
PATH=$newPATH:$PATH
export PATH
# echo PATH=$PATH

psbbin=$PIPELINE_ROOT/bin
# Skip a line
echo
# Make _SURE_ that user has Perl 5 in their environment
perl --version | grep -q "perl 5"
if [ $? -ne 0 ]; then
  echo "Perl 5 is a pre-requisite for the psb pipeline"
  echo "Pipeline is certain to fail with your version of perl"
else
echo "Perl 5 detected"
# Configure the perl environment for the Ensembl API and the Variant Effect Predictor
# OPT=/dors/capra_lab/opt
# newPERL5LIB=$OPT/vcftools_0.1.12b/perl
newPERL5LIB=$PSBADMIN_PERL_ROOT/localperl/lib/perl5
# newPERL5LIB=$newPERL5LIB:$PSBADMIN_PERL_ROOT/localperl/lib/perl5/x86_64-linux/
newPERL5LIB=$newPERL5LIB:$PSBADMIN_PERL_ROOT/ensembl/ensembl/modules
newPERL5LIB=$newPERL5LIB:$PSBADMIN_PERL_ROOT/ensembl/ensembl-variation/modules
newPERL5LIB=$newPERL5LIB:$PSBADMIN_PERL_ROOT/ensembl/ensembl-compara/modules
newPERL5LIB=$newPERL5LIB:$PSBADMIN_PERL_ROOT/ensembl/ensembl-funcgen/modules
newPERL5LIB=$newPERL5LIB:$PSBADMIN_PERL_ROOT/ensembl/ensembl-tools/modules
newPERL5LIB=$newPERL5LIB:$PSBADMIN_PERL_ROOT/ensembl/ensembl-io/modules
newPERL5LIB=$newPERL5LIB:$PSBADMIN_PERL_ROOT/ensembl/BioPerl-1.6.924
unset OPT
echo -ne "Replacing PERL5LIB:\n\t"
echo $newPERL5LIB | sed 's/:/\n\t/g'
# export PERL5LIB=$newPERL5LIB:${PERL5LIB} # We used to prepend - but perl is just crazy now with Anaconda
export PERL5LIB=$newPERL5LIB

echo
newPYTHONPATH=$PIPELINE_ROOT/pdbmap:$PIPELINE_ROOT/bin
echo -ne "Prepending to PYTHONPATH:\n\t"
echo $newPYTHONPATH | sed 's/:/\n\t/g'
export PYTHONPATH=$newPYTHONPATH:$PYTHONPATH


# Prevent loading of an individual users ~/.local type python - really a mess for pipeline!
export PYTHONNOUSERSITE=x

# Don't buffer output, to get better errors if things go south
export PYTHONUNBUFFERED=x

echo

if $isProduction; then
export UDN=/dors/capra_lab/projects/psb_collab/UDN
foreGround=15
backGround=4
if [ `whoami` == 'psbadmin' ]
then
backGround=1
fi
export PS1='\r\n\[$(/usr/bin/tput setaf '$foreGround'; /usr/bin/tput setab '$backGround')\]   PSB Pipeline    host:\h  user:\u         \[$(/usr/bin/tput sgr0)\]\r\n\w\$ '
else
export UDN=/dors/capra_lab/users/$USER/UDNtests
# export PS1='\[$(/usr/bin/tput setaf 15; /usr/bin/tput setab 1)\] DEV \[$(/usr/bin/tput sgr0)\] '$PS1
export PS1='DEV '$PS1
fi
echo "\$UDN=$UDN"

fi
if [[ $genome == 'GRCh37' ]]; then
echo *******  WARNING GRCh37 selected ${ENSEMBL_REGISTRY}
fi
echo "\$ENSEMBL_REGISTRY=${ENSEMBL_REGISTRY}"
