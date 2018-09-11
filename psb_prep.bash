echo "Preparing PDB Pipeline PATH and environment settings"
echo "This file should only be sourced from a bash or sh shell"
echo 
# Directories that must be added to your path
newPATH=/dors/capra_lab/users/psbadmin/anaconda2/bin
newPATH=$newPATH:/dors/capra_lab/users/psbadmin/bin
newPATH=$newPATH:/dors/capra_lab/users/psbadmin/pdbmap
newPATH=$newPATH:/dors/capra_lab/users/psbadmin/pathprox
newPATH=$newPATH:/dors/capra_lab/opt/ensembl-tools-release-87/scripts/variant_effect_predictor
newPATH=$newPATH:/dors/capra_lab/opt/ensembl-tools/release-87/scripts/id_history_converter
newPATH=$newPATH:/dors/capra_lab/bin/vcftools/bin
newPATH=$newPATH:/dors/capra_lab/bin/vcftools/perl
echo -ne "Prepending to PATH:\n\t"
echo $newPATH | sed 's/:/\n\t/g'
PATH=$newPATH:$PATH
export PATH
# echo PATH=$PATH

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
OPT=/dors/capra_lab/opt
newPERL5LIB=$OPT/vcftools_0.1.12b/perl
newPERL5LIB=$newPERL5LIB:$OPT/bioperl-live
newPERL5LIB=$newPERL5LIB:$OPT/src/ensembl/modules
newPERL5LIB=$newPERL5LIB:$OPT/src/ensembl-compara/modules
newPERL5LIB=$newPERL5LIB:$OPT/src/ensembl-variation/modules
newPERL5LIB=$newPERL5LIB:$OPT/src/ensembl-funcgen/modules
echo -ne "Prepending to PERL5LIB:\n\t"
echo $newPERL5LIB | sed 's/:/\n\t/g'
export PERL5LIB=$newPERL5LIB:${PERL5LIB}

echo
# newPYTHONPATH=/dors/capra_lab/users/psbadmin/pdbmap:/dors/capra_lab/psbadmin/bin
newPYTHONPATH=/dors/capra_lab/users/psbadmin/pdbmap:/dors/capra_lab/users/psbadmin/bin
echo -ne "Prepending to PYTHONPATH:\n\t"
echo $newPYTHONPATH | sed 's/:/\n\t/g'
export PYTHONPATH=$newPYTHONPATH:$PYTHONPATH

echo
export UDN=/dors/capra_lab/projects/psb_collab/UDN
echo "\$UDN set to $UDN"

export PS1='\r\n\[$(/usr/bin/tput setaf 15; /usr/bin/tput setab 4)\]   PSB Pipeline    host:\h  user:\u         \[$(/usr/bin/tput sgr0)\]\r\n\w\$ '

fi
