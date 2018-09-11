echo "Preparing PDB Pipeline PATH and environment settings"
echo "This file should only be sourced from a tcsh or csh shell"
echo "Your shell is "$SHELL
echo 
# Directories that must be added to your path
setenv newPATH /dors/capra_lab/users/psbadmin/anaconda2/bin
setenv newPATH $newPATH\:/dors/capra_lab/users/psbadmin/bin
setenv newPATH $newPATH\:/dors/capra_lab/users/psbadmin/pdbmap
setenv newPATH $newPATH\:/dors/capra_lab/users/psbadmin/pathprox
setenv newPATH $newPATH\:/dors/capra_lab/opt/ensembl-tools-release-87/scripts/variant_effect_predictor
setenv newPATH $newPATH\:/dors/capra_lab/opt/ensembl-tools/release-87/scripts/id_history_converter
setenv newPATH $newPATH\:/dors/capra_lab/bin/vcftools/bin
setenv newPATH $newPATH\:/dors/capra_lab/bin/vcftools/perl
echo -n "Prepending to PATH:\n\t"
echo $newPATH | sed 's/:/\n\t/g'
setenv PATH $newPATH\:$PATH
echo
# Make _SURE_ that user has Perl 5 in their environment
perl --version | grep -q "perl 5"
if ( $? != 0 ) then
  echo "Perl 5 is a pre-requisite for the psb pipeline"
  echo "Pipeline is certain to fail with your version of perl"
else
echo "Perl 5 detected"
# Configure the perl environment for the Ensembl API and the Variant Effect Predictor
set OPT=/dors/capra_lab/opt
setenv newPERL5LIB $OPT/vcftools_0.1.12b/perl
setenv newPERL5LIB $newPERL5LIB\:$OPT/bioperl-live
setenv newPERL5LIB $newPERL5LIB\:$OPT/src/ensembl/modules
setenv newPERL5LIB $newPERL5LIB\:$OPT/src/ensembl-compara/modules
setenv newPERL5LIB $newPERL5LIB\:$OPT/src/ensembl-variation/modules
setenv newPERL5LIB $newPERL5LIB\:$OPT/src/ensembl-funcgen/modules
echo -n "Prepending to PERL5LIB:\n\t"
echo $newPERL5LIB | sed 's/:/\n\t/g'
setenv PERL5LIB $newPERL5LIB\:$PERL5LIB
echo
# setenv newPYTHONPATH /dors/capra_lab/users/psbadmin/pdbmap\:/dors/capra_lab/psbadmin/bin
setenv newPYTHONPATH /dors/capra_lab/users/psbadmin/pdbmap
set newPYTHONPATH=/dors/capra_lab/users/psbadmin/pdbmap
echo -n "Prepending to PYTHONPATH:\n\t"
echo $newPYTHONPATH | sed 's/:/\n\t/g'
setenv PYTHONPATH $newPYTHONPATH\:$PYTHONPATH
set PYTHONPATH=$newPYTHONPATH\:$PYTHONPATH

echo
setenv UDN /dors/capra_lab/projects/psb_collab/UDN
set UDN=/dors/capra_lab/projects/psb_collab/UDN
echo \$UDN set to $UDN

set prompt="\n%S   PSB Pipeline  host:%m  user:%n   %s\n%/% "

endif
