#!/bin/bash

# Failure of bash script should cause it to halt
set -e

# The PSB Pipeline requires a number of add-on modules to the 
# Anaconda base install

# To avoid LOTS of trouble:
# pip install --upgrade setuptools
# pip install --upgrade pip

# While anaconda includes numy and pandas default, raw OS installs of python3
# typically don't
pip install numpy
pip install pandas

# apScheduler allows the flask web front end to perform web updates 
# at fixed intervals
pip install apscheduler

# Biopython requires, and will include, numpy library
# Note that biopython may require manual patches due to bugs
pip install biopython

pip install mysql-connector-python
pip install mysqlclient
pip install argparse
pip install WeasyPrint

echo 'CAUTION.  Installation of PyVCF at this point in container construction'
echo 'May cause problems with python "import vcf" from container'
echo 'We are using a local PyVCF Python3 copy in pdbmap until figured out'
# pip install PyVCF

# Newer versions of networkx have incompatible data file formats
# Await news from Souhrid pip install networkx==1.9
pip install fuzzywuzzy

# python-Levenshtein may be already installed - no worries
pip install python-Levenshtein

# Seems jinja2 also now installed already - no worries
pip install jinja2
pip install tabulate

# Digenic analysis code requires plotly
pip install plotly

# We use openpyxl to support generation of case summary spreadsheet
pip install openpyxl
