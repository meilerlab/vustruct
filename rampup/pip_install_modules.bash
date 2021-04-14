#!/bin/bash

# The PSB Pipeline requires a number of add-on modules to the 
# Anaconda base install

# To avoid LOTS of trouble:
pip install --upgrade setuptools
pip install --upgrade pip

# As of March 2019 python 3, pandas seems to have been loaded by default
# pip install pandas

# Biopython requires, and will include, numpy library
# Note that biopython may require manual patches due to bugs
pip install biopython

pip install mysql-connector-python
pip install mysqlclient
pip install argparse
pip install WeasyPrint
pip install PyVCF

# Newer versions of networkx have incompatible data file formats
# Await news from Souhrid pip install networkx==1.9
pip install fuzzywuzzy

# python-Levenshtein may be already installed - no worries
pip install python-Levenshtein

# Seems jinja2 also now installed already - no worries
pip install jinja2
pip install tabulate
