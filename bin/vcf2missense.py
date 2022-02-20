#!/usr/bin/env python
#
# Project        : PDBMap
# Filename       : PDBMap.py
# Author         : R. Michael Sivley - refactored by Chris Moth 2019 November
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu  chris.moth@vanderbilt.edu
# Date           : 2014-02-17                  2019-11-19
#
#=============================================================================#
"""PDBMap is a command-line interface for loading data
and providing access to the PDBMap* library modules.
PDBMap library modules map codons in the human exome to
their corresponding amino acids in known and predicted
protein structures. Using this two-way mapping, PDBMap
allows for the mapping of genomic annotations onto protein
structures and the mapping of protein structural annotation 
onto genomic elements. The methods to build and utilize this
tool may be called from this master class."""

# See main check for cmd line parsing
import argparse,configparser
import traceback
import sys,os,csv,time,pdb,glob,gzip,shutil,getpass
import subprocess as sp
from multiprocessing import cpu_count
import pandas as pd
import numpy as np
# from warnings import filterwarnings,resetwarnings
from lib import PDBMapVEP
# from lib import PDBMapVCF
from lib import PDBMapProtein
from lib import PDBMapAlignment,PDBMapData
from lib import PDBMapModel, PDBMapSwiss
from lib.PDBMapVisualize import PDBMapVisualize
from lib import amino_acids
import logging
from logging.handlers import RotatingFileHandler
from logging import handlers
from psb_shared import psb_config



sh = logging.StreamHandler()
LOGGER = logging.getLogger()
LOGGER.addHandler(sh)

log_format_string = '%(asctime)s %(levelname)-4s [%(filename)16s:%(lineno)d] %(message)s'
date_format_string = '%H:%M:%S'
log_formatter = logging.Formatter(log_format_string,date_format_string)

LOGGER.setLevel(logging.DEBUG)
sh.setLevel(logging.INFO)
sh.setFormatter(log_formatter)

cmdline_parser = psb_config.create_default_argument_parser(__doc__,os.path.dirname(os.path.dirname(__file__)))
cmdline_parser.add_argument("vcffile",type=str,metavar="FILE",help="filename in vcf format",default=os.path.basename(os.getcwd())+".vcf",nargs='?')

args = cmdline_parser.parse_args()

required_config_items = ['vep','vep_cache_dir','idmapping']
config,config_dict = psb_config.read_config_files(args,required_config_items)
PDBMapProtein.load_idmapping(config_dict['idmapping'])
pdbmap_vep = PDBMapVEP(config_dict)

vcf_reader =  pdbmap_vep.vcf_reader_from_file_supplemented_with_vep_outputs(args.vcffile)

df = pd.DataFrame(columns=['gene','chrom','pos','transcript','unp','refseq','mutation'])

for vcf_record in pdbmap_vep.yield_completed_vcf_records(vcf_reader):
    for CSQ in vcf_record.CSQ:
        unp = PDBMapProtein.enst2unp(CSQ['Feature'])
        if type(unp) is list: # << This is typical
            unp = unp[0]
        refseq = "NA"
        if not unp:
            LOGGER.warning("No uniprot ID for Ensembl transcript %s"%CSQ['Feature'])
            continue
        if unp:
            refseq = PDBMapProtein.unp2refseqNT(unp)
        if type(refseq) is list and refseq:
            refseq = refseq[0]
        else:
            refseq = "NA"
        variant_aa ="%s%s%s"%(CSQ['Ref_AminoAcid'],CSQ['Protein_position'],CSQ['Alt_AminoAcid'])
        if CSQ['Ref_AminoAcid'] == CSQ['Alt_AminoAcid']:
            LOGGER.warning("Excluding %s VEP-predicted synonymous variant %s"%(unp,variant_aa))
            continue
        df = df.append(
            {'gene': CSQ['SYMBOL'],
             'chrom': vcf_record.CHROM,
             'pos': vcf_record.POS,
             'transcript': CSQ['Feature'],
             'unp': unp,
             'refseq': refseq,
             'mutation': variant_aa
            },ignore_index=True)

df.to_csv(args.vcffile.split('.')[0]+'_missense.csv_withduplicates',sep=',')

df_without_duplicates = df.drop_duplicates(['gene','unp','refseq','mutation'])
df = df.set_index(['gene','unp','refseq','mutation']).sort_index()
for index,row in df_without_duplicates.iterrows():
    variant_index = (row['gene'],row['unp'],row['refseq'],row['mutation'])
    rows_with_various_transcripts = df.loc[[variant_index]]
    transcript_list = rows_with_various_transcripts['transcript'].tolist()
    # print("for variant_index %s rows are %s"%(str(variant_index),str(rows_with_various_transcripts)))
    # import pdb; pdb.set_trace()
    df_without_duplicates.at[index,'transcript'] = ';'.join(transcript_list)

df_without_duplicates.to_csv(args.vcffile.split('.')[0]+'_missense.csv',sep=',')
# df = pd.DataFrame(columns=['gene','chrom','pos','transcript','unp','refseq','mutation'])
