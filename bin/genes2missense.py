#!/usr/bin/env python
#
# Project        : PDB Pipeline
# Filename       : genes2missense.py
# Author         : Chris Moth 2021 March
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : chris.moth@vanderbilt.edu
# Date           : 2021-03-08
#
#=============================================================================#
"""Parse a simple comma-delimted file in format:
Gene,Variant
Gene1,A1234B
Gene2,C4567D"""


# See main check for cmd line parsing
import argparse,configparser
import traceback
import sys,os,csv,time,pdb,glob,gzip,shutil,getpass
import subprocess as sp
from multiprocessing import cpu_count
import pandas as pd
import numpy as np
from warnings import filterwarnings,resetwarnings
from lib import PDBMapProtein
from lib import PDBMapTranscriptUniprot
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
cmdline_parser.add_argument("genefile",type=str,metavar="FILE",help="filename containing rows of gene,variant",default=os.path.basename(os.getcwd())+"_genes.csv",nargs='?')

args = cmdline_parser.parse_args()

required_config_items = ['idmapping']
config,config_dict = psb_config.read_config_files(args,required_config_items)
PDBMapProtein.load_idmapping(config_dict['idmapping'])

gene_df = pd.read_csv(args.genefile)
assert len(gene_df) > 1,"No data rows were found in %s"%args.genefile
assert 'gene' in gene_df.columns,"A header row in %s is required, and must contain 'gene' and 'variant'"
assert 'variant' in gene_df.columns,"A header row in %s is required, and must contain 'gene' and 'variant'"

missense_df = pd.DataFrame(columns=['gene','unp','refseq','mutation'])

for index,row in gene_df.iterrows():
    gene_id = PDBMapProtein.gene2gene_id(row['gene'])
    if not gene_id:
        fail_message = "Gene %s in your input file lacks a corresponding uniprot identifier in the idmapping file"%row['gene']
        LOGGER.critical(fail_message)
        sys.exit(fail_message)
    canonical_unp = PDBMapProtein.hgnc2unp(row['gene'])
    LOGGER.info("Gene %s maps to gene_id: %s and canonical unp: %s",row['gene'],gene_id,canonical_unp)
    canonical_unp_transcript = PDBMapTranscriptUniprot(canonical_unp)
    LOGGER.info("Transcript for %s has len %d",canonical_unp_transcript.id,len(canonical_unp_transcript.aa_seq))

    try:
        ref_aa = row['variant'][0]
        alt_aa = row['variant'][-1]
        variant_pos = int(row['variant'][1:-1])
    except:
        sys.exit("Variant amino acid format must be A1234B.  You have [%s] in row index %d"%(row['variant'],variant_pos))

    assert ref_aa == canonical_unp_transcript.aa_seq[variant_pos-1],"You provided variant %s but position %d is %s not %s"%(
        row['variant'],variant_pos,canonical_unp_transcript.aa_seq[variant_pos-1],ref_aa)

    refseq = PDBMapProtein.unp2refseqNT(canonical_unp)

    if type(refseq) is list and refseq:
        refseq = refseq[0]
    else:
        refseq = "NA"

    missense_df = missense_df.append(
        {'gene': row['gene'],
         'transcript': PDBMapProtein.unp2enst(canonical_unp)[0],
         'unp': canonical_unp,
         'refseq': refseq,
         'mutation': row['variant']
        },ignore_index=True)

missense_filename = os.path.basename(os.getcwd())+"_missense.csv"
LOGGER.info("Writing %d rows of missense variants to %s",len(missense_df),missense_filename)
missense_df.to_csv(missense_filename)

