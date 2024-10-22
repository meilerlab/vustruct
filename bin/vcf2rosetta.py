#!/usr/bin/env python
#
# Project        : PSB Pipeline - support for MutationExplorer Server
                   Deployed by Rene Staritzbichler Danie Wiegreffe et al
                   https://www.biorxiv.org/content/10.1101/2023.03.23.533926v1
# Filename       : vcf2rosetta.py
# Author         : Chris Moth
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : chris.moth@vanderbilt.edu
# Date           : 2022 December
#
#=============================================================================#
"""vcf2rosetta.py is similar to the vcf2missense.py used generally in the
VUstruct pipeline.  For the MutationExplorer in Leipzig, I modifed
vcf2missense.py to output alpha-fold model names, focus on the canoncial 
uniprot transcript, and few other things quickly revealed by
$ diff vcf2rosetta.py vcf2missense.py
"""

# See main check for cmd line parsing
import argparse,configparser
import traceback
import sys,os,csv,time,pdb,glob,gzip,shutil,getpass
import subprocess as sp
from multiprocessing import cpu_count
import pandas as pd
import numpy as np
import copy
import vcf
# from warnings import filterwarnings,resetwarnings
from lib import PDBMapVEP
# from lib import PDBMapVCF
from lib import PDBMapProtein
from lib import PDBMapAlignment,PDBMapData
from lib import PDBMapModel
from lib import PDBMapAlphaFold
from lib.PDBMapVisualize import PDBMapVisualize
from lib import amino_acids
import logging
from logging.handlers import RotatingFileHandler
from logging import handlers
from psb_shared import psb_config
import pandas as pd


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
cmdline_parser.add_argument("--liftover",action='store_true',help="Pre-convert VCF file to GRCh38 from GRCh37 - hardcoded for capra lab")
args = cmdline_parser.parse_args()

required_config_items = ['vep','vep_cache_dir','idmapping']
config,config_dict = psb_config.read_config_files(args,required_config_items)
PDBMapProtein.load_idmapping(config_dict['idmapping'])

pdbmap_AlphaFold = PDBMapAlphaFold(config_dict)

pdbmap_vep = PDBMapVEP(config_dict)

vcf_input_filename = args.vcffile # Will change if we do a liftover

if args.liftover:
    LOGGER.info("--liftover set on command line")
    vcf_file_compressed = None  # <- This works great with vcf.Reader() calls involving .gz (compressed automatically) and all else (not compressed)
    if args.vcffile.endswith(
            ".bgz"):  # Why does Gnomad use .bgz for .gz files?  So irritating!  But, we override and all works out
        vcf_file_compressed = True



    grch37_original_records = []
    # Annoyingly, at apppears vcf.Reader never closes the file....
    grch37_vcf_reader =  vcf.Reader(filename=args.vcffile,
                      prepend_chr=False, #,  # (not 'liftover' in vcf_filename),
                      compressed=vcf_file_compressed,
                      encoding='UTF-8')

    liftover_input_df = pd.DataFrame(columns=['chrom','pos'])
    for vcf_record in grch37_vcf_reader:
        grch37_original_records.append(vcf_record)
        record_as_series = pd.Series({
            'chrom': 'chr' + str(vcf_record.CHROM), # vcf inputs lack CHR prefix.  .BED requires it for liftover
            'pos': int(vcf_record.POS)

        })
        liftover_input_df = pd.concat([liftover_input_df,record_as_series.to_frame().T],ignore_index=True)


    liftover_input_df.insert(1, 'bed_zerobased_pos', liftover_input_df["pos"] - 1)
    liftover_input_df.insert(3, 'original_lineno', liftover_input_df.index)

    vcffile_without_extension = os.path.splitext(args.vcffile)[0]
    liftover_bed_input = vcffile_without_extension + "_grch37_input.bed"
    liftover_output_lifted = vcffile_without_extension + "_grch38_lifted.bed"
    liftover_output_unlifted = vcffile_without_extension + "_grch38_unlifted.bed"
    liftover_input_df.to_csv(liftover_bed_input,sep='\t',header=False,index=False)
    import subprocess as sp
    # Ugly - but we need to find liftover either in a docker container or in the Vandy accre filesystem
    liftover_files = ('/psbadmin/liftover/liftOver','/psbadmin/liftover/hg19ToHg38.over.chain.gz')
    if not os.path.exists(liftover_files[0]) or not os.path.exists(liftover_files[1]):
        liftover_files = ('/dors/capra_lab/bin/liftOver','/dors/capra_lab/data/ucsc/liftOver/hg19ToHg38.over.chain.gz')
    liftover_arg_list = [liftover_files[0],
         liftover_bed_input,
         liftover_files[1],
         liftover_output_lifted,
         liftover_output_unlifted]

    LOGGER.info("Running liftOver to convert to GRCh38 positions:\n%s",str(' '.join(liftover_arg_list)))

    completed_process = sp.run(liftover_arg_list,
         capture_output=True,check=True,text=True)

    completed_process.check_returncode()

    try:
        df_lifted_grch38 = pd.read_csv(liftover_output_lifted,header=None,sep='\t',names=["chrom","bed_zerobased_pos_hg38","pos",'original_lineno'])
        LOGGER.info("%d genomic postions lifted to GRCh38",df_lifted_grch38.shape[0])
        # df_lifted_grch38 = df_lifted_grch38.drop('bed_zerobased_pos_hg38',axis='columns').set_index('original_lineno',drop=False)
    except pd.errors.EmptyDataError:
        LOGGER.critical("Liftover to GRCh38 created ntull lifted output")

    try:
        df_unlifted = pd.read_csv(liftover_output_unlifted,header=None,sep='\t',comment='#')
        LOGGER.critical("%d genomic postions could not be lifted to GRCh38",df_unlifted.shape[0])
    except pd.errors.EmptyDataError:
        LOGGER.info("All positions lifted to GRCh38 successfully")

    vcf_input_filename = vcffile_without_extension + "_grch38_lifted.vcf"


    with open(vcf_input_filename,'w') as vcf_lifted_to_Grch38:
        vcf_writer = vcf.Writer(vcf_lifted_to_Grch38,grch37_vcf_reader)
        for original_lineno,row in df_lifted_grch38.iterrows():
            record = copy.deepcopy(grch37_original_records[original_lineno])
            assert record.CHROM == row['chrom'][3:] # Skip opening CHR text
            record.POS = row['bed_zerobased_pos_hg38'] + 1
            vcf_writer.write_record(record)



vcf_reader =  pdbmap_vep.vcf_reader_from_file_supplemented_with_vep_outputs(vcf_input_filename)

raw_missense_df = pd.DataFrame(columns=['gene','chrom','pos','change','transcript','unp','refseq','mutation','alphafold_id'])

for vcf_record in pdbmap_vep.yield_completed_vcf_records(vcf_reader):
    for CSQ in vcf_record.CSQ:
        unp = PDBMapProtein.enst2unp(CSQ['Feature'])
        if unp and type(unp) is list: # << This is typical
            unp = unp[0]
        refseq = "NA"
        if not unp:
            LOGGER.warning("%-10s No uniprot ID for Ensembl transcript %s",
                           CSQ['SYMBOL'] + ':' if CSQ['SYMBOL'] else '', CSQ['Feature'])
            continue
        if not PDBMapProtein.isCanonicalByUniparc(unp):
            LOGGER.info("Skipping non-canonical unp %s" % unp)
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

        model_seq_start, model_seq_end, alphafold_id = pdbmap_AlphaFold.best_covering_model(unp,2700,int(CSQ['Protein_position']))
        next_df_row = pd.DataFrame.from_dict(
            [{'gene': CSQ['SYMBOL'],
             'chrom': vcf_record.CHROM,
             'pos': vcf_record.POS,
             'change': vcf_record.REF + '/' + vcf_record.ALT[0].sequence,
             'transcript': CSQ['Feature'],
             'unp': unp,
             'refseq': refseq,
             'mutation': variant_aa,
             'alphafold_id': alphafold_id
             }]
            )
        raw_missense_df = pd.concat([raw_missense_df,next_df_row],ignore_index=True)

raw_missense_df.to_csv(args.vcffile.split('.')[0]+'_missense.csv_withduplicates',sep=',')
df_without_duplicates = raw_missense_df.drop_duplicates(['gene','unp','refseq','mutation'])
raw_missense_df = raw_missense_df.set_index(['gene','unp','refseq','mutation']).sort_index()
for index,row in df_without_duplicates.iterrows():
    variant_index = (row['gene'],row['unp'],row['refseq'],row['mutation'])
    rows_with_various_transcripts = raw_missense_df.loc[[variant_index]]
    transcript_list = rows_with_various_transcripts['transcript'].tolist()
    # print("for variant_index %s rows are %s"%(str(variant_index),str(rows_with_various_transcripts)))
    # import pdb; pdb.set_trace()
    df_without_duplicates.at[index,'transcript'] = ';'.join(transcript_list)

final_missense_filename = args.vcffile.split('.')[0]+'_missense.csv'
df_without_duplicates.to_csv(final_missense_filename,sep=',')
# df = pd.DataFrame(columns=['gene','chrom','pos','transcript','unp','refseq','mutation'])
LOGGER.info("%d rows written to %s" , len(df_without_duplicates) , final_missense_filename)
