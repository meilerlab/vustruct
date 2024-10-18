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
import errno
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
from lib import PDBMapModel, PDBMapSwiss
from lib.PDBMapVisualize import PDBMapVisualize
from lib import amino_acids
import logging
from logging.handlers import RotatingFileHandler
from logging import handlers
from psb_shared import psb_config
import pandas as pd
from vustruct import VUStruct

cmdline_parser = psb_config.create_default_argument_parser(__doc__,os.path.dirname(os.path.dirname(__file__)))
cmdline_parser.add_argument("vcffile",type=str,metavar="FILE",help="filename in vcf format",default=os.path.basename(os.getcwd())+".vcf",nargs='?')
cmdline_parser.add_argument("--liftover",action='store_true',help="Pre-convert VCF file to GRCh38 from GRCh37 - hardcoded for capra lab")
cmdline_parser.add_argument("project", type=str, help="Project ID (ex. UDN124356)",
                            default=os.path.basename(os.getcwd()), nargs='?')
args = cmdline_parser.parse_args()

# Now that we've added streamHandler, basicConfig will not add another handler (important!)
vustruct = VUStruct('preprocess', args.project,  __file__)
vustruct.stamp_start_time()
vustruct.initialize_file_and_stderr_logging(__file__)
LOGGER = logging.getLogger();

# Prior to getting going, we save the vustruct filename
# and then _if_ we die with an error, at least there is a record
# and psb_rep.py should be able to create a web page to that effect
vustruct.exit_code = 1
vustruct.input_filename = os.path.join('./',args.vcffile) # will change if we do a liftover?
vustruct.write_file()

required_config_items = ['vep','vep_cache_dir','idmapping']
config,config_dict = psb_config.read_config_files(args,required_config_items)
udn_root_directory = os.path.join(config_dict['output_rootdir'], config_dict['collaboration'])
collaboration_dir = os.path.join(udn_root_directory, args.project)

PDBMapProtein.load_idmapping(config_dict['idmapping'])
pdbmap_vep = PDBMapVEP(config_dict)


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

    # In bed format https://en.wikipedia.org/wiki/BED_(file_format)
    # the 4th column is a name of the line.  So, we want to simply use the 0-based
    # row index - but to get that output to the .bed file in the rightmost column
    # means we cannot use the usual to_csv (index...) option
    liftover_input_df.insert(3, 'original_lineno', liftover_input_df.index)

    vcffile_without_extension = os.path.splitext(args.vcffile)[0]
    liftover_bed_input = vcffile_without_extension + "_grch37_input.bed"
    liftover_output_lifted = vcffile_without_extension + "_grch38_lifted.bed"
    liftover_output_unlifted = vcffile_without_extension + "_grch38_unlifted.bed"

    # The .bed input file that is output will have chr*, 0-based chromStart, and chromEnd (non-inclusive)
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

        # Critically, make sure that the index is NOT the 0..Count(lifted) - but rather the ORIGINAL
        # subscripts from the input file.  Thus, there will be skips in the index from records that were
        # NOT lifted
        df_lifted_grch38.set_index('original_lineno',inplace=True)
        # df_lifted_grch38 = df_lifted_grch38.drop('bed_zerobased_pos_hg38',axis='columns').set_index('original_lineno',drop=False)
    except pd.errors.EmptyDataError:
        LOGGER.critical("Liftover to GRCh38 created ntull lifted output")

    try:
        df_unlifted = pd.read_csv(liftover_output_unlifted,header=None,sep='\t',names=["chrom","bed_zerobased_pos_hg38","pos",'original_lineno'])
        # df_unlifted = pd.read_csv(liftover_output_unlifted,header=None,sep='\t',comment='#')
        LOGGER.critical("%d genomic postions could not be lifted to GRCh38",df_unlifted.shape[0])
        df_unlifted.set_index('original_lineno',inplace=True)
    except pd.errors.EmptyDataError:
        LOGGER.info("All positions lifted to GRCh38 successfully")

    vustruct.input_filename = vcffile_without_extension + "_grch38_lifted.vcf"


    # Now open the final GRCh38 vcf file which will be fed through the VEP
    with open(vustruct.input_filename,'w') as vcf_lifted_to_Grch38:
        vcf_writer = vcf.Writer(vcf_lifted_to_Grch38,grch37_vcf_reader)

        # Liftover can send back all kinds of nonsense chromosomes
        # 2024-10-08 I give up on trying to process that with vep
        # and we stick to what we know must work
        valid_chrom_strings = set()
        for chr in range(1,24):
            valid_chrom_strings.add("chr%s" % chr)
        valid_chrom_strings.add("chrX")
        valid_chrom_strings.add("chrY")

        # We could sort by CHROM, otherwise VEP gets upset when processing the output
        # df_lifted_grch38_sorted_by_chrom = df_lifted_grch38.sort_values('chrom')
        # But when vep reuiqres sorting that really makes for problems with matching up liftover records
        for original_lineno,lifted_row in df_lifted_grch38.iterrows():
            record = copy.deepcopy(grch37_original_records[original_lineno])

            lifted_chrom = lifted_row['chrom']            
            pos_from_lifted = lifted_row['bed_zerobased_pos_hg38'] + 1
            if lifted_chrom not in valid_chrom_strings:
                LOGGER.warning("Liftover returned invalid chrom string %s:%s from %s:%s bed lineno=%s" % (lifted_chrom, pos_from_lifted, record.CHROM, record.POS, original_lineno))
                continue

            chrom_from_lifted = lifted_row['chrom'][3:]

            """ if lifted_row['chrom'].startswith('chrUn'):
                LOGGER.warning("Unknown new CHROM %s:%s lifted from %s:%s will not be saved" % (
                    lifted_row['chrom'], pos_from_lifted, record.CHROM, record.POS))
                continue

            # A weird thing is that sometimes you get a very strange
            # strange new chr entry from liftover.  For example
            # chrUn1_KI270766v1_alt.   I'm not sure - but I now
            # trum everthing at _ and following in these cases.
            chrom_underscore_pos = chrom_from_lifted.find('_')
            if chrom_underscore_pos > 0:
                # chrom_from_lifted = chrom_from_lifted[0:chrom_underscore_pos]
                # In these weird examples, the CHR can get reassigned
                # record.CHROM = 'chr%s' % chrom_from_lifted[3:]
                LOGGER.warning("New CHROM %s:%s replacig %s:%s" % (chrom_from_lifted, pos_from_lifted, record.CHROM, record.POS))
                record.CHROM = chrom_from_lifted
            elif record.CHROM != chrom_from_lifted:
                LOGGER.warning("New CHROM %s:%s replacig %s:%s" % (chrom_from_lifted, pos_from_lifted, record.CHROM, record.POS))
                record.CHROM = chrom_from_lifted
            """

            record.CHROM = chrom_from_lifted
            record.POS = pos_from_lifted
            vcf_writer.write_record(record)


vcf_reader =  pdbmap_vep.vcf_reader_from_file_supplemented_with_vep_outputs(vustruct.input_filename,coding_only = False)

raw_vustruct_df = pd.DataFrame(columns=['gene','chrom','pos','change','transcript','unp','refseq','mutation'])

vep_consequences_df = pd.DataFrame(columns=['chrom','pos','ref','alt','consequence','transcript','unp','mutation','gene','BIOTYPE','SIFT','PolyPhen','refseq'])

for vcf_record in pdbmap_vep.yield_completed_vcf_records(vcf_reader):
    for CSQ in vcf_record.CSQ:
        ensembl_transcript_id_from_vep = None
        if 'Feature' in CSQ and CSQ['Feature']:
            ensembl_transcript_id_from_vep = CSQ['Feature']
            unp = PDBMapProtein.enst2unp(ensembl_transcript_id_from_vep)
        else:
            unp = ""
        if unp and type(unp) is list: # << This is typical
            unp = unp[0]
        if not unp: # Don't leave unp set to []
            unp = ''
        refseq = "NA"

         

        if 'Alt_AminoAcid' in CSQ and CSQ['Alt_AminoAcid']:
            variant_aa ="%s%s%s"%(CSQ['Ref_AminoAcid'],CSQ['Protein_position'],CSQ['Alt_AminoAcid']) 
        else :
            variant_aa = ''

        # In our missense file, we ideally want ALL The cross-referenced Ensembl transcript IDs.
        # So that is to say...  
        # CHROM/POS -> ONE ensembl transcript ID -> Uniprot -> LIST of ENSEMBL transcript IDs with version #s
        # Important - from uniprot cross-ref, all ENSEMBL transcript IDs with the same uniprot sequence
        # We can start this list with the versionless 
        # ensembl_transcript_ids_from_uniprot = [ensembl_transcript_id_from_vep]
        if unp:
            # ensembl_transcript_ids_from_uniprot = PDBMapProtein.unp2enst(unp)
            refseq = PDBMapProtein.unp2refseqNT(unp)
        if type(refseq) is list and refseq:
            refseq = refseq[0]
        else:
            refseq = "NA"

        # Let's presume vep has emitted a missense variant, although it many cases it will not have
        add_to_vustruct = True

        vep_results_row_dict = {
             'gene': CSQ['SYMBOL'],
             'chrom': vcf_record.CHROM,
             'pos': vcf_record.POS,
             'ref': vcf_record.REF,
             'alt': vcf_record.ALT[0].sequence,
             'consequence': CSQ['Consequence'],
             'transcript': CSQ['Feature'],
             'unp': unp,
             'mutation': variant_aa,
             'BioType': CSQ['BIOTYPE'],
             'SIFT': CSQ['SIFT'],
             'PolyPhen': CSQ['PolyPhen'],
             'refseq': refseq
             }

        if add_to_vustruct and not unp:
            vep_results_row_dict['consequence'] = 'Skipped: Unable to cross-references %s to Uniprot ID. %s' % (
                CSQ['Feature'],
                CSQ['Consequence']
                )
            LOGGER.warning("%s %s %s %s: Excluding VEP-predicted %s %s.  No cross-reference from %s to Uniprot ID"% (
                vcf_record.CHROM, 
                vcf_record.POS, 
                vcf_record.REF, 
                vcf_record.ALT[0].sequence, 
                CSQ['Feature'],
                variant_aa,
                CSQ['Consequence']
                ))
            LOGGER.warning("%s", vep_results_row_dict['consequence'])
            add_to_vustruct = False


        if add_to_vustruct and 'Consequence' in CSQ and CSQ['Consequence'] != 'missense_variant':
            LOGGER.warning("%s %s %s %s: Excluding %s VEP-predicted %s %s"%(vcf_record.CHROM, vcf_record.POS, vcf_record.REF, vcf_record.ALT[0].sequence, unp,CSQ['Consequence'], variant_aa))
            vep_results_row_dict['consequence'] = 'Skipped: %s' % CSQ['Consequence']
            add_to_vustruct = False


        if add_to_vustruct and CSQ['Ref_AminoAcid'] == CSQ['Alt_AminoAcid']:
            LOGGER.warning("Excluding %s VEP-predicted synonymous variant %s"%(unp,variant_aa))
            vep_results_row_dict['consequence'] = "Skipping synonymous"
            add_to_vustruct = False

        # We store all consequences from vep here, even if we skip a missense row
        vep_consequences_df = pd.concat([vep_consequences_df,pd.DataFrame.from_dict([vep_results_row_dict])],ignore_index=True)

        if add_to_vustruct:
            next_vustruct_df_row = pd.DataFrame.from_dict(
                [{'gene': CSQ['SYMBOL'],
                 'chrom': vcf_record.CHROM,
                 'pos': vcf_record.POS,
                 'change': vcf_record.REF + '/' + vcf_record.ALT[0].sequence,
                 'transcript': PDBMapProtein.versioned_ensembl_transcript(ensembl_transcript_id_from_vep),
                 'unp': unp,
                 'refseq': refseq,
                 'mutation': variant_aa
                 }]
                )
            raw_vustruct_df = pd.concat([raw_vustruct_df,next_vustruct_df_row],ignore_index=True)

raw_vustruct_df.to_csv(args.vcffile.split('.')[0]+'_vustruct.csv_withduplicates',sep=',')
df_without_duplicates = raw_vustruct_df.drop_duplicates(['gene','unp','refseq','mutation'])
raw_vustruct_df = raw_vustruct_df.set_index(['gene','unp','refseq','mutation']).sort_index()
for index,row in df_without_duplicates.iterrows():
    variant_index = (row['gene'],row['unp'],row['refseq'],row['mutation'])
    rows_with_various_transcripts = raw_vustruct_df.loc[[variant_index]]
    transcript_list = rows_with_various_transcripts['transcript'].tolist()
    # print("for variant_index %s rows are %s"%(str(variant_index),str(rows_with_various_transcripts)))
    # import pdb; pdb.set_trace()
    df_without_duplicates.at[index,'transcript'] = ';'.join(transcript_list)

final_vustruct_filename = args.vcffile.split('.')[0]+'_vustruct.csv'

# Renumber to get rid of gaps in the numbering which result from transcript matches
df_without_duplicates = df_without_duplicates.reset_index(drop=True)
# For writing final dataframe, index from 1 (not zero)
df_without_duplicates.index = df_without_duplicates.index + 1
df_without_duplicates.to_csv(final_vustruct_filename,sep=',', index_label='index')
vep_consequences_df.to_csv('vep_consequences_df.csv',sep=',', index_label='index')
# df = pd.DataFrame(columns=['gene','chrom','pos','transcript','unp','refseq','mutation'])
LOGGER.info("%d rows written to %s" , len(df_without_duplicates) , final_vustruct_filename)
vustruct.exit_code = 0
vustruct.stamp_end_time()
vustruct.write_file()
