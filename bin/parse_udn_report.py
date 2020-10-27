#!/usr/bin/env python
"""\
This script will parse a UDN Patient Report (excel format) and 
save all of the mutations to a structured .csv file which is
formatted for input to add_udn_report.py\
"""

import logging


import sys
import os

import re
import string
import pandas as pd
import json
from lib.amino_acids import longer_names
from lib import PDBMapProtein
from lib import PDBMapVEP
import unicodedata

# logging.basicConfig(format='%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s',
#    datefmt='%d-%m-%Y:%H:%M:%S',)
logging.basicConfig(format='%(levelname)-8s [%(filename)s:%(lineno)-3d] %(message)s', )

# Send INFO and higher logging information to stderr
# WARNING and higher log details are echoed to the log file
root = logging.getLogger()
root.setLevel(logging.INFO)

LOGGER = logging.getLogger()


def remove_unicode_control_characters(s):
    return "".join(ch for ch in s if unicodedata.category(ch)[0] != "C")


def remove_ascii_control_characters(s):
    return ''.join(c for c in s if (32 <= ord(c) < 128))


from psb_shared import psb_config

cmdline_parser = psb_config.create_default_argument_parser(__doc__,
                                                           os.path.dirname(os.path.dirname(__file__)))
cmdline_parser.add_argument("project", type=str, help="Project ID (ex. UDN124356)",
                            default=os.path.basename(os.getcwd()), nargs='?')

# parser.add_argument("udn_excel",type=str,help="Raw input UDN patient report (format: xls/xlsx)")
# parser.add_argument("udn_csv",type=str,help="Parsed output pipeline filename (e.g. filename.csv)")

args = cmdline_parser.parse_args()

required_items = ['output_rootdir', 'collaboration']

config, config_dict = psb_config.read_config_files(args, required_items)

udn_root_directory = os.path.join(config_dict['output_rootdir'], config_dict['collaboration'])
collaboration_dir = os.path.join(udn_root_directory, args.project)

udn_excel_filename = os.path.join(udn_root_directory, args.project, args.project + ".xlsx")
missense_csv_filename = os.path.join(udn_root_directory, args.project, args.project + "_missense.csv")
genes_txt_filename = os.path.join(udn_root_directory, args.project, args.project + "_genes.txt")
genes_json_filename = os.path.join(udn_root_directory, args.project, args.project + "_genes.json")
LOGGER.info('Parsing Excel file: %s' % udn_excel_filename)
LOGGER.info(' Writing csv file %s ' % missense_csv_filename)
LOGGER.info(' Writing txt file %s ' % genes_txt_filename)
LOGGER.info('Writing json file %s ' % genes_json_filename)

df = pd.read_excel(udn_excel_filename,header=0)
# Now read it again and map everything to strings
df = pd.read_excel(udn_excel_filename,header=0,converters={k: str for k in range(df.shape[1])})
# print df

LOGGER.info('Excel file read successfully with dimensions (pandas shape) of %s' % str(df.shape))
dfRows = df.shape[0]

genes = {}
gene_patterns = {'\u25cf\u25cb': 'Het',
                 '\u25cb\u25cb': 'WT',
                 '\u25cf\u25cf': 'Homo',
                 '\u25cfy': 'Het/Y',
                 '\u25cby': 'WT/Y',
                 '\u25cf?': 'Het/?'}

# Replaced with line below by ChrisMoth
# PDBMapProtein.load_idmapping("/dors/capra_lab/data/uniprot/idmapping/HUMAN_9606_idmapping_sprot.dat.gz")
# This really needs to come out of config file!!!!
idmappingFilename = config_dict[
    'idmapping']  # "/dors/capra_lab/users/mothcw/mydata/idmapping/HUMAN_9606_idmapping_sprot.dat.gz"
LOGGER.info("Loading uniprot idmapping from: %s" % idmappingFilename)
PDBMapProtein.load_idmapping(idmappingFilename)

GeneWordEncountered = False
if df.columns.values[0].strip() == 'Gene':
    logging.critical('The word "Gene" must not be in the top left header of the spreadsheet.' +
                     'It must appear later, after the initial row which includes a case name.')
    sys.exit(1)

# The strings Gene, Chr, Change, Effect, Proband, Mother, Father 
# Will likely be the "headers" for the gene dictionary.
# However, the structure is a little flexible to allow for UDN slop
# See below where this important list is initialized
HeadersList = []

# This loop is not Pythonic - but the input .xls file is very unstructured.  We need to skip rows
i = 0
csv_rows = []

genome = 'hg38'
while i < dfRows:
    row = df.iloc[i]
    # print i, row[0]
    try:
        if type(row[0]) == str and row[0].find('hg19') != -1:
            genome = 'hg19'
        if type(row[0]) == str and row[0].find('Gene') != -1:
            position_right1_down1 = str(df.iloc[i+1][1])
            if 'Position' in position_right1_down1 and 'hg19' in position_right1_down1:
               genome = 'hg19'
            GeneWordEncountered = True
            HeadersList = [str(x).strip().split()[0] for x in row]
            i += 1
            continue

    except (TypeError,IndexError):
        pass
    
    if not GeneWordEncountered:
        i += 1
        continue
    
    if GeneWordEncountered:  # Until we see "Gene" - we just keep looping
        try:
            possible_gene = row[0].strip()  # .encode('utf-8')
        except (TypeError, AttributeError):
            i += 1
            continue

        # if isinstance(possible_gene,str):
        #  cleaned_possible_gene = remove_ascii_control_characters(possible_gene)
        #  possible_gene = cleaned_possible_gene
        if isinstance(possible_gene, str):
            cleaned_possible_gene = remove_unicode_control_characters(possible_gene)
            possible_gene = cleaned_possible_gene

        possible_gene_split = possible_gene.split()
        # Skip out if the text is obviously single letter/spurious
        if len(possible_gene_split[0]) < 2:
            i += 1
            continue

        # Sometimes there is no space between Gene and NM_.  In these cases, try regular expression
        if len(possible_gene_split) == 1:
            NM_underscore = possible_gene_split[0].find('NM_')
            if NM_underscore > 1:
                possible_gene_split = [possible_gene_split[0][0:NM_underscore], possible_gene_split[0][NM_underscore:]]

        if possible_gene_split[0][0:3] == 'NM_':
            i += 1
            continue

        unp = None
        refseq = ''
        # Get the first white-space delimited substring, and remove any stray commas etc
        # Previous idea reg = re.compile("([0-9A-Z_]+)(orf)(0,1)([0-9A-Z_]+)[#@](0,1)")
        # Based on analysis of Jonathan's downloaded hgnc_complete_set.txt file...
        # Gene names are upper case, or numbers, underscore, dash with
        # possibly a lower case orf, then back to more numbers, upper case, and
        # possibly a terminating # or @ character
        gene_name_reg = re.compile("([0-9A-Z_-]+(orf)?[0-9A-Z_-]*[#@]?)$")

        mat = gene_name_reg.match(possible_gene_split[0])
        if mat:
            gene = mat.group(1)
        else:
            gene = possible_gene_split[0].translate(str.maketrans('', '', string.punctuation))

        # Skip out if our gene is common text that we immediately recognize as uninteresting
        if gene in ['Gene', 'Secondary', 'Heterozygous', 'Homozygous', 'Compound', 'None', 'Structural', 'DeNovo', 'De',
                    'Medically', 'Primary', 'Primary/Diagnostic', 'PrimaryDiagnostic', 'Table']:
            i += 1
            continue

        # LOGGER.info("Gene Candidate %s",gene)
        if not PDBMapProtein.ishgnc(gene):
            LOGGER.info("Row %3d: %-8s is not a recognized human gene" % (i, gene))
        else:
            effect = row[3]
            if isinstance(effect, str):
                cleaned_effect = remove_ascii_control_characters(effect)
                effect = cleaned_effect
            elif isinstance(effect, str):
                cleaned_effect = remove_unicode_control_characters(effect)
                effect = cleaned_effect
            else:
                effect = str(row[3]).strip().encode('utf-8')

            if gene not in genes:
                genes[gene] = {}

            for j in range(len(HeadersList)):
                if row[j] in gene_patterns:
                    rel = HeadersList[j]
                    status = gene_patterns[row[j]]
                    if rel not in genes[gene]:
                        genes[gene][rel] = [status]
                    elif status not in genes[gene][rel]:
                        genes[gene][rel].append(status)

            if not "missense" in effect.lower():
                LOGGER.info("Row %3d: %-8s skipped, non mis-sense mutation(%s)" % (i, gene, effect.replace('\n', '')))
                i += 2
                continue
            refseq_reg = re.compile("([XN]M_.*[0-9.]{2,30})")


            def refseq_match(_possible_refseq: str) -> str:
                _mat = refseq_reg.match(_possible_refseq.split()[0])
                return _mat.group(1) if _mat else None


            if len(possible_gene_split) > 1:
                refseq = None
                word_count = len(possible_gene_split)
                while not refseq and word_count > 1:
                    refseq = refseq_match(possible_gene_split[word_count - 1])
                    word_count -= 1
            else:  # Try to fish NM_ out of the cell immediately below Gene cell.  Sometimes the clinic sticks things there
                nextrow = df.iloc[i + 1]
                possible_refseq = str(nextrow[0]).strip().encode('utf-8')
                if possible_refseq and type(possible_refseq) == str and refseq_match(possible_refseq):
                    refseq = possible_refseq
                else:
                    nextrow = df.iloc[i + 2]
                    possible_refseq = str(nextrow[0]).strip().encode('utf-8')
                    if possible_refseq and type(possible_refseq) == str and refseq_match(possible_refseq):
                        refseq = possible_refseq

            if len(refseq) > 0:
                unpsForRefseq = PDBMapProtein.refseqNT2unp(refseq)
                if len(unpsForRefseq):
                    unp = ','.join(PDBMapProtein.refseqNT2unp(refseq))

            # if spreadsheet is missing the refseq, or refseq not found in idmapping file...
            # Then go with gene as last resort

            if unp == None:
                if refseq:
                    LOGGER.warning(
                        "Could not map refseq input [%s] uniprot ID (or it was missing).  Gene_refseq input in .xls file is: %s" % (
                            refseq, gene))
                else:
                    LOGGER.warning(
                        "Could not map refseq input to uniprot ID (or it was missing).  Gene_refseq input in .xls file is: %s" % gene)
                refseq = "RefSeqNotFound_UsingGeneOnly"
                unp = PDBMapProtein.hgnc2unp(gene)

            try:
                raw_mut = df.iloc[i + 2][2]
                mut = raw_mut.replace("p.", "")
            except IndexError:
                LOGGER.warning("Row %3d: %-8s  Could not read mutation info from spreadsheet" % (i, gene))
                i += 2
                continue

            # Convert 3-letter amino acid codes to 1-letter codes
            try:
                int(mut[1:-1])
            except (TypeError, IndexError, ValueError):
                # Convert mutation argument to 1-letter codes
                try:
                    mut = ''.join([longer_names[mut[:3].upper()], mut[3:-3], longer_names[mut[-3:].upper()]])
                except IndexError:
                    LOGGER.warning("Row %3d: %-8s  Failed to parse mutation %s" % (i, gene, raw_mut))
                    i += 2
                    continue  # Do not add this mutation

            try:
                chrom = df.iloc[i][1]     # Example 'chr1' in Spreadsheet column B
                pos = int(df.iloc[i+1][1]) # Example 150915463 below chr1 in column B
                change = df.iloc[i][2].strip()  # Example A->G to right of chr1 (not using c.809A>G to right of pos in column C)
                change = change[0] + '/' + change[-1] # Change format to A/G
              
            except (TypeError, IndexError):
                LOGGER.warning("Row %3d: Failed to parse Chr/Position/Change from columns B and C"%i)
       

            if genome != 'hg19':
                LOGGER.warning("GENOME IS SET TO %s for Row %d",genome,i)
            mutation_info = [gene,genome,chrom,pos,change,refseq,mut,unp]
            LOGGER.info("Row %3d: %-8s Adding %s" % (i, gene, " ".join([str(x) for x in mutation_info])))
            csv_rows.append(mutation_info)

            i += 1
    else:
        try:
            if row[0].find('Gene') != -1:
                GeneWordEncountered = True
                HeadersList = [str(x).strip().split()[0] for x in row]
        except TypeError:
            i += 1
            continue
    i += 1

with open(genes_json_filename, 'w') as fp:
    json.dump(genes, fp)
LOGGER.info("%d genes written to %s" % (len(genes), genes_json_filename))

with open(genes_txt_filename, 'w') as fp:
    for gene in genes:
        fp.write(gene + '\n')
LOGGER.info("%d genes written to %s" % (len(genes), genes_txt_filename))

if csv_rows:
    # Create a bed file of all our grch37 positions
    # And lift that over to GRCh38
    df = pd.DataFrame(csv_rows,columns=["gene","genome","chrom","pos","change","refseq","mutation","unp"])
    liftover_input_df = df[["chrom","pos"]]
    liftover_input_df.insert(1,'bed_zerobased_pos',liftover_input_df["pos"] - 1)
    liftover_input_df.insert(3,'uniquekey',liftover_input_df.index)
    liftover_bed_input = args.project + "_grch37_input.bed"
    liftover_output_lifted = args.project + "_grch38_lifted.bed"
    liftover_output_unlifted = args.project + "_grch38_unlifted.bed"

    liftover_input_df.to_csv(liftover_bed_input,sep='\t',header=False,index=False)
    import subprocess as sp
    liftover_arg_list = ['/dors/capra_lab/bin/liftOver',
         liftover_bed_input,
         '/dors/capra_lab/data/ucsc/liftOver/hg19ToHg38.over.chain.gz',
         liftover_output_lifted,
         liftover_output_unlifted]

    LOGGER.info("Running liftOver to convert to GRCh38 positions:\n%s",str(''.join(liftover_arg_list)))

    completed_process = sp.run(liftover_arg_list,
         capture_output=True,check=True,text=True)

    completed_process.check_returncode()

    try:
        df_lifted_grch38 = pd.read_csv(liftover_output_lifted,header=None,sep='\t',names=["chrom_hg38","bed_zerobased_pos_hg38","pos_hg38",'uniquekey'])
        LOGGER.info("%d genomic postions lifted to GRCh38",df_lifted_grch38.shape[0])
        df_lifted_grch38 = df_lifted_grch38.drop('bed_zerobased_pos_hg38',axis='columns').set_index('uniquekey')
    except pd.io.common.EmptyDataError:
        LOGGER.critical("Liftover to GRCh38 created null lifted output")

    try:
        df_unlifted = pd.read_csv(liftover_output_unlifted,header=None,sep='\t')
        LOGGER.critical("%d genomic postions could not be lifted to GRCh38",df_unlifted.shape[0])
    except pd.io.common.EmptyDataError:
        LOGGER.info("All positions lifted to GRCh38 successfully")

    if len(df_lifted_grch38): # Create a hg38 .vcf file that we can run vcf2missense.py on to create a new case
        df_vcf_hg38 = pd.DataFrame(
            columns = ["CHROM","POS","ID","REF","ALT"])
        df_vcf_hg38["CHROM"] = [chrom[3:] for chrom in df_lifted_grch38['chrom_hg38']]
        df_vcf_hg38["POS"] = [pos for pos in df_lifted_grch38['pos_hg38']]
        df_vcf_hg38["ID"] = ['.' for chrom in df_lifted_grch38['chrom_hg38']]
        df_vcf_hg38["REF"] = [df.loc[uniquekey]['change'][0] for uniquekey in df_lifted_grch38.index]
        df_vcf_hg38["ALT"] = [df.loc[uniquekey]['change'][-1] for uniquekey in df_lifted_grch38.index]
        hg38_vcffile = args.project+"_hg38.vcf"
        with open(hg38_vcffile,'w') as f:
            f.write("#")
            df_vcf_hg38.to_csv(f,sep='\t',index=False)

        LOGGER.info("Supplement hg38 liftover with VEP run")
        pdbmap_vep = PDBMapVEP(config_dict)

        vcf_reader =  pdbmap_vep.vcf_reader_from_file_supplemented_with_vep_outputs(hg38_vcffile)

        df_missense_grch38 = pd.DataFrame(columns=['gene','chrom','pos','transcript','unp','refseq','mutation'])

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
                df_missense_grch38 = df_missense_grch38.append(
                    {'gene': CSQ['SYMBOL'],
                     'chrom': vcf_record.CHROM,
                     'pos': vcf_record.POS,
                     'transcript': CSQ['Feature'],
                     'unp': unp,
                     'refseq': refseq,
                     'mutation': "%s%s%s"%(CSQ['Ref_AminoAcid'],CSQ['Protein_position'],CSQ['Alt_AminoAcid'])
                    },ignore_index=True)
        
        try:
            df_missense_grch38.to_csv(missense_csv_filename+".with_VEP_duplicates", header=True, encoding='ascii', sep=',')

        except Exception as e:
            LOGGER.critical("Unable to write %s: %s",missense_csv_filename+".with_VEP_duplicates", str(e))
            sys.exit(1)
    
        df_without_duplicates = df_missense_grch38.drop_duplicates(['gene','unp','refseq','mutation'])
        df = df_missense_grch38.set_index(['gene','unp','refseq','mutation'])
        for index,row in df_without_duplicates.iterrows():
            variant_index = (row['gene'],row['unp'],row['refseq'],row['mutation'])
            rows_with_various_transcripts = df.loc[variant_index]
            transcript_list = rows_with_various_transcripts['transcript'].tolist()
        # print("for variant_index %s rows are %s"%(str(variant_index),str(rows_with_various_transcripts)))
        # import pdb; pdb.set_trace()
        df_without_duplicates.at[index,'transcript'] = ';'.join(transcript_list)
    
        df_without_duplicates.to_csv(missense_csv_filename, header=True, encoding='ascii',sep=',')
        sys.exit("That's all for now")

    # Add the grch38 positions to the original missense case dataframe, as informational fields
    df = df.join(df_lifted_grch38,how='left',lsuffix='_l',rsuffix='_r')

    # Run the variant effect predictor on the GRCh19 positions to create a more robust transcript list
    pdbmap_vep = PDBMapVEP(config_dict)

    hg19_vcffile = args.project+"_hg19.vcf"
    hg19_vcf_df = df[['chrom','pos','change']]
    hg19_vcf_df['chrom'] = [chrom[3:] for chrom in hg19_vcf_df['chrom']]
    # pull apart the c.nnnnALT>REF format of the change field
    hg19_vcf_df.insert(3,'REF',[change[0] for change in hg19_vcf_df['change']])
    hg19_vcf_df.insert(4,'ALT',[change[-1] for change in hg19_vcf_df['change']])
    hg19_vcf_df = hg19_vcf_df.drop('change',axis='columns').rename(columns={'chrom': 'CHROM', 'pos': 'POS'})
    hg19_vcf_df.insert(2,'ID','.')
    with open(hg19_vcffile,'w') as f:
        f.write("#")
        hg19_vcf_df.to_csv(f,sep='\t',index=False)

    """ This might be a good idea below - EXCEPT that we are moving to GRCh38 - needs thought
    vcf_reader =  pdbmap_vep.vcf_reader_from_file_supplemented_with_vep_outputs(hg19_vcffile)

    unp = PDBMapProtein.enst2unp("ENST00000498193");

    hg19_missense_file = args.project + "_hg19_missense.csv"
    with open(hg19_missense_file,"w") as hg19_missense:
        hg19_missense.write(",gene,chrom,pos,transcript,unp,refseq,mutation\n")
        counter = 0
        for vcf_record in pdbmap_vep.yield_completed_vcf_records(vcf_reader):
            for CSQ in vcf_record.CSQ:
                unp = PDBMapProtein.enst2unp(CSQ['Feature'])
                if type(unp) is list: # << This is typical
                    unp = unp[0]
                refseq = "NA"
                if not unp:
                    LOGGER.warning("No Swiss-curated uniprot ID for transcript %s Gene %s"%(CSQ['Feature'],CSQ['SYMBOL']))
                    continue
                if unp:
                    refseq = PDBMapProtein.unp2refseqNT(unp)
                if type(refseq) is list and refseq:
                    refseq = refseq[0]
                else:
                    refseq = "NA"
                counter += 1
                hg19_missense.write("%d,%s,%s,%s,%s,%s,%s,%s%s%s\n"%(
                    counter,
                    CSQ['SYMBOL'],  # Gene
                    vcf_record.CHROM,
                    vcf_record.POS,
                    CSQ['Feature'], # ENST Transcript ID
                    unp,            # Isoform specific Uniprot ID
                    refseq,
                    CSQ['Ref_AminoAcid'],
                    CSQ['Protein_position'],
                    CSQ['Alt_AminoAcid']))
    """


    try:
        df.to_csv(missense_csv_filename, header=True, encoding='ascii')
    except Exception as e:
        LOGGER.critical("Unable to write %s: %s",missense_csv_filename, str(e))
        sys.exit(1)

    LOGGER.info("%d rows written to %s successfully", df.shape[0], missense_csv_filename)

else:
    print("No suitable mutations were found in the source excel file")

# Exit 0 (default)
