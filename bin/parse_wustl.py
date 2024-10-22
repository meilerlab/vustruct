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
# =============================================================================#
"""Parse a simple comma-delimted file in format:
Gene,Variant
Gene1,A1234B
Gene2,C4567D"""

# See main check for cmd line parsing
import argparse, configparser
import traceback
import sys
import os
import gzip
import re
import pandas as pd
from lib import PDBMapProtein
from lib import PDBMapSQLdb
from lib import PDBMapTranscriptUniprot
import logging
from psb_shared import psb_config

sh = logging.StreamHandler()
LOGGER = logging.getLogger()
LOGGER.addHandler(sh)

log_format_string = '%(asctime)s %(levelname)-4s [%(filename)16s:%(lineno)d] %(message)s'
date_format_string = '%H:%M:%S'
log_formatter = logging.Formatter(log_format_string, date_format_string)

LOGGER.setLevel(logging.DEBUG)
sh.setLevel(logging.INFO)
sh.setFormatter(log_formatter)

cmdline_parser = psb_config.create_default_argument_parser(__doc__, os.path.dirname(os.path.dirname(__file__)))
cmdline_parser.add_argument("wustlfile",
                            type=str,
                            metavar="FILE",
                            help="WashU spreadsheet input",
                            default=os.path.basename(os.getcwd()) + ".xlsx",
                            nargs='?')

args = cmdline_parser.parse_args()

required_config_items = ['idmapping', 'refseq_genomic_gff']
config, config_dict = psb_config.read_config_files(args, required_config_items)
PDBMapProtein.load_idmapping(config_dict['idmapping'])
PDBMapSQLdb.set_access_dictionary(config_dict)
PDBMapSQLdb.open_global_connection(config_dict)

df = pd.read_excel(args.wustlfile, header=0)
df = pd.read_excel(args.wustlfile, header=0, converters={k: str for k in range(df.shape[1])})
print("%d rows read from %s" % (len(df), args.wustlfile))
print("Columns = %s" % df.columns)

columns_from_original = ['gene', 'refseq', 'exon', 'c', 'prot']

# This will ultimately provide all the columns for the vcf output at the end
# Some additional columns, lowe rcase, are kept around to supply as INFO
genomic_df = pd.DataFrame(
    columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'exon', 'cpos', 'refseq', 'wustl_c',
             'wustl_prot', 'xref_unp'])
final_df = pd.DataFrame(columns=columns_from_original + ['unp'])
genomic_simple_change_regex = re.compile("^c.[ACTG]([0-9]+)[ACTG]$")  # Could be a missense
protein_variant_regex = re.compile("^p.[A-Z]([0-9]+)[A-Z]$")

complmentary_nucleic_acid = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A'
}

for index, row in df.iterrows():
    all_text = str(row['AAChange.refGene'])
    split_text = re.split('[,;]', all_text)
    for substring in split_text:
        colon_delimited_strings = substring.split(':')
        if len(colon_delimited_strings):
            if len(colon_delimited_strings) == 5:
                as_dictionary = dict(zip(columns_from_original, colon_delimited_strings))

                genomic_simple_change_match = genomic_simple_change_regex.match(as_dictionary['c'])
                if not genomic_simple_change_match:
                    print("Genomic change %s cannot be a missense variant" % as_dictionary['c'])
                    continue

                variant_match = protein_variant_regex.match(as_dictionary['prot'])

                if not variant_match:
                    print("Protein string %s not a missense variant" % as_dictionary['prot'])
                    continue

                genomic_df_row = {
                    'CHROM': None,  # We shall fill this in later
                    'POS': None,  # We shall fill this in later
                    'ID': as_dictionary['refseq'] + '_' + as_dictionary['c'],
                    # WUSTL giving us RNA nucleotides in their c.XnnnY format
                    'REF': as_dictionary['c'][2],
                    'ALT': as_dictionary['c'][-1],
                    'QUAL': '.',  # NULL value for VCF FILE
                    'FILTER': '.',  # NULL value for VCF FILE
                    'INFO': '.',  # NULL value for VCF FILE
                    'exon': int(as_dictionary['exon'][4:]),
                    'cpos': int(genomic_simple_change_match.group(1)),
                    'refseq': as_dictionary['refseq'],
                    'wustl_c': as_dictionary['c'],
                    'wustl_prot': as_dictionary['prot'],
                    'xref_unp': None
                }

                variant_resnum = int(variant_match.group(1))

                print("Protein string %s is a missense variant %d" % (as_dictionary['prot'], variant_resnum))

                unp = None
                if 'refseq' in as_dictionary and as_dictionary['refseq']:
                    unps = PDBMapProtein.refseqNT2unp(as_dictionary['refseq'])
                    if unps:
                        if not type(unps) is list:
                            unps = [unps]
                        for unp in unps:
                            transcript = PDBMapTranscriptUniprot(unp)
                            if len(transcript.aa_seq) < variant_resnum:
                                print("REJECTION %s because variant # %d is not in transcript" % (
                                    str(as_dictionary),
                                    variant_resnum))
                                continue
                            if as_dictionary['prot'][2] != transcript.aa_seq[variant_resnum - 1]:
                                print("REJECTION %s because amino acid is wrong - must be %s not %s" % (
                                    str(as_dictionary),
                                    transcript.aa_seq[variant_resnum - 1],
                                    as_dictionary['prot'][2]))
                                continue
                            genomic_df_row['xref_unp'] = unp
                            as_dictionary['unp'] = unp
                            break

                genomic_df = pd.concat([genomic_df, pd.DataFrame.from_dict([genomic_df_row])],
                                       ignore_index=True)
                final_df = pd.concat([final_df, pd.DataFrame.from_dict([as_dictionary])],
                                     ignore_index=True)
            else:
                print("Skipping substring which lacks 5 elements: %s %s" % (
                    len(colon_delimited_strings), colon_delimited_strings))

genomic_df_filename = os.path.splitext(args.wustlfile)[0] + "_exons.csv"
genomic_df.to_csv(genomic_df_filename)
LOGGER.info("%d rows written to %s", len(genomic_df), genomic_df_filename)

# Gather up the NC_ (chromosome) ids and exonic information for each reseq ID, also chromosomes

refseq_chroms = {}

# Initiialize empty dictionary of exon locations for each refseq ID.
# We'll populate these as we work through all of this
refseq_exons = {}
for index, genomic_row in genomic_df.iterrows():
    assert genomic_row['refseq'].find('.') == -1, \
        "Halting as the wustl file has a versioned refseq %s" % genomic_row['refseq']
    refseq_exons[genomic_row['refseq']] = {}

# For each refseq ID, also mine a list of (start, end, strand) CDS entries
# These are simply in order of encounter as we drill through the file
refseq_cds_entries = {}

LOGGER.info("Scanning huge refseq 'gff' file: %s", config_dict['refseq_genomic_gff'])
# Parse through the huge gff file from refseq and parse out what we need
with gzip.open(config_dict['refseq_genomic_gff'], 'rt') as gff:
    for line in gff:
        tabsplit = line.split('\t')
        if tabsplit[0] == '#':
            continue
        if len(tabsplit) < 9:
            continue
        (seqname, source, feature, start, end, score, strand, frame, attribute) = tabsplit
        seqname_prefix = seqname[0:3]
        if seqname_prefix == 'NC_':
            if feature == 'exon':
                attribute_split = attribute.split(';')
                # You will have a string like
                # ID=exon-XM_017000207.3-16; and you need to get the refseq id and the exon
                # number (16) out of it...
                if (source == 'Gnomon' and attribute_split[0].startswith('ID=exon-XM')) or \
                        (source == 'BestRefSeq' and attribute_split[0].startswith('ID=exon-NM')):
                    try:
                        _, refseq, exon = attribute_split[0].split('-')
                    except ValueError:
                        LOGGER.warning("Can't unpack attribute_split %s", attribute_split)
                    # If this refseq  not in our WUSTL format file, don't bother with it.
                    # Noe that the refseq file has version suffixes.  Our WUST file does not
                    refseq_without_version = refseq.split('.')[0]
                    if refseq_without_version in refseq_exons:
                        refseq_exons[refseq_without_version][int(exon)] = {
                            'chrom': refseq_chroms[seqname],
                            'start': int(start),
                            'end': int(end),
                            'strand': strand
                        }
            elif feature == 'CDS':
                # The marking of the coding regions happens in order we encouter them in the file
                attribute_split = attribute.split(';')
                if (source == 'Gnomon' and attribute_split[0].startswith('ID=cds-XP')) or \
                        (source == 'BestRefSeq' and attribute_split[0].startswith('ID=cds-NP')):
                    # We have to iterate to find the Parent-rna-entry to get our refseq
                    refseq = None
                    for sub_attribute in attribute_split[1:]:
                        if sub_attribute.startswith("Parent=rna-"):
                            refseq = sub_attribute[11:]
                            break
                    assert refseq is not None, "Could not find refseq ID in line: %s" % line

                    refseq_without_version = refseq.split('.')[0]

                    # Is our refseq ID one that we initially (from input file) care about?
                    # If so then store it!
                    if refseq_without_version in refseq_exons:
                        _cds_list_for_refseq = refseq_cds_entries.get(refseq_without_version, [])
                        _cds_list_for_refseq.append({
                            'chrom': refseq_chroms[seqname],
                            'start': int(start),
                            'end': int(end),
                            'strand': strand
                        }
                        )
                        refseq_cds_entries[refseq_without_version] = _cds_list_for_refseq

            elif feature == 'region':
                if source == 'RefSeq':
                    # Line Example: NC_000016.10	RefSeq	region	.. ID=....;chromosome=15
                    attribute_split = attribute.split(';')
                    chromosome_found = False
                    for attribute_fragment in attribute_split:
                        if attribute_fragment.startswith('chromosome='):
                            # Don't forget that chromosome can be X, Y
                            refseq_chroms[seqname] = attribute_fragment.split('=')[1]
                            chromosome_found = True
                            break
                    if not chromosome_found:  # Perhaps it is mitocondrial
                        for attribute_fragment in attribute_split:
                            if attribute_fragment.startswith('genome=mitochondrion'):
                                # Don't forget that chromosome can be X, Y
                                refseq_chroms[seqname] = "MT"
                                chromosome_found = True
                                break

# Now fill the chrom, pos, etc positions n the genomic dataframe
# Make sure index here is 0..N-1 rows
genomic_df.reset_index()
for i in range(len(genomic_df)):
    genomic_df_row = genomic_df.iloc[i]
    refseq = genomic_df_row['refseq']
    if refseq not in refseq_cds_entries:
        LOGGER.warning("CDS entries for Refseq ID %s were not found in %s", refseq, config_dict['refseq_genomic_gff'])
        continue

    # Cpos is offset in all the coding DNA - not just one transcript
    # It is the nnn from the original c.XmmmX string we parsed in originally
    cpos = genomic_df_row['cpos']

    for cds_entry in refseq_cds_entries[refseq]:
        cds_length = 1 + int(cds_entry['end']) - int(cds_entry['start'])
        if cpos <= cds_length:  # Great - we've found our exon
            genomic_df.at[i, 'CHROM'] = cds_entry['chrom']
            if cds_entry['strand'] == '+':
                genomic_df.at[i, 'POS'] = cds_entry['start'] + (cpos - 1)
            elif cds_entry['strand'] == '-':  # < Thanks David Rinker
                genomic_df.at[i, 'POS'] = cds_entry['end'] - (cpos - 1)
                # It appears that when we have a "backwards" transcript, then we need to flip
                # The bases too
                genomic_df.at[i, 'ALT'] = complmentary_nucleic_acid[genomic_df.iloc[i]['ALT']]
                genomic_df.at[i, 'REF'] = complmentary_nucleic_acid[genomic_df.iloc[i]['REF']]
            else:
                assert False, "cds_entry %s must be + or - for refseq %s" % (str(cds_entry), refseq)

            cpos = 0
            break
        else:
            cpos -= cds_length

    # *** THIS CODE BLOCK ON HOLD BECAUSE USING THE EXONS is not what we need.  We need the CPOS
    # exon = 1
    # matched = False
    # while exon in refseq_exons[refseq]:
    #    exoninfo = refseq_exons[refseq][exon]
    #    exon_length = 1+int(exoninfo['end']) -int(exoninfo['start'])
    #    # assert exon_length % 3 == 0, "%d-%d=%d"%(int(exoninfo['end']),int(exoninfo['start']),int(exoninfo['end']) -int(exoninfo['start']))
    #    if cpos < exon_length: # Great - we've found our exon
    #        genomic_df.at[i, 'CHROM'] = exoninfo['chrom']
    #        genomic_df.at[i, 'POS'] = int(exoninfo['start']) + cpos
    #        matched=True
    #        break
    #    cpos -= exon_length
    #    exon += 1
    # assert matched

    # if int(genomic_df_row['exon']) not in refseq_exons[refseq]:
    #     LOGGER.warning("Refseq ID %s exon %d was not found in %s", refseq, int(genomic_df_row['exon']),config_dict['refseq_genomic_gff'])
    #     continue
    exoninfo = refseq_exons[refseq][int(genomic_df_row['exon'])]

genomic_df.dropna(subset='CHROM', inplace=True)

final_vcf_filename = os.path.splitext(args.wustlfile)[0] + ".vcf"
with open(final_vcf_filename, 'w') as f_vcf:
    f_vcf.write("##fileformat=VCFv4.2\n")
    f_vcf.write("##source=%s  processed through %s\n" % (args.wustlfile, __file__))
    f_vcf.write('##INFO=<ID=RS,Number=1,Type=String,Description="Refseq ID from input">\n')
    f_vcf.write('##INFO=<ID=EX,Number=1,Type=Integer,Description="Exon from input">\n')
    genomic_df['INFO'] = genomic_df.apply(
        lambda row: ("RS=%s;EX=%s" % (row['refseq'], row['exon'])), axis=1)

    genomic_vcf_df = genomic_df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']]

    f_vcf.write('#')  # To left of CHROM on header line we need a #
    genomic_vcf_df.to_csv(f_vcf, index=False, sep='\t')
print("%d rows written to %s" % (len(genomic_vcf_df), final_vcf_filename))

sys.exit(0)

final_df_filename = sys.argv[1] + '.csv'
final_df.to_csv(final_df_filename)

assert len(gene_df) > 1, "No data rows were found in %s" % args.genefile
assert 'gene' in gene_df.columns, "A header row in %s is required, and must contain 'gene' and 'variant'"
assert 'variant' in gene_df.columns, "A header row in %s is required, and must contain 'gene' and 'variant'"

missense_df = pd.DataFrame(columns=['gene', 'unp', 'refseq', 'mutation'])

for index, row in gene_df.iterrows():
    gene_id = PDBMapProtein.gene2gene_id(row['gene'])
    if not gene_id:
        fail_message = "Gene %s in your input file lacks a corresponding uniprot identifier in the idmapping file" % \
                       row['gene']
        LOGGER.critical(fail_message)
        sys.exit(fail_message)
    canonical_unp = PDBMapProtein.hgnc2unp(row['gene'])
    LOGGER.info("Gene %s maps to gene_id: %s and canonical unp: %s", row['gene'], gene_id, canonical_unp)
    canonical_unp_transcript = PDBMapTranscriptUniprot(canonical_unp)
    LOGGER.info("Transcript for %s has len %d", canonical_unp_transcript.id, len(canonical_unp_transcript.aa_seq))

    try:
        ref_aa = row['variant'][0]
        alt_aa = row['variant'][-1]
        variant_pos = int(row['variant'][1:-1])
    except:
        sys.exit(
            "Variant amino acid format must be A1234B.  You have [%s] in row index %d" % (row['variant'], variant_pos))

    assert ref_aa == canonical_unp_transcript.aa_seq[
        variant_pos - 1], "You provided variant %s but position %d is %s not %s" % (
        row['variant'], variant_pos, canonical_unp_transcript.aa_seq[variant_pos - 1], ref_aa)

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
         }, ignore_index=True)

missense_filename = os.path.basename(os.getcwd()) + "_missense.csv"
LOGGER.info("Writing %d rows of missense variants to %s", len(missense_df), missense_filename)
missense_df.to_csv(missense_filename)
