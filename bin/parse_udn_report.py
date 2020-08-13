#!/usr/bin/env python
"""\
This script will parse a UDN Patient Report (excel format) and 
save all of the mutations to a structured .csv file which is
formatted for input to add_udn_report.py\
"""

import logging
# logging.basicConfig(format='%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s',
#    datefmt='%d-%m-%Y:%H:%M:%S',)
logging.basicConfig(format='%(levelname)-8s [%(filename)s:%(lineno)-3d] %(message)s',)


import sys,csv,os

import re,string
import subprocess as sp
import pandas as pd
import argparse,configparser
import json
from lib.amino_acids import longer_names
from lib import PDBMapProtein,PDBMapModel,PDBMapSwiss
# Send INFO and higher logging information to stderr
# WARNING and higher log details are echoed to the log file
root  = logging.getLogger()
root.setLevel(logging.INFO)

logger = logging.getLogger()

import unicodedata
def remove_unicode_control_characters(s):
    return "".join(ch for ch in s if unicodedata.category(ch)[0]!="C")
def remove_ascii_control_characters(s):
    return ''.join(c for c in s if (ord(c) >= 32 and ord(c) < 128))
  
from psb_shared import psb_config
cmdline_parser = psb_config.create_default_argument_parser(__doc__,os.path.dirname(os.path.dirname(__file__)))
cmdline_parser.add_argument("project",type=str,help="Project ID (ex. UDN124356)",default=os.path.basename(os.getcwd()),nargs='?')

# parser.add_argument("udn_excel",type=str,help="Raw input UDN patient report (format: xls/xlsx)")
# parser.add_argument("udn_csv",type=str,help="Parsed output pipeline filename (e.g. filename.csv)")

args = cmdline_parser.parse_args()

required_items = ['output_rootdir','collaboration']

config,config_dict = psb_config.read_config_files(args,required_items)

udn_root_directory = os.path.join(config_dict['output_rootdir'],config_dict['collaboration'])
collaboration_dir = os.path.join(udn_root_directory,args.project)

udn_excel_filename = os.path.join(udn_root_directory,args.project,args.project + ".xlsx")  
missense_csv_filename = os.path.join(udn_root_directory,args.project,args.project + "_missense.csv")
genes_txt_filename = os.path.join(udn_root_directory,args.project,args.project + "_genes.txt")
genes_json_filename = os.path.join(udn_root_directory,args.project,args.project + "_genes.json")
logger.info('Parsing Excel file: %s'%udn_excel_filename)
logger.info(' Writing csv file %s '%missense_csv_filename)
logger.info(' Writing txt file %s '%genes_txt_filename)
logger.info('Writing json file %s '%genes_json_filename)

df = pd.read_excel(udn_excel_filename,header=0,iconverters={'Gene':str})
# print df

logger.info('Excel file read successfully with dimensions (pandas shape) of %s'%str(df.shape))
dfRows = df.shape[0]

genes = {}
gene_patterns={'\u25cf\u25cb':'Het',
   '\u25cb\u25cb':'WT',
   '\u25cf\u25cf':'Homo',
   '\u25cfy':'Het/Y',
   '\u25cby':'WT/Y',
   '\u25cf?':'Het/?'}

# Replaced with line below by ChrisMoth PDBMapProtein.load_idmapping("/dors/capra_lab/data/uniprot/idmapping/HUMAN_9606_idmapping_sprot.dat.gz")
# This really needs to come out of config file!!!!
idmappingFilename = config_dict['idmapping'] # "/dors/capra_lab/users/mothcw/mydata/idmapping/HUMAN_9606_idmapping_sprot.dat.gz"
logger.info("Loading uniprot idmapping from: %s"%idmappingFilename)
PDBMapProtein.load_idmapping(idmappingFilename)


GeneWordEncountered = False
if df.columns.values[0].strip() == 'Gene':
  logging.critical('The word "Gene" must not be in the top left header of the spreadsheet.  It must appear later, after the initial row which includes a case name.')
  sys.exit(1)

# The strings Gene, Chr, Change, Effect, Proband, Mother, Father 
# Will likely be the "headers" for the gene dictionary.
# However, the structure is a little flexible to allow for UDN slop
# See below where this important list is initialized
HeadersList = []

# This loop is not Pythonic - but the input .xls file is very unstructured.  We need to skip rows
i = 0
csv_rows = []
while i < dfRows:
  # import pdb; pdb.set_trace()
  row = df.iloc[i]
  # print i, row[0]
  if GeneWordEncountered: # Until we see "Gene" - we just keep looping
    try:
      possible_gene = row[0].strip() # .encode('utf-8')
    except:
      i += 1
      continue

    # if isinstance(possible_gene,str):
    #  cleaned_possible_gene = remove_ascii_control_characters(possible_gene)
    #  possible_gene = cleaned_possible_gene
    if isinstance(possible_gene,str):
      cleaned_possible_gene = remove_unicode_control_characters(possible_gene)
      possible_gene = cleaned_possible_gene

    possible_gene_split = possible_gene.split()
    # Skip out if the text is obviously single letter/spurious
    if len(possible_gene_split[0]) < 2:
      i += 1
      continue

    # Sometimes there is no space between Gene and NM_.  In these cases, try regular expression
    if len(possible_gene_split) == 1:
      NM_underscore = possible_gene_split[0].find('NM_');
      if (NM_underscore > 1):
        possible_gene_split = [possible_gene_split[0][0:NM_underscore],possible_gene_split[0][NM_underscore:]]

    if possible_gene_split[0][0:3] == 'NM_':
      i += 1
      continue
   
    # import pdb; pdb.set_trace() 
    unp = None
    refseq = ''
    # Get the first white-space delimited substring, and remove any stray commas etc
    # Previous idea reg = re.compile("([0-9A-Z_]+)(orf)(0,1)([0-9A-Z_]+)[#@](0,1)")
    # Based on analysis of Jonathan's downloaded hgnc_complete_set.txt file...
    # Gene names are upper case, or numbers, underscore, dash with 
    # possibly a lower case orf, then back to more numbers, upper case, and
    # possibly a terminating # or @ character
    gene_name_reg = re.compile("([0-9A-Z_\-]+(orf){0,1}[0-9A-Z_\-]*[#@]{0,1})$")

    mat = gene_name_reg.match(possible_gene_split[0])
    if mat:
      gene = mat.group(1)
    else:
      gene = possible_gene_split[0].translate(str.maketrans('','', string.punctuation))

    # Skip out if our gene is common text that we immediately recognize as uninteresting
    if gene in ['Gene','Secondary','Heterozygous','Homozygous','Compound','None','Structural','DeNovo', 'De', 'Medically', 'Primary', 'Primary/Diagnostic','PrimaryDiagnostic','Table']:
      i += 1
      continue

    # logger.info("Gene Candidate %s",gene)
    if not PDBMapProtein.ishgnc(gene):
      logger.info("Row %3d: %-8s is not a recognized human gene"%(i,gene))
    else: 
      effect = row[3];
      if isinstance(effect,str):
         cleaned_effect = remove_ascii_control_characters(effect)
         effect = cleaned_effect
      elif isinstance(effect,str):
         cleaned_effect = remove_unicode_control_characters(effect)
         effect = cleaned_effect
      else:
         effect = str(row[3]).strip().encode('utf-8')

      if not gene in genes:
        genes[gene] = {}

      for j in range(len(HeadersList)):
        if row[j] in gene_patterns:
          rel = HeadersList[j]
          status = gene_patterns[row[j]]
          if rel not in genes[gene]:
            genes[gene][rel] = [status]
          elif not status in genes[gene][rel]:
            genes[gene][rel].append(status)

      if not "missense" in effect.lower():
        logger.info("Row %3d: %-8s skipped, non mis-sense mutation(%s)"%(i,gene,effect.replace('\n','')))
        i += 2
        continue
      refseq_reg = re.compile("([XN]M_.*[0-9.]{2,30})")
      def refseq_match(possible_refseq: str) -> bool:
          mat = refseq_reg.match(possible_refseq.split()[0])
          return mat.group(1) if mat else None

      if len(possible_gene_split) > 1:
          refseq = None
          word_count = len(possible_gene_split)
          while not refseq and word_count > 1:
              refseq = refseq_match(possible_gene_split[word_count-1])
              word_count -= 1
      else: # Try to fish NM_ out of the cell immediately below Gene cell.  Sometimes the clinic sticks things there
        nextrow = df.iloc[i+1]
        possible_refseq = str(nextrow[0]).strip().encode('utf-8')
        if possible_refseq and type(possible_refseq) == str and refseq_match(possible_refseq):
          refseq = possible_refseq
        else:
          nextrow = df.iloc[i+2]
          possible_refseq = str(nextrow[0]).strip().encode('utf-8')
          if possible_refseq and type(possible_refseq) == str and refseq_match(possible_refseq):
            refseq = possible_refseq

      if len(refseq) > 0:
        unpsForRefseq = PDBMapProtein.refseqNT2unp(refseq)
        if len(unpsForRefseq):
          unp    = ','.join(PDBMapProtein.refseqNT2unp(refseq))

      #if spreadsheet is missing the refseq, or refseq not found in idmapping file...
      #Then go with gene as last resort

      if unp == None:
        if refseq:
          logger.warning("Could not map refseq input [%s] uniprot ID (or it was missing).  Gene_refseq input in .xls file is: %s"%(refseq,gene))
        else:
          logger.warning("Could not map refseq input to uniprot ID (or it was missing).  Gene_refseq input in .xls file is: %s"%gene)
        refseq = "RefSeqNotFound_UsingGeneOnly" 
        unp    = PDBMapProtein.hgnc2unp(gene)


      try:
        raw_mut = df.iloc[i+2][2]
        mut = raw_mut.replace("p.","")
      except:
        logger.warning("Row %3d: %-8s  Could not read mutation info from spreadsheet"%(i,gene))
        i += 2
        continue

      # Convert 3-letter amino acid codes to 1-letter codes
      try:
        int(mut[1:-1])
      except:
        # Convert mutation argument to 1-letter codes
        try:
          mut = ''.join([longer_names[mut[:3].upper()],mut[3:-3],longer_names[mut[-3:].upper()]])
        except:
          logger.warning("Row %3d: %-8s  Failed to parse mutation %s"%(i,gene,raw_mut))
          i += 2
          continue # Do not add this mutation

      mutation_info = [gene,refseq,mut,unp]
      logger.info("Row %3d: %-8s Adding %s"%(i,gene," ".join([str(x) for x in mutation_info])))
      csv_rows.append(mutation_info)
        
      i += 1
  else:
    try:
      if row[0].find('Gene') != -1:
        GeneWordEncountered = True
        HeadersList = [str(x).strip().split()[0] for x in row]
    except:
      i += 1
      continue
  i += 1

with open(genes_json_filename, 'w') as fp:
    json.dump(genes, fp)
logger.info("%d genes written to %s"%(len(genes),genes_json_filename))

with open(genes_txt_filename, 'w') as fp:
  for gene in genes:
    fp.write(gene+'\n')
logger.info("%d genes written to %s"%(len(genes),genes_txt_filename))

if csv_rows:    
  df = pd.DataFrame(csv_rows,columns=["gene","refseq","mutation","unp"])
  try:
  	 df.to_csv(missense_csv_filename,header=True,encoding='ascii')
  except Exception as e:
    logger.critical("Unable to write %s: %s"%(missense_csv_filename,str(e)))
    sys.exit(1)

  logger.info("%d rows written to %s successfully"%(df.shape[0],missense_csv_filename))
      
else:
   print("No suitable mutations were found in the source excel file")

# Exit 0 (default)      
