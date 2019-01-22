#!/usr/bin/env python2.7
#
# Project        : PSB Pipeline
# Filename       : undetermined.py
# Author         : Souhrid Muhkerjee
# Project        : PSB Pipeline
# Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : souhrid.mukherjee@vanderbilt.edu
# Date           : 2018-08-01
# Description    : Generate gene level dictionaries from UDN case mutation lists and phenotypes.
#=============================================================================#

"""\
Create gene level dictionaries containing Network Info. for an entire case
(ex. UDN123456 on command ine).

Configuration options must be provided with
   -c global.config file
   -u user.config overrides
   
"""

import logging,os,pwd,sys,grp,stat,string
import inspect # For status updates

capra_group = grp.getgrnam('capra_lab').gr_gid

def set_capra_group_sticky(dirname):
  try:
    os.chown(dirname, -1, capra_group)
  except:
    pass

  # Setting the sticky bit on directories also fantastic
  try:
    os.chmod(dirname, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IXGRP | stat.S_ISGID);
  except:
    pass

import time, datetime
import argparse,ConfigParser
import pprint

from logging.handlers import RotatingFileHandler
from logging import handlers
sh = logging.StreamHandler()
logger = logging.getLogger()
logger.addHandler(sh)

# Now that we've added streamHandler, basicConfig will not add another handler (important!)
formatter = logging.Formatter('%(asctime)s %(levelname)-4s [%(filename)20s:%(lineno)d] %(message)s',datefmt="%H:%M:%S")

logger.setLevel(logging.DEBUG)
sh.setLevel(logging.INFO)
sh.setFormatter(formatter)

import json
import sys
import os
import argparse
import cPickle as pickle
import networkx as nx
#import pygraphviz as pgv
import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as pd
import seaborn
from scipy.stats import hypergeom

# Setting up Arg Parser

default_global_config = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))) , "config" , "global.config")


cmdline_parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
cmdline_parser.add_argument("-c","--config",
help="PDBMap configuration profile for database access\n(default: %(default)s)", required=False,metavar="FILE",default=default_global_config)
cmdline_parser.add_argument("-u","--userconfig",
help="User specific settings and configuration profile overrides", required=True,metavar="FILE")
cmdline_parser.add_argument("-v","--verbose",
help="Include routine info log entries on stderr", action = "store_true")
cmdline_parser.add_argument("-d","--debug",
help="Include routine info AND 'debug' log entries on stderr", action = "store_true")

cmdline_parser.add_argument('-g', '--genes', type=str,
                    help='List of genes mutated in the patient (.txt or .json file) (optional: if not provided, will use gene list with same ID)',
                    required=False)
cmdline_parser.add_argument('-p', '--phenotypes', type=str,
                    help='List of patient phenotypes (.txt file) (optional: if not provided, will use gene list with same ID)',
                    required=False)
cmdline_parser.add_argument("--outdir",type=str,
                      help="Directory to place output and results")
cmdline_parser.add_argument("--uniquekey",type=str,default='',
                      help="A gene/refseq/mutation/structure/chain/flavor unique identifer for this particular job")

cmdline_parser.add_argument("project",type=str,help="The project ID, Example: UDN123456")
args,remaining_argv = cmdline_parser.parse_known_args()                    

args = cmdline_parser.parse_args()

fh = None # Lame but we can inspect this to see if our global logger is in place or not

def __initFileLogger(outdir):
  global fh
  assert not fh
  log_filename = os.path.join(outdir,'psb_genedicts.log')
  # Initiate logging in the output directory
  needRoll = os.path.isfile(log_filename)

  fh = RotatingFileHandler(log_filename, maxBytes=(1048576*5), backupCount=7)
  # formatter = logging.Formatter('%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s')
  # fh.setFormatter(log_formatter)
  fh.setLevel(logging.INFO)
  logger.addHandler(fh)

  if needRoll:
    fh.doRollover()
  logger.info("Log file is %s"%log_filename)

if args.debug:
  infoLogging = True
  sh.setLevel(logging.DEBUG)
elif args.verbose:
  infoLogging = True
  sh.setLevel(logger.info)

if args.outdir: # Great - start logging, else start below after config files read
  __initFileLogger(args.outdir)

# Load the structure details of all jobs that were not dropped

# print "Command Line Arguments"
logger.info("Command Line Arguments:\n%s"%pprint.pformat(vars(args)))


if not os.path.exists(args.config):
  msg = "Global config file not found: " + args.config
  logger.critical(msg)
  sys_exit_failure(msg)
  sys.exit(1)

config = ConfigParser.SafeConfigParser(allow_no_value=True)
configFilesRead = config.read([args.config,args.userconfig])
if len(configFilesRead) == 0:
  msg = "Unable to open config files: %s or %s  Exiting"%(args.config,args.userconfig)
  logger.critical(msg)
  sys_exit_failure(msg)
  sys.exit(1)

logger.info("Successful read of config files: %s"%str(configFilesRead))
config_dict = dict(config.items("Genome_PDB_Mapper")) # item() returns a list of (name, value) pairs

udn_root_directory = os.path.join(config_dict['output_rootdir'],config_dict['collaboration'])
collaboration_dir = os.path.join(udn_root_directory,args.project)

if not args.outdir: # We have not already done this - so do it now
  args.outdir = os.path.join(collaboration_dir,"casewide","Network_Analysis")

if not os.path.exists(args.outdir):
    logger.info("Creating args.outdir=%s"%args.outdir)
    os.makedirs(args.outdir)

set_capra_group_sticky(args.outdir)

if not fh:
  __initFileLogger(args.outdir)
  logger.info("The --outdir parameter is missing.  Set to %s"%args.outdir)



output_rootdir_shadow = '/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/UDN_Results/'

case_dir_shadow = os.path.join(output_rootdir_shadow, args.project)
pheno_dir = os.path.join(collaboration_dir, 'Phenotypes')

if not os.path.exists(collaboration_dir):
    msg = "{} directory has not been created.  Terminating".format(collaboration_dir)
    logger.critical(msg)
    sys_exit_failure(msg)
    sys.exit(1)


if not os.path.exists(case_dir_shadow):
    logger.info("Creating directory %s"%case_dir_shadow)
    os.makedirs(case_dir_shadow)
set_capra_group_sticky(case_dir_shadow)

statusdir =  os.path.abspath(os.path.join(args.outdir,"status"))
logger.info("Job status directory: %s"%statusdir)
if not os.path.exists(statusdir):
  try:
    os.makedirs(statusdir)
  except:
    pass

set_capra_group_sticky(statusdir)

for the_file in os.listdir(statusdir):
  file_path = os.path.join(statusdir, the_file)
  try:
    if os.path.isfile(file_path):
      os.unlink(file_path)
  except Exception as e:
    msg = "Unable to delete file %s in status directory"%file_path
    logger.execption(msg)
    sys.exit(1)

def __write_info_progress(info,progress):
  with open('%s/info.new'%statusdir,'w') as f:
    f.write(info + '\n')
  os.rename('%s/info.new'%statusdir,'%s/info'%statusdir)
  with open('%s/progress.new'%statusdir,'w') as f:
    f.write(progress)
  os.rename('%s/progress.new'%statusdir,'%s/progress'%statusdir)

def sys_exit_failure(info):
  __write_info_progress(info,"%s: %s\n"%(__file__,inspect.currentframe().f_back.f_lineno))
  # Mark this job as failed
  open("%s/FAILED"%statusdir,'w').close()
  logger.critical("Terminating in sys_exit_failure(): %s"%info)

  sys.exit(info)

def statusdir_info(info):
  __write_info_progress(info,"%s: %s\n"%(__file__,inspect.currentframe().f_back.f_lineno))


statusdir_info('Begun')

gene_list = []
# Loading genes
if not args.genes:
    #    genes_dict = json.load(open('{}/{case_name}_genes.json'.
    #                               format(onfig_dict['output_rootdir']gene_dir=gene_dir, case_name=name)))
    args.genes = os.path.join(collaboration_dir,"%s_genes.json"%args.project)
    logger.info("Genes filename derived from case ID and set to %s"%args.genes)

if args.genes.endswith('.json'):
    with open(args.genes) as json_genes:
        genes_dict = json.load(json_genes)
        logger.info('User supplied JSON file used to get {} genes: \'{}\''.format(len(genes_dict),args.genes))
elif args.genes.endswith('.txt'):
    with open(args.genes) as txt_genes:
        gene_list = [line for line in txt_genes.read().split('\n') if line]
    logger.info("Read {} genes from '{}'".format(len(gene_list),args.genes))
else:
    msg = "Invalid genes filename {}.  Filename must end with .txt or .json suffix".format(args.genes)
    logger.critical(msg)
    sys_exit_failure(msg)
    sys.exit(1)

#    logger.info('Gene list file derived from case ID \nID:{id}\nFile: \'{file}\''.format(
#        id=args.project,
#        file='{gene_dir}/{case_name}_genes.json'.
#            format(gene_dir=gene_dir, case_name=name)))

# Loading Phenotypes
if not args.phenotypes:
    args.phenotypes = os.path.join(collaboration_dir,"Phenotypes","%s_phenotypes.txt"%args.project)
    logger.info("Phenotypes list file derived from case ID and set to %s"%args.phenotypes)

with open(args.phenotypes) as txt_phenotypes:
    list_phenos_txt = [line for line in txt_phenotypes.read().split('\n') if line]
    logger.info("Read {} genes from '{}'".format(len(list_phenos_txt),args.phenotypes))

family_members_with_variant = {}
for gene in genes_dict:
    family_members_with_variant[gene] = set()     # Default is empty set of relatives with variant
    for family_member in genes_dict[gene]:
        if 'Het' in genes_dict[gene][family_member]:  # Match not only "Het" but also "Het/Y" and "Het/?" zygosities
            family_members_with_variant[gene].add(family_member)

variant_inheritance_info = {}
for gene in family_members_with_variant:
    if len(family_members_with_variant[gene]) > 1:  # Don't re-count "Proband"
        if set(['Mother', 'Father']).issubset(family_members_with_variant[gene]):
            variant_inheritance_info[gene] = ['Compound Het']
            if len(family_members_with_variant[gene]) > 3:
                variant_inheritance_info[gene].append('in sibling')

        elif 'Mother' in family_members_with_variant[gene]:
            # if i not in variant_inheritance_info:
            #     variant_inheritance_info[i] = []
            variant_inheritance_info[gene] = ['in Mom']

            if len(family_members_with_variant[gene]) > 2:
                variant_inheritance_info[gene].append('in sibling')

        elif 'Father' in family_members_with_variant[gene]:
            variant_inheritance_info[gene] = ['in Dad']
            if len(family_members_with_variant[gene]) > 2:
                variant_inheritance_info[gene].append('in sibling')

    # If no Mother/Father entry in the list of relatives with the variant, THEN
    # at future point, make inferences from brothers, half-sisters, etc.
    # For today, Jan 22 2019, we declare these "De Novo" when NEITHER parent
    # OR oarent sequencing is unavailable
    if not gene in variant_inheritance_info: # No record of (or in) Mother, Father
        variant_inheritance_info[gene] = ['De Novo']  # We only have Proband


positive_inheritance = []

for i in genes_dict:
    for j in genes_dict:
        if i != j:
            inheritance_union = list(set(variant_inheritance_info[i]).union(variant_inheritance_info[j]))
            inheritance_common = list(set(variant_inheritance_info[i]).intersection(variant_inheritance_info[j]))

            if inheritance_common == []:
                if (i, j) not in positive_inheritance and (j, i) not in positive_inheritance:
                    positive_inheritance.append((i, j))

            elif 'Compound Het' in inheritance_union or 'De Novo' in inheritance_union:
                if (i, j) not in positive_inheritance and (j, i) not in positive_inheritance:
                    positive_inheritance.append((i, j))

if gene_list == []:
    gene_list = sorted([str(i) for i in genes_dict])
# logging.debug('This message should go to the log file')
# logging.info
logger.info('ID:\'{}\''.format(args.project))
logger.info('Genes:{}'.format(gene_list))
logger.info('Inheritance patterns of variants:{}'.format(str(variant_inheritance_info)))
logger.info('Phenotypes:{}'.format(list_phenos_txt))
logger.info('Directory:\'{}\''.format(udn_root_directory))

logger.info('Importing Packages...')
logger.info('Importing Packages DONE!!')

logger.info('Importing Phenotype/Pathway Analysis modules...')
from Digenic_Classifier.HPO_for_gene import *
from Digenic_Classifier.Pathway_enrichment_analysis import *
from Digenic_Classifier.Regulation_analysis import *

logger.info('Importing Phenotype/Pathway Analysis modules  DONE!!')

logger.info('Importing Network Data...')
from Digenic_Classifier.load_UCSC_network_info import *

logger.info('Network Data improting..DONE!!')

logger.info('Importing Coexpression data...')
from Digenic_Classifier.coexpression_analysis import *

logger.info('Coexpression data Importing..DONE!!')

hpo_udn, hpo_udn_report = map_patient_phenotype_terms_to_HPO_codes(phenotype_list=list_phenos_txt)

# hpo_udn = open(sys.argv[3]).read().split('\n)[:-1]

logger.info("""
        Shorter_case_dictionary: 
            (case_name).json :
             [gene1][gene2]:{[Direct_Interaction]: True/False,
                             [Distance] : float(),
                             [Common_pathways] : True/False, 
                             [Coexpression] : float(),
                             [Common_phenotypes] : True/False, 
                             [Patient_Phenotype_Overlap] : 0-2,
                             [Inheritance]: True/False}
        Detailed_case_dictionary:
            (case_name)_detail.json:
             [gene1][gene2]:{[Direct_Interaction]: [details..],
                             [Distance] : [details..],
                             [Common_pathways] : [details..], 
                             [Coexpression] : [details..],
                             [Common_phenotypes] : [details..], 
                             [Patient_Phenotype_Overlap] : [details..],
                             [Inheritance]: [details..]}
             """)
dict_case = {}
dict_details = {}

for g in gene_list:
    gene = g
    ppi_con = []
    pwy_con = []
    txt_con = []

    if gene in G_ppi:
        ppi_con = [i for i in list(set(G_ppi.neighbors(gene)).intersection(gene_list)) if i != gene]
        for i in ppi_con:
            int_details = G_ppi.get_edge_data(gene, i)

            if gene not in dict_case:
                dict_case[gene] = {}
            if i not in dict_case[gene]:
                dict_case[gene][i] = {'Direct_Interaction': False, 'Distance': float('nan'),
                                      'Common_pathways': False, 'Coexpression': float('nan'),
                                      'Common_phenotypes': False, 'Patient_Phenotype_Overlap': 0, 'Inheritance': False}

            dict_case[gene][i]['Direct_Interaction'] = True
            logger.info('Direct PPI Interaction found between {gene1} and {gene2} ;\tDetails={d}'.
                         format(gene1=gene, gene2=i, d=str(int_details)))

            if gene not in dict_details:
                dict_details[gene] = {}
            if i not in dict_details[gene]:
                dict_details[gene][i] = {'Direct_Interaction': [], 'Distance': [],
                                         'Common_pathways': [], 'Coexpression': [],
                                         'Common_phenotypes': [], 'Patient_Phenotype_Overlap': [], 'Inheritance': []}

            dict_details[gene][i]['Direct_Interaction'].append(str(int_details))

            if (gene, i) in positive_inheritance or (i, gene) in positive_inheritance:
                dict_case[gene][i]['Inheritance'] = True
                inh_info = {gene: variant_inheritance_info[gene], i: variant_inheritance_info[i]}
                logger.info(
                    'Variants in {gene1} and {gene2} are not inherited from the same parent!\tDetails:{d}'.format(
                        gene1=gene, gene2=i, d=str(inh_info)))
            else:
                dict_case[gene][i]['Inheritance'] = False
                inh_info = {gene: variant_inheritance_info[gene], i: variant_inheritance_info[i]}
                logger.info(
                    'Variants in {gene1} and {gene2} do not fulfil Inheritance criteria\tDetails:{d}'.format(
                        gene1=gene, gene2=i, d=str(inh_info)))

            dict_details[gene][i]['Inheritance'] = [str(inh_info)]


    if gene in G_pwy:
        pwy_con = [i for i in list(set(nx.all_neighbors(G_pwy, gene)).intersection(gene_list)) if i != gene]
        for i in pwy_con:
            if G_pwy.has_edge(gene, i):
                int_details = G_pwy.get_edge_data(gene, i)
            else:
                int_details = G_pwy.get_edge_data(i, gene)

            if gene not in dict_case:
                dict_case[gene] = {}
            if i not in dict_case[gene]:
                dict_case[gene][i] = {'Direct_Interaction': False, 'Distance': float('nan'),
                                      'Common_pathways': False, 'Coexpression': float('nan'),
                                      'Common_phenotypes': False, 'Patient_Phenotype_Overlap': 0, 'Inheritance': False}

            dict_case[gene][i]['Direct_Interaction'] = True
            logger.info('Direct PWY Interaction found between {gene1} and {gene2} ;\tDetails={d}'.
                         format(gene1=gene, gene2=i, d=str(int_details)))

            if gene not in dict_details:
                dict_details[gene] = {}
            if i not in dict_details[gene]:
                dict_details[gene][i] = {'Direct_Interaction': [], 'Distance': [],
                                         'Common_pathways': [], 'Coexpression': [],
                                         'Common_phenotypes': [], 'Patient_Phenotype_Overlap': [], 'Inheritance': []}

            dict_details[gene][i]['Direct_Interaction'].append(str(int_details))

            if (i, gene) in positive_inheritance or (gene, i) in positive_inheritance:
                dict_case[gene][i]['Inheritance'] = True
                inh_info = {gene: variant_inheritance_info[gene], i: variant_inheritance_info[i]}
                logger.info(
                    'Variants in {gene1} and {gene2} are not inherited from the same parent!\tDetails:{d}'.format(
                        gene1=gene, gene2=i, d=str(inh_info)))
            else:
                dict_case[gene][i]['Inheritance'] = False
                inh_info = {gene: variant_inheritance_info[gene], i: variant_inheritance_info[i]}
                logger.info(
                    'Variants in {gene1} and {gene2} do not fulfil Inheritance criteria\tDetails:{d}'.format(
                        gene1=gene, gene2=i, d=str(inh_info)))

            dict_details[gene][i]['Inheritance'] = [str(inh_info)]

    if gene in G_txt:
        txt_con = [i for i in list(set(nx.all_neighbors(G_txt, gene)).intersection(gene_list)) if i != gene]
        for i in txt_con:
            if G_txt.has_edge(gene, i):
                int_details = G_txt.get_edge_data(gene, i)
            else:
                int_details = G_txt.get_edge_data(i, gene)

            if gene not in dict_case:
                dict_case[gene] = {}
            if i not in dict_case[gene]:
                dict_case[gene][i] = {'Direct_Interaction': False, 'Distance': float('nan'),
                                      'Common_pathways': False, 'Coexpression': float('nan'),
                                      'Common_phenotypes': False, 'Patient_Phenotype_Overlap': 0, 'Inheritance': False}

            dict_case[gene][i]['Direct_Interaction'] = True
            logger.info('Direct Txt Interaction found between {gene1} and {gene2} ;\tDetails={d}'.
                         format(gene1=gene, gene2=i, d=str(int_details)))

            if gene not in dict_details:
                dict_details[gene] = {}
            if i not in dict_details[gene]:
                dict_details[gene][i] = {'Direct_Interaction': [], 'Distance': [],
                                         'Common_pathways': [], 'Coexpression': [],
                                         'Common_phenotypes': [], 'Patient_Phenotype_Overlap': [], 'Inheritance': []}

            dict_details[gene][i]['Direct_Interaction'].append(str(int_details))

            if (i, gene) in positive_inheritance or (gene, i) in positive_inheritance:
                dict_case[gene][i]['Inheritance'] = True
                inh_info = {gene: variant_inheritance_info[gene], i: variant_inheritance_info[i]}
                logger.info(
                    'Variants in {gene1} and {gene2} are not inherited from the same parent!\tDetails:{d}'.format(
                        gene1=gene, gene2=i, d=str(inh_info)))
            else:
                dict_case[gene][i]['Inheritance'] = False
                inh_info = {gene: variant_inheritance_info[gene], i: variant_inheritance_info[i]}
                logger.info(
                    'Variants in {gene1} and {gene2} do not fulfil Inheritance criteria\tDetails:{d}'.format(
                        gene1=gene, gene2=i, d=str(inh_info)))

            dict_details[gene][i]['Inheritance'] = [str(inh_info)]

    ints = {}
    ppi_int = []
    pwy_int = []
    txt_int = []

    if gene in dists_ppi:
        ppi_int = [i for i in list(set(dists_ppi[gene]).intersection(gene_list)) if i != gene]
        for i in ppi_int:
            if i not in ints:
                ints[i] = [int(dists_ppi[gene][i])]
            else:
                ints[i] = list(set(ints[i]).union([int(dists_ppi[gene][i])]))

    if gene in dists_pwy:
        pwy_int = [i for i in list(set(dists_pwy[gene]).intersection(gene_list)) if i != gene]
        for i in pwy_int:
            if i not in ints:
                ints[i] = [int(dists_pwy[gene][i])]
            else:
                ints[i] = list(set(ints[i]).union([int(dists_pwy[gene][i])]))
    if gene in dists_txt:
        txt_int = [i for i in list(set(dists_txt[gene]).intersection(gene_list)) if i != gene]
        for i in txt_int:
            if i not in ints:
                ints[i] = [int(dists_txt[gene][i])]
            else:
                ints[i] = list(set(ints[i]).union([int(dists_txt[gene][i])]))
    paths = {}
    for i in ints:
        paths[i] = []
        if gene in G_ppi and i in G_ppi and nx.has_path(G_ppi, source=gene, target=i):
            dist = nx.shortest_path_length(G_ppi, source=gene, target=i)
            if dist == min(ints[i]):
                paths[i].append(('ppi', nx.shortest_path(G_ppi, source=gene, target=i)))
        if gene in G_pwy and i in G_pwy:
            if nx.has_path(G_pwy, source=gene, target=i):
                dist = nx.shortest_path_length(G_pwy, source=gene, target=i)
                if dist == min(ints[i]):
                    paths[i].append(('pwy', nx.shortest_path(G_pwy, source=gene, target=i)))
            if nx.has_path(G_pwy, source=i, target=gene):
                dist = nx.shortest_path_length(G_pwy, source=i, target=gene)
                if dist == min(ints[i]):
                    paths[i].append(('pwy', nx.shortest_path(G_pwy, source=i, target=gene)))
        if gene in G_txt and i in G_txt:
            if nx.has_path(G_txt, source=gene, target=i):
                dist = nx.shortest_path(G_txt, source=gene, target=i)
                if dist == min(ints[i]):
                    paths[i].append(('txt', nx.shortest_path(G_txt, source=gene, target=i)))
            if nx.has_path(G_txt, source=i, target=gene):
                dist = nx.shortest_path_length(G_txt, source=i, target=gene)
                if dist == min(ints[i]):
                    paths[i].append(('txt', nx.shortest_path(G_txt, source=i, target=gene)))

        if gene not in dict_case:
            dict_case[gene] = {}
        if i not in dict_case[gene]:
            dict_case[gene][i] = {'Direct_Interaction': False, 'Distance': float('nan'),
                                  'Common_pathways': False, 'Coexpression': float('nan'),
                                  'Common_phenotypes': False, 'Patient_Phenotype_Overlap': 0, 'Inheritance': False}

        dict_case[gene][i]['Distance'] = min(ints[i])
        logger.info('Min. distance between {gene1} and {gene2} is {value}; \tPath:{path}'.format(
            gene1=gene,
            gene2=i,
            value=min(ints[i]),
            path=paths[i]
        ))

        if gene not in dict_details:
            dict_details[gene] = {}
        if i not in dict_details[gene]:
            dict_details[gene][i] = {'Direct_Interaction': [], 'Distance': [],
                                     'Common_pathways': [], 'Coexpression': [],
                                     'Common_phenotypes': [], 'Patient_Phenotype_Overlap': [], 'Inheritance': []}

        dict_details[gene][i]['Distance'] = [str(paths[i])]

        if (i, gene) in positive_inheritance or (gene, i) in positive_inheritance:
            dict_case[gene][i]['Inheritance'] = True
            inh_info = {gene: variant_inheritance_info[gene], i: variant_inheritance_info[i]}
            logger.info(
                'Variants in {gene1} and {gene2} are not inherited from the same parent!\tDetails:{d}'.format(
                    gene1=gene, gene2=i, d=str(inh_info)))
        else:
            dict_case[gene][i]['Inheritance'] = False
            inh_info = {gene: variant_inheritance_info[gene], i: variant_inheritance_info[i]}
            logger.info(
                'Variants in {gene1} and {gene2} do not fulfil Inheritance criteria\tDetails:{d}'.format(
                    gene1=gene, gene2=i, d=str(inh_info)))

        dict_details[gene][i]['Inheritance'] = [str(inh_info)]

    genes_w_common_pathways = []
    common_paths = {}
    if gene in reactome_gene_to_path_codes:
        for i in reactome_gene_to_path_codes[gene]:
            common = [j for j in list(set(reactome_path_to_genes[i]).intersection(gene_list)) if j != gene]
            genes_w_common_pathways = list(set(genes_w_common_pathways).union(common))

            for j in common:
                if j not in common_paths:
                    common_paths[j] = [reactome_path_code_to_name[i]]
                else:
                    common_paths[j] = list(set(common_paths[j]).union([reactome_path_code_to_name[i]]))

    if gene in kegg_gene_to_path_codes:
        for i in kegg_gene_to_path_codes[gene]:
            common = [j for j in list(set(kegg_path_to_genes[i]).intersection(gene_list)) if j != gene]
            genes_w_common_pathways = list(set(genes_w_common_pathways).union(common))

            for j in common:
                if j not in common_paths:
                    common_paths[j] = [kegg_path_code_to_name[i]]
                else:
                    common_paths[j] = list(set(common_paths[j]).union([kegg_path_code_to_name[i]]))

    for i in genes_w_common_pathways:

        if gene not in dict_case:
            dict_case[gene] = {}
        if i not in dict_case[gene]:
            dict_case[gene][i] = {'Direct_Interaction': False, 'Distance': float('nan'),
                                  'Common_pathways': False, 'Coexpression': float('nan'),
                                  'Common_phenotypes': False, 'Patient_Phenotype_Overlap': 0, 'Inheritance': False}

        dict_case[gene][i]['Common_pathways'] = True
        logger.info('{gene1} and {gene2} are in common pathway(s);\tPathway(s):{pathway}'.format(
            gene1=gene,
            gene2=i,
            pathway=common_paths[i]
        ))

        if gene not in dict_details:
            dict_details[gene] = {}
        if i not in dict_details[gene]:
            dict_details[gene][i] = {'Direct_Interaction': [], 'Distance': [],
                                     'Common_pathways': [], 'Coexpression': [],
                                     'Common_phenotypes': [], 'Patient_Phenotype_Overlap': [], 'Inheritance': []}

        dict_details[gene][i]['Common_pathways'] = [str(common_paths[i])]

        if (i, gene) in positive_inheritance or (gene, i) in positive_inheritance:
            dict_case[gene][i]['Inheritance'] = True
            inh_info = {gene: variant_inheritance_info[gene], i: variant_inheritance_info[i]}
            logger.info(
                'Variants in {gene1} and {gene2} are not inherited from the same parent!\tDetails:{d}'.format(
                    gene1=gene, gene2=i, d=str(inh_info)))
        else:
            dict_case[gene][i]['Inheritance'] = False
            inh_info = {gene: variant_inheritance_info[gene], i: variant_inheritance_info[i]}
            logger.info(
                'Variants in {gene1} and {gene2} do not fulfil Inheritance criteria\tDetails:{d}'.format(
                    gene1=gene, gene2=i, d=str(inh_info)))

        dict_details[gene][i]['Inheritance'] = [str(inh_info)]

    coex_con = []
    if gene in coexpress_dict:
        coex_con = list(set([i for i in coexpress_dict[gene] if coexpress_dict[gene][i] < 50.]).intersection(gene_list))
        for i in coex_con:

            if gene not in dict_case:
                dict_case[gene] = {}
            if i not in dict_case[gene]:
                dict_case[gene][i] = {'Direct_Interaction': False, 'Distance': float('nan'),
                                      'Common_pathways': False, 'Coexpression': float('nan'),
                                      'Common_phenotypes': False, 'Patient_Phenotype_Overlap': 0, 'Inheritance': False}

            dict_case[gene][i]['Coexpression'] = coexpress_dict[gene][i]
            logger.info('{gene1} is coexpressed with {gene2} ;\tCoexpression Rank:{coex_rank}'.format(
                gene1=gene,
                gene2=i,
                coex_rank=coexpress_dict[gene][i]
            ))

            if gene not in dict_details:
                dict_details[gene] = {}
            if i not in dict_details[gene]:
                dict_details[gene][i] = {'Direct_Interaction': [], 'Distance': [],
                                         'Common_pathways': [], 'Coexpression': [],
                                         'Common_phenotypes': [], 'Patient_Phenotype_Overlap': [], 'Inheritance': []}

            dict_details[gene][i]['Coexpression'] = [str(coexpress_dict[gene][i])]

            if (i, gene) in positive_inheritance or (gene, i) in positive_inheritance:
                dict_case[gene][i]['Inheritance'] = True
                inh_info = {gene: variant_inheritance_info[gene], i: variant_inheritance_info[i]}
                logger.info(
                    'Variants in {gene1} and {gene2} are not inherited from the same parent!\tDetails:{d}'.format(
                        gene1=gene, gene2=i, d=str(inh_info)))
            else:
                dict_case[gene][i]['Inheritance'] = False
                inh_info = {gene: variant_inheritance_info[gene], i: variant_inheritance_info[i]}
                logger.info(
                    'Variants in {gene1} and {gene2} do not fulfil Inheritance criteria\tDetails:{d}'.format(
                        gene1=gene, gene2=i, d=str(inh_info)))

            dict_details[gene][i]['Inheritance'] = [str(inh_info)]

    genes_w_common_phenotypes = []
    common_phenos = {}
    if gene in hpo_gene_to_code:
        for j in hpo_gene_to_code[gene]:
            common = [i for i in list(set(hpo_code_to_gene[j]).intersection(gene_list)) if i != gene]
            genes_w_common_phenotypes = list(set(genes_w_common_phenotypes).union(common))

            for i in common:
                if i not in common_phenos:
                    common_phenos[i] = [hpo_code_to_name[j]]
                else:
                    common_phenos[i] = list(set(common_phenos[i]).union([hpo_code_to_name[j]]))

    for i in genes_w_common_phenotypes:

        if gene not in dict_case:
            dict_case[gene] = {}
        if i not in dict_case[gene]:
            dict_case[gene][i] = {'Direct_Interaction': False, 'Distance': float('nan'),
                                  'Common_pathways': False, 'Coexpression': float('nan'),
                                  'Common_phenotypes': False, 'Patient_Phenotype_Overlap': 0, 'Inheritance': False}

        dict_case[gene][i]['Common_phenotypes'] = True
        logger.info('{gene1} and {gene2} share phenotypes ;\tCommon phenotypes:{phenos}'.format(
            gene1=gene,
            gene2=i,
            phenos=common_phenos[i]
        ))

        if gene not in dict_details:
            dict_details[gene] = {}
        if i not in dict_details[gene]:
            dict_details[gene][i] = {'Direct_Interaction': [], 'Distance': [],
                                     'Common_pathways': [], 'Coexpression': [],
                                     'Common_phenotypes': [], 'Patient_Phenotype_Overlap': [], 'Inheritance': []}

        dict_details[gene][i]['Common_phenotypes'] = [str(common_phenos[i])]

        if (i, gene) in positive_inheritance or (gene, i) in positive_inheritance:
            dict_case[gene][i]['Inheritance'] = True
            inh_info = {gene: variant_inheritance_info[gene], i: variant_inheritance_info[i]}
            logger.info(
                'Variants in {gene1} and {gene2} are not inherited from the same parent!\tDetails:{d}'.format(
                    gene1=gene, gene2=i, d=str(inh_info)))
        else:
            dict_case[gene][i]['Inheritance'] = False
            inh_info = {gene: variant_inheritance_info[gene], i: variant_inheritance_info[i]}
            logger.info(
                'Variants in {gene1} and {gene2} do not fulfil Inheritance criteria\tDetails:{d}'.format(
                    gene1=gene, gene2=i, d=str(inh_info)))

        dict_details[gene][i]['Inheritance'] = [str(inh_info)]

        pheno_overlap1 = list(set(hpo_gene_to_code[gene]).intersection(hpo_udn))
        pheno_overlap2 = list(set(hpo_gene_to_code[i]).intersection(hpo_udn))

        if len(pheno_overlap1) > 0 or len(pheno_overlap2) > 0:
            if len(pheno_overlap1) > 0 and len(pheno_overlap2) > 0:

                dict_case[gene][i]['Patient_Phenotype_Overlap'] = 2
                logger.info('{gene1} and {gene2} connected with patient phenotypes ;\tPatient phenotypes overlap for \
                {gene1}:{phenos1}\tPatient phenotypes overlap for {gene2}:{phenos2}'.format(
                    gene1=gene,
                    gene2=i,
                    phenos1=[hpo_code_to_name[j] for j in pheno_overlap1],
                    phenos2=[hpo_code_to_name[j] for j in pheno_overlap2]
                ))
            else:

                dict_case[gene][i]['Patient_Phenotype_Overlap'] = 1
                if len(pheno_overlap1) > 0:
                    logger.info(
                        '{gene1} connected with patient phenotypes ;\tPatient phenotypes overlap:{phenos1}'.format(
                            gene1=gene,
                            phenos1=[hpo_code_to_name[j] for j in pheno_overlap1]
                        ))
                else:
                    logger.info(
                        '{gene2} connected with patient phenotypes ;\tPatient phenotypes:{phenos2}'.format(
                            gene2=i,
                            phenos2=[hpo_code_to_name[j] for j in pheno_overlap2]
                        ))
        else:

            dict_case[gene][i]['Patient_Phenotype_Overlap'] = 0

# write final directories to case and shadow directories
for dict_fname_base in [os.path.join(args.outdir,args.project), os.path.join(case_dir_shadow,args.project)]:
    dict_case_fname = dict_fname_base + ".json"
    
    with open(dict_case_fname,'w') as fp:
        json.dump(dict_case, fp)
        logger.info("dict_case with len=%d added to %s"%(len(dict_case),dict_case_fname))
    
    dict_details_fname = dict_fname_base + "_detail.json"
    with open(dict_details_fname,'w') as fp:
        json.dump(dict_details, fp)
        logger.info("dict_details with len=%d added to %s"%(len(dict_details),dict_details_fname))


# Mark this program as complete
statusdir_info('Success.  dict_case len=%d  dict_details len=%d'%(len(dict_case),len(dict_details)))
open("%s/complete"%statusdir,'w').close()
