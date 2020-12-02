#!/usr/bin/env python
#
# Project        : PathProx
# Filename       : pathprox.py  
#                  (pathprox2.py reflects new logging, configuration, output structure of pipelineV2 architecture)
#                  (pathprox3.py processes complex and ties into various PDBMap* updates for improved isoform treatment)
# Authors        : R. Michael Sivley
#                : Xiaoyi Dou
#                : Chris Moth
# Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2016-11-04 to 2019-09
# Description    : Predicts the pathogenicity of missense variants by their
#                : spatial position relative to known pathogenic and neutral
#                : variation in protein structure.
#=============================================================================#
"""Pathprox predicts (scores) the pathogenicity of a missense variant by its
spatial position relative to known pathogenic and neutral variation in its 
3D protein structure."""

## Package Dependenecies ##
import sys
import os
import shutil
import logging
import gzip
import lzma
import csv
import platform
import grp
import stat
import json
import re
import string
import subprocess as sp
from time import strftime
from itertools import combinations
from collections import OrderedDict
import pprint
from subprocess import Popen, PIPE

from typing import Dict

## cache_dir = '/tmp/sqlcache'  # Must be overridden early


import configparser

# Numerical
import pandas as pd
import numpy as np
import math
TOL = 1e-5 # zero-tolerancection with a given name return the same logger instance. This means that logger instances never need to be passed between d

# Stats
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.stats import norm, percentileofscore, mannwhitneyu
from scipy.spatial import KDTree
# Import our PDBMap library
from lib import PDBMapGlobals
from lib import PDBMapAlignment
from lib import PDBMapProtein
from lib import PDBMapSQLdb
from lib import PDBMapTranscriptBase
from lib import PDBMapTranscriptUniprot
from lib import PDBMapTranscriptFasta
from lib import PDBMapTranscriptEnsembl
from lib import PDBMapGenomeVariants
from lib import PDBMapGnomad
from lib import PDB36Parser
from lib.PDBMapAlignment import sifts_best_unps

# BioPython routines to deal with pdb files
from Bio import PDB 
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio import *

if __name__ == "__main__":
    # Logging setup can vary a lot depending on whether pathprox is used as import module or mainline route
    ch = logging.StreamHandler()
    LOGGER = logging.getLogger()
    LOGGER.addHandler(ch)

    log_format_string = '%(asctime)s %(levelname)-4s [%(filename)16s:%(lineno)d] %(message)s'
    date_format_string = '%H:%M:%S'
    log_formatter = logging.Formatter(log_format_string,date_format_string)
    ch.setFormatter(log_formatter)

    LOGGER.setLevel(logging.DEBUG)
    ch.setLevel(logging.INFO)

    from psb_shared import psb_config,psb_progress,psb_capra_group


    # Setup the Command Line Argument Parser
    cmdline_parser = psb_config.create_default_argument_parser(__doc__,os.path.dirname(os.path.dirname(__file__)),add_help=True)


    cmdline_parser.add_argument("-o","--outdir",type=str,
         help="Directory to use for output and results")
    cmdline_parser.add_argument("--uniquekey",type=str,required=False,
         help="A gene/refseq/mutation/structure/chain/flavor unique identifer for this pathprox job")
    cmdline_parser.add_argument("--variants",type=str,
                      help="Comma-separated list of amino acids in {transcript|chain}:HGVS or \
                            a filename containing one identifier per line")
    # Input parameters
    group = cmdline_parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--pdb',type=str,
         help="4 character PDB ID with optional .chain suffix")
    group.add_argument('--biounit',type=str,
         help="4 character PDB ID with optional .chain suffix")
    group.add_argument('--modbase',type=str,
         help="Modbase 20 model ID with optional .chain suffix")
    group.add_argument('--swiss',type=str,
         help="Swissmodel ID with optional .chain suffix")
    group.add_argument('--usermodel',type=str,metavar='FILE',
         help="Filename of a user model.  Requires explicit transcript specifications")
    #cmdline_parser.add_argument("entity",type=str,
    #                   help="Gene ID, UniProt AC, PDB ID, or PDB filename")
    args,args_remaining = cmdline_parser.parse_known_args()

    if args.debug:
        LOGGER.setLevel(logging.DEBUG)
    if not args.outdir:
        if args.uniquekey:
            args.outdir = args.uniquekey
        else:
            args.outdir = 'results'
        timestamp = strftime("%Y-%m-%d")
        if not args.no_timestamp and timestamp not in args.outdir:
            args.outdir += "_%s"%timestamp
        LOGGER.info("The --outdir parameter is missing.  Set to %s"%args.outdir)

    psb_capra_group.makedirs_capra_lab(args.outdir,'Main outdir creation')
       
    # cache_dir = "/tmp/cache"
    log_filename = os.path.join(args.outdir,args.uniquekey + ".log")

    LOGGER.info("Log file is %s"%log_filename)
    needRoll = os.path.isfile(log_filename)

    # Write current version of this script to output directory
    # pw_name = pwd.getpwuid( os.getuid() ).pw_name # example jsheehaj or mothcw
    from logging.handlers import RotatingFileHandler
    from logging import handlers

    fh = RotatingFileHandler(log_filename, maxBytes=(1048576*5), backupCount=7)
    # formatter = logging.Formatter('%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s')
    fh.setFormatter(log_formatter)
    fh.setLevel(logging.INFO)
    LOGGER.addHandler(fh)

    if needRoll:
      fh.doRollover()
else:
    LOGGER = logging.getLogger(__name__)

import inspect # Update status with stack trace info if we crash

# Warnings
from warnings import filterwarnings,resetwarnings

# Removed October -> because built-in zscore chokes on sigma=0 from scipy.stats.mstats import zscore

# Plotting
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("white")

# Machine Learning
filterwarnings('ignore',category=RuntimeWarning)
from sklearn.linear_model import LogisticRegression as LR

from  sklearn.preprocessing import scale
import  sklearn.metrics 
# from  sklearn.metrics import roc_curve 
# from sklearn.metrics import precision_recall_curve
# from sklearn.metrics import average_precision_score
resetwarnings()

from Bio.SubsMat.MatrixInfo import blosum100

def normalize_reflect_and_symmetrize_blosum100_matrix():
    # I do not linke that we are modifying a resource from external library
    blosum100_values_as_list = list(blosum100.values())
    min_blosum100 = np.min(blosum100_values_as_list)
    max_blosum100 = np.max(blosum100_values_as_list)
    # The blosum100 keys are tuples of amino acid single letter codes
    for aa_pair,blosum_score in list(blosum100.items()):
        # Transform to 0..1, then invert so that severe mutations get high scores
        blosum100[aa_pair] = 1.0-((float(blosum_score)-min_blosum100)/(max_blosum100-min_blosum100))
        # Store redundant copy with reversed amino acids, for easier queries
        blosum100[tuple(reversed(aa_pair))] = blosum100[aa_pair] 


normalize_reflect_and_symmetrize_blosum100_matrix()
LOGGER.info("Normalzation, reflection, and symmetrization of blosum100 complete");

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# PDBMap
from lib import PDBMapSwiss
from lib import PDBMapModbase2020
# from lib import PDBMapModbase2013
# from lib import PDBMapModbase2016
from lib import PDBMapProtein
# from lib.amino_acids import longer_names


#=============================================================================#
## Function Definitions ##
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

def makedirs_capra_lab(DEST_PATH, module_name):
    if not os.path.exists(DEST_PATH):
        os.makedirs(DEST_PATH)
    try:
        assert os.path.exists(DEST_PATH)
    except:
       sys_exit_failure("Fatal: Function %s failed to create destination path %s" % (module_name,DEST_PATH))
    set_capra_group_sticky(DEST_PATH)

def load_structure(id,coord_filename):
    tryParser36 = False
    with gzip.open(coord_filename,'rt') if coord_filename.split('.')[-1]=="gz" else \
         lzma.open(coord_filename,'rt') if coord_filename.split('.')[-1]=='xz' else \
             open(coord_filename,'r') as fin:
        filterwarnings('ignore',category=PDBConstructionWarning)
        try:
            structure = PDBParser().get_structure(id,fin)
        except ValueError:
            tryParser36 = True # We will make a last ditch effort to read this because of alpha in int columns
            structure = None

        resetwarnings()

    if tryParser36:
        LOGGER.critical("Attempting Hybrid36 PDB Parser for %s"%coord_filename)
        with gzip.open(coord_filename,'rt') if coord_filename.split('.')[-1]=="gz" else \
             lzma.open(coord_filename,'rt') if coord_filename.split('.')[-1]=='xz' else \
                 open(coord_filename,'r') as fin:
            filterwarnings('ignore',category=PDBConstructionWarning)
            structure = PDB36Parser().get_structure(id,fin)
            resetwarnings()

    return structure


def PDBMapComplex_load_pdb_align_chains(pdb_id,try_biounit_first,chain_to_transcript):
    """From a pdb ID, load the structure, and sifts-align its chains to the best uniprot IDs
       (as returned by sifts"""
    chain_to_alignment = {}
    is_biounit = 0
    structure = None
    if try_biounit_first:
        coord_filename = os.path.join(config_dict['pdb_dir'],"biounit","mmCIF","divided",pdb_id.lower()[1:3],"%s-assembly1.cif.gz"%pdb_id.lower())
        if os.path.exists(coord_filename): # Great we will load the biounit from mmCIF file
            try:
                mmCIF_parser = MMCIFParser(QUIET=True)
                LOGGER.debug("Loading %s",coord_filename)
                with gzip.open(coord_filename,'rt') as structure_fin:
                    structure =  mmCIF_parser.get_structure(pdb_id.lower(),structure_fin)
                is_biounit = 1
            except:
                LOGGER.info("The biounit file %s was not found.  Attempting non-biounit pdb"%os.path.basename(coord_filename))
        else:  # Typically, we load a biounit from PDB file - but good to have tried CIF first
            coord_filename = os.path.join(config_dict['pdb_dir'],"biounit","PDB","divided",pdb_id.lower()[1:3],"%s.pdb1.gz"%pdb_id.lower())
            try:
                structure = load_structure(pdb_id,coord_filename)
                is_biounit = 1
            except:
                LOGGER.info("The biounit file %s was not found.  Attempting normal pdb"%os.path.basename(coord_filename))
  
    if not structure:  # Always the case if try_biounit_first is false
        coord_filename = os.path.join(config_dict['pdb_dir'],"structures","divided","mmCIF",pdb_id.lower()[1:3],"%s.cif.gz"%pdb_id.lower())
        try:
            mmCIF_parser = MMCIFParser(QUIET=True)
            LOGGER.debug("Loading %s",coord_filename)
            with gzip.open(coord_filename,'rt') as structure_fin:
                structure =  mmCIF_parser.get_structure(pdb_id.lower(),structure_fin)
        except:
            LOGGER.info("CIF structure file %s was not found or not loaded."%os.path.basename(coord_filename))

    if not structure:
        coord_filename = os.path.join(config_dict['pdb_dir'],"structures","divided","pdb",pdb_id.lower()[1:3],"pdb%s.ent.gz"%pdb_id.lower())
        LOGGER.warning("Reverting to .pdb format: %s",coord_Filename)
        structure = load_structure(pdb_id,coord_filename)
        LOGGER.info("Success loading %s",coord_filename)

    sifts_chain_to_best_unp = sifts_best_unps(structure)
    for chain_letter in sifts_chain_to_best_unp:
        if chain_letter not in structure[0]:
            continue  # It is A-OK to have a biounit which lacks chains from sifts
        if chain_letter in chain_to_transcript:
            best_unp_transcript = chain_to_transcript[chain_letter] # Use the transcript inited from command line
        else:
            # Before we get too excited about an isoform specific ID
            # If there is no cross-ref'd ENST transcript from idmapping, then revert to canonical
            ensemble_transcript_ids = PDBMapProtein.unp2enst(sifts_chain_to_best_unp[chain_letter])
            if not ensemble_transcript_ids:
                LOGGER.critical("Unp %s is best for chain %s. HOWEVER it lacks ENST*.. so reverting to canonical",sifts_chain_to_best_unp[chain_letter],chain_letter)
                canonical_unp = PDBMapProtein.best_unp(sifts_chain_to_best_unp[chain_letter].split('-')[0])
                # If this unp is NOT in our curated/reviewed set from uniprot, we need to skip this chain
                if PDBMapProtein.unp2uniparc(canonical_unp) is None:
                    LOGGER.warning("Chain %s has canonical unp %s, but it has not been reviewed by uniprot.  Skipping"%(chain_letter,canonical_unp))
                    continue
                sifts_chain_to_best_unp[chain_letter] = canonical_unp
                LOGGER.critical("For chain %s, unp now=%s",chain_letter,sifts_chain_to_best_unp[chain_letter])

            best_unp_transcript = PDBMapTranscriptUniprot(sifts_chain_to_best_unp[chain_letter])

        alignment = PDBMapAlignment()
        (success,message) = alignment.align_trivial(best_unp_transcript,structure,chain_id=chain_letter)
        if success:
            LOGGER.info("%s to chain %s Trivial align successful for residues in [%d,%d]\n%s"%(best_unp_transcript.id,chain_letter,next(iter(alignment.seq_to_resid)),next(reversed(alignment.seq_to_resid)),alignment.aln_str))
        else:
            LOGGER.info("Unable to trivially align transcript %s to %s.%s"%(best_unp_transcript.id,structure.id,chain_letter))
            is_canonical = PDBMapProtein.isCanonicalByUniparc(best_unp_transcript.id)
            if is_canonical:
                (success,message) = alignment.align_sifts_canonical(best_unp_transcript,structure,chain_id=chain_letter)
            if not success: 
                mmCIF_parser = MMCIFParser(QUIET=True)
                mmcif_structure_filename = os.path.join(config_dict['pdb_dir'],'structures','divided','mmCIF',pdb_id.lower()[1:3],"%s.cif.gz"%pdb_id)
                if os.path.exists(mmcif_structure_filename):
                    LOGGER.debug("Loading %s",mmcif_structure_filename)
                    if mmcif_structure_filename.endswith(".gz"):
                        structure_fin = gzip.open(mmcif_structure_filename,'rt')
                    else:
                        structure_fin = open(mmcif_structure_filename,'rt')
                    temp_structure =  mmCIF_parser.get_structure(pdb_id.lower(),structure_fin)
                    mmCIF_dict = mmCIF_parser._mmcif_dict # << This is a little dangerous but the get_structure code above looks quite clean
                    pdb_seq_resid_xref = PDBMapAlignment.create_pdb_seq_resid_xref(mmCIF_dict)

                    (success,message) = alignment.align_sifts_isoform_specific(best_unp_transcript,temp_structure,pdb_seq_resid_xref,chain_id=chain_letter,is_canonical=is_canonical)
                if not success:
                    LOGGER.warning("Unable to align with sifts: %s",message)
                    LOGGER.warning("Attempting to align chain %s with biopython call",chain_letter)
                    success,message = alignment.align_biopython(best_unp_transcript,structure,chain_letter)
                    if not success:
                        LOGGER.critical("Also Unable to align with biopython: %s",message)


        assert success,message
        chain_to_alignment[chain_letter] = alignment
        chain_to_transcript[chain_letter] = best_unp_transcript
    return structure,is_biounit,chain_to_transcript,chain_to_alignment
        


def structure_lookup(io,sid,bio=True,chain=None):
  """ Returns coordinate files for a PDB ID """
  q    = "SELECT DISTINCT biounit FROM Chain "
  q   += "WHERE label=%s AND structid=%s AND biounit>0 "
  if chain:
    q += "AND chain=%s"
    res  = [r[0] for r in io.secure_cached_query(cache_dir,q,(io.slabel,sid.upper(),chain),
                                          cursorclass='Cursor')]
  else:
    res  = [r[0] for r in io.secure_cached_query(cache_dir,q,(io.slabel,sid.upper()),
                                          cursorclass='Cursor')]
  if bio and res: # Biological assemblies were requested and found
    LOGGER.info("Using the first biological assembly for %s."%sid)
    flist = []
    loc   = "%s/biounits/%s.pdb%d.gz"
    for b in res:
      f = loc%(config_dict['pdb_dir'],sid.lower(),b)
      if not os.path.exists(f):
        msg  = "Coordinate file missing for %s[%s]\n"%(sid,0)
        msg += "Expected: %s\n"%f
        raise Exception(msg)
      flist.append((sid,b,f))
    return flist
  else:   # No biological assemblies were found; Using asymmetric unit
    if bio:
      LOGGER.info("Using the asymmetric unit for %s."%sid)
    loc = "%s/structures/pdb%s.ent.gz"
    f   = loc%(config_dict['pdb_dir'],sid.lower())
    if not os.path.exists(f):
      msg  = "Coordinate file missing for %s[%s]\n"%(sid,0)
      msg += "Expected: %s\n"%f
      raise Exception(msg)
    return [(sid,0,f)]

#def model_lookup(io,mid):
#  """ Returns coordinate files for a ModBase ID """
#  PDBMapModel.load_modbase(config_dict['modbase2016_dir'],config_dict['modbase2016_summary'])
#  PDBMapModel.load_modbase(config_dict['modbase2013_dir'],config_dict['modbase2013_summary'])
#  f = PDBMapModel.get_coord_file(mid.upper())
#  LOGGER.info("File location for %s: %s"%(mid.upper(),f))
#  #f = "%s/Homo_sapiens_2016/model/%s.pdb.gz"%(args.modbase_dir,mid.upper())
#  if not os.path.exists(f):
#    try:
#      cmd = ["xz","-d",'.'.join(f.split('.')[:-1])+'.xz']
#      sp.check_call(cmd)
#    except:
#      msg  = "Coordinate file missing for %s\n"%mid
#      msg += "Expected: %s\n"%f
#      raise Exception(msg)
#    cmd = ["gzip",'.'.join(f.split('.')[:-1])]
#    sp.check_call(cmd)
#  return [(mid,0,f)]

def swiss_lookup(io,model_id):
  """ Returns coordinate files for a Swiss ID """
  PDBMapSwiss.load_swiss_INDEX_JSON(config_dict['swiss_dir'],config_dict['swiss_summary']);
  f = PDBMapSwiss.get_coord_file(model_id)
  LOGGER.info("File location for %s: %s"%(model_id,f))
  if not os.path.exists(f):
    msg  = "Coordinate file missing for %s\n"%model_id
    msg += "Expected: %s\n"%f
    logging.getLogger(__name__).critical(msg)
    raise Exception(msg)
  return [(model_id,0,f)]

#def uniprot_lookup(io,ac):
#  """ Returns coordinate files for a UniProt AC """
#  entities = io.load_unp(ac)
#  if not entities:
#    msg = "No structures or models associated with %s\n"%ac
#    LOGGER.warning(msg)
#    return None
#  flist = []
#  for etype,ac in entities:
#    if etype == "structure":
#      flist.extend(structure_lookup(io,ac))
#    elif etype == "model":
#      flist.extend(model_lookup(io,ac))
#    elif etype == "swiss":
#      flist.extend(swiss_lookup(io,ac))
#  return flist

def structure_renumber_per_alignments(structure):
    from Bio.PDB.Structure import Structure
    from Bio.PDB.Model import Model
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Residue import Residue
    renumbered_structure = Structure(structure.id)
    for model in structure:
        renumbered_model = Model(model.id)
        for chain in model:
            if chain.id not in chain_to_alignment:
                # This is ugly - but we are just pointing at the old chain
                # It solves the problem of having disordered residues that are not brought over by .copy()
                renumbered_model.add(chain)
            else:
                alignment = chain_to_alignment[chain.id]
                renumbered_chain = Chain(chain.id)
                for residue in chain:
                    if residue.id in alignment.resid_to_seq:
                        renumbered_residue = residue.copy()
                        renumbered_residue.id = (' ',alignment.resid_to_seq[residue.id],' ')
                        # The above copy routines fail to copy disorder properly - so just wipe out all notion of disorder
                        for atom in renumbered_residue:
                            atom.disordered_flag = 0
                        renumbered_residue.disordered = 0
                        renumbered_chain.add(renumbered_residue)
                renumbered_model.add(renumbered_chain)

        renumbered_structure.add(renumbered_model)
    return renumbered_structure


def default_var_query():
  """ Default string for variant queries """
  s  = "SELECT distinct protein_pos,ref_amino_acid,alt_amino_acid,chain FROM GenomicConsequence b "
  s += "INNER JOIN GenomicIntersection a ON a.gc_id=b.gc_id "
  w  = "WHERE a.slabel in ('pdb','modbase','swiss') and b.label=%s AND consequence LIKE '%%missense_variant%%' "
  w += "AND a.structid=%s and length(ref_amino_acid)=1 and length(alt_amino_acid)=1 "
  return s,w

def sequence_var_query():
  """ Alternate string for custom protein models """
  s  = "SELECT distinct protein_pos,ref_amino_acid,alt_amino_acid FROM GenomicConsequence b "
  w  = "WHERE label=%s and consequence LIKE '%%missense_variant%%' and length(ref_amino_acid)=1 and length(alt_amino_acid)=1 "
  if args.isoform:
    LOGGER.info("Querying GenomicConsequence for transcript=%s",args.isoform)
    w += "AND transcript='%s'"%args.isoform
  else:
    LOGGER.warning("No --isoform  Querying GenomicConsequence for non-specific uniprot ID=%s",args.unp)
    w += "AND uniprot=%s"%args.unp

    if not args.fasta:
      msg  = "\nWARNING: Reference isoform was not specified with --fasta. \n"
      msg  = "       : Are there multiple isoforms for this protein?\n"
      msg += "       : Explictly declare isoform to avoid mis-alignments\n\n"
      LOGGER.warning(msg)
  return s,w

def query_1kg(io,sid,refid=None,chains=None,indb=False):
  """ Query natural variants (1000G) from PDBMap """
  if indb:
    s,w = default_var_query()
    if chains:
      w += "AND chain in (%s) "%','.join(["'%s'"%c for c in chains])
    f   = ("1kg3",sid)
    c   = ["unp_pos","ref","alt","chain"]
  else:
    s,w = sequence_var_query()
    f   = ("1kg3",refid)
    c   = ["unp_pos","ref","alt"]
  q   = s+w
  res = [list(r) for r in io.secure_cached_query(cache_dir,q,f,cursorclass="Cursor")]
  # If user-specified model...
  if not indb:
    # Add chains IDs
    res = [r+[c] for r in res for c in chains]
    # # Add null pdb_pos (not yet aligned)
    # res = [[None]+r for r in res]
  return res

def query_exac(io,sid,unp=None,chains=None,indb=False):
  """ Query natural variants (ExAC) from PDBMap """
  if indb:
    LOGGER.info("Known structure, querying pre-intersected variants...")
    s,w = default_var_query()
    if chains:
      w += "AND chain in (%s) "%','.join(["'%s'"%c for c in chains])
    f   = ("exac",sid)
    # c   = ["unp_pos","ref","alt","chain"]
  else:
    LOGGER.info("Unknown structure, querying all variants for %s..."%unp)
    s,w = sequence_var_query()
    f   = ["exac"]
    # c   = ["unp_pos","ref","alt"]
  q   = s+w
  res = [list(r) for r in io.secure_cached_query(cache_dir,q,f,cursorclass="Cursor")]
  # If user-specified model...
  if not indb:
    # Add chains IDs
    res = [r+[c] for r in res for c in chains]
    # # Add null pdb_pos (not yet aligned)
    # res = [[None]+r for r in res]
  return res

def query_gnomad(io,sid,refid=None,chains=None,indb=False):
  """ Query natural variants (gnomAD) from PDBMap """
  """if indb:
    s,w = default_var_query()
    if chains:
      w += "AND chain in (%s) "%','.join(["'%s'"%c for c in chains])
    f   = ("gnomad",sid)
    c   = ["unp_pos","ref","alt","chain"]
  else:"""
  s,w = sequence_var_query()
  f   = ["gnomad"]
  c   = ["unp_pos","ref","alt"]
  q   = s+w
  res = [list(r) for r in io.secure_cached_query(cache_dir,q,f,cursorclass="Cursor")]
  # If user-specified model...
  # if not indb:
  # Add chains IDs
  res = [r+[c] for r in res for c in chains]
  # # Add null pdb_pos (not yet aligned)
  # res = [[None]+r for r in res]
  return res

def query_benign(io,sid,refid=None,chains=None,indb=False):
  """ Query benign variants (ClinVar) from PDBMap """
  if indb:
    s,w = default_var_query()
    if chains:
      w += "AND chain in (%s) "%','.join(["'%s'"%c for c in chains])
    f   = ("clinvar",sid)
    c   = ["unp_pos","ref","alt","chain"]
  else:
    s,w = sequence_var_query()
    f   = ("clinvar",refid)
    c   = ["unp_pos","ref","alt"]
  s  += "INNER JOIN clinvar d "
  s  += "ON b.chr=d.chr and b.start=d.start "
  b   = "AND pathogenic=0"
  q   = s+w+b
  res = [list(r) for r in io.secure_cached_query(cache_dir,q,f,cursorclass="Cursor")]
  # If user-specified model...
  if not indb:
    # Add chains IDs
    res = [r+[c] for r in res for c in chains]
    # # Add null pdb_pos (not yet aligned)
    # res = [[None]+r for r in res]
  return res

def query_pathogenic(io,sid,refid=None,chains=None,indb=False):
  """ Query pathogenic variants (ClinVar) from PDBMap """
  if indb:
    s,w = default_var_query()
    if chains:
      w += "AND chain in (%s) "%','.join(["'%s'"%c for c in chains])
    f   = ("clinvar",sid)
    c   = ["unp_pos","ref","alt","chain"]
  else:
    s,w = sequence_var_query()
    f   = ["clinvar"]
    c   = ["unp_pos","ref","alt"]
  s  += "INNER JOIN clinvar d "
  s  += "ON b.chr=d.chr and b.start=d.start "
  p   = "AND pathogenic=1"
  q   = s+w+p
  res = [list(r) for r in io.secure_cached_query(cache_dir,q,f,cursorclass="Cursor")]
  # If user-specified model...
  if not indb:
    # Add chains IDs
    res = [r+[c] for r in res for c in chains]
    # # Add null pdb_pos (not yet aligned)
    # res = [[None]+r for r in res]
  return res

def query_drug(io,sid,refid=None,chains=None,indb=False):
  """ Query drug response-affecting variants from PDBMap """
  if indb:
    s,w = default_var_query()
    if chains:
      w += "AND chain in (%s) "%','.join(["'%s'"%c for c in chains])
    f   = ("clinvar",sid)
    c   = ["unp_pos","ref","alt","chain"]
  else:
    s,w = sequence_var_query()
    f   = ("clinvar",refid)
    c   = ["unp_pos","ref","alt"]
  s  += "INNER JOIN clinvar d "
  s  += "ON b.chr=d.chr and b.start=d.start "
  d   = "AND d.clnsig like '%7%'"
  q   = s+w+d
  res = [list(r) for r in io.secure_cached_query(cache_dir,q,f,cursorclass="Cursor")]
  # If user-specified model...
  if not indb:
    # Add chains IDs
    res = [r+[c] for r in res for c in chains]
    # # Add null pdb_pos (not yet aligned)
    # res = [[None]+r for r in res]
  return res

def query_cosmic_obsolete_code(io,sid,refid=None,chains=None,indb=False):
    pass
    ##""" Query somatic variants (COSMIC) from PDBMap """
    ##if indb:
    ##  s,w = default_var_query()
    ##  if chains:
    ##    w += "AND chain in (%s) "%','.join(["'%s'"%c for c in chains])
    ##  f   = ("cosmic",sid)
    ##  c   = ["unp_pos","ref","alt","chain"]
    ##else:
    ##  s,w = sequence_var_query()
    ##  f   = ["cosmic"]
    ##  c   = ["unp_pos","ref","alt"]
    ##m   = "INNER JOIN cosmic d "
    ##m  += "ON b.chr=d.chr and b.start=d.start "
    ##m  += "AND d.cnt>1 "
    ##q   = s+m+w
    ##res = [list(r) for r in io.secure_cached_query(cache_dir,q,f,cursorclass="Cursor")]
    ### If user-specified model...
    ##if not indb:
    ##  # Add chains IDs
    ##  res = [r+[c] for r in res for c in chains]
    ##  # # Add null pdb_pos (not yet aligned)
    ##  # res = [[None]+r for r in res]
    ##return res

def query_tcga(io,sid,refid=None,chains=None,indb=False):
  """ Query somatic variants (TCGA) from PDBMap """
  """ 2019-Feb 11- Get Away from the default query mode which requires pre-intersection """
  """if indb:
    s,w = default_var_query()
    if chains:
      w += "AND chain in (%s) "%','.join(["'%s'"%c for c in chains])
    f   = ("tcga",sid)
    c   = ["unp_pos","ref","alt","chain"]
  else:"""
  s,w = sequence_var_query()
  f   = ["tcga"]
  c   = ["unp_pos","ref","alt"]
  m   = "INNER JOIN tcga d "
  m  += "ON b.chr=d.chr and b.start=d.start "
  # m  += "AND d.cnt>1 "  For tcga, the count is actually 1 for all of them 
  q   = s+m+w
  res = [list(r) for r in io.secure_cached_query(cache_dir,q,f,cursorclass="Cursor")]
  # If user-specified model...
  # if not indb:  ACTUALLY - since we are doing the sequence_var_query, add back the chains...
  # Add chains IDs
  res = [r+[c] for r in res for c in chains]
  return res

# def get_coord_files(entity,io):
#   """ Returns the relevant coordinate file or None if no file found """
#   # Check if the entity is a filename. If so, assume PDB/ENT format
#   global unp_flag
#   if os.path.isfile(entity):
#     exts = 1 + int(os.path.basename(entity).split('.')[-1]=='gz')
#     sid  = '.'.join(os.path.basename(entity).split('.')[:-exts]) 
#     return [(sid,-1,entity)]
#   else:
#     LOGGER.info("Attempting to identify entity: %s..."%entity)
#     etype = io.detect_entity_type(entity)
#     if   etype == "structure":
#       LOGGER.info("Input is a PDB ID. Locating coordinate file...")
#       return structure_lookup(io,entity,args.chain is None)
#     elif etype == "model":
#       LOGGER.info("Input is a ModBase model ID. Locating coordinate file...")
#       return model_lookup(io,entity)
#     elif etype == "swiss":
#       LOGGER.info("Input is a Swissmodel ID. Locating coordinate file...")
#       return swiss_lookup(io,entity)
#     elif etype == "unp":
#       unp_flag = True
#       LOGGER.info("Input is a UniProt AC. Identifying all relevant PDB structures and ModBase models...")
#       return uniprot_lookup(io,entity)
#     elif etype: # HGNC returned a UniProt ID
#       unp_flag = True
#       LOGGER.info("Input is an HGNC gene name. Identifying all relevant PDB structures and ModBase models...")
#       LOGGER.info("HGNC: %s => UniProt: %s..."%(entity,etype))
#       return uniprot_lookup(io,etype)
#     else:
#       msg = "Could not identify %s."%entity
#       LOGGER.critical(msg);
#       sys_exit_failure(msg)

          
# The parser has to cope with AnnnB, AaaNNNBbb, and up to chain_id|transcript_id:p.AaaNNNBbb (official) input formats
# See http://varnomen.hgvs.org/
# Strings ending with .X are consired to terminate with chain indicators
from Bio.SeqUtils import seq1
class PDBMapVariantSet():
    """Manage a set of variants, i.e. pathogenic or benign,
    in the context of their transcript or perhaps specific
    chain identifier"""
    def __init__(self,variants_flavor = None):
        self._chain_id = None  # Chain letter code from a PDB entry
        self._unp_id = None    # Uniprot ID for these variants
        self._enst_id = None   # Ensembl Transcript ID for these variants
        self._gene_id = None   # Gene ID for these variants
        self._variants = []
        self._variants_flavor = variants_flavor

    @property
    def variant_list(self):
        return self._variants

    @property
    def variants_flavor(self):
        return self._variants_flavor

    def __repr__(self):
        
        if self._unp_id:
            repr_string = "Uniprot " + self._unp_id
        elif self._enst_id:
            repr_string = self._enst_id # Will be ENST.....
        elif self._gene_id:
            repr_string = "Gene " + self._gene_id 
        else:
            repr_string = "Unknown"

        return self.variants_flavor + " " + repr_string + ":\n" + str(self._variants)

    _one_letter_matcher = re.compile(r"^(?P<refAA>[A-Z])(?P<position>[0-9]+)(?P<altAA>[A-Z])$")
    assert _one_letter_matcher

    _three_letter_matcher = re.compile(r"^(?P<refAA>[A-Z]{3})(?P<position>[0-9]+)(?P<altAA>[A-Z]{3})$")
    assert _three_letter_matcher

    def query_and_extend(variants_flavor: str, id_to_variant_set,variant_set_id,ENST_transcripts,label): # ,**kwargs):
        """ Default string for variant queries """

        # Label and transcript are fast indices into genomic consequence
        if not ENST_transcripts:
            LOGGER.critical("Cannot execute %s query without ENST transcripts"%label)

        # when we're done we'll have nice "in" string so that we get ALL the transcripts in one query
        transcript_in_str = 'transcript in ('
        first_one = True
        ENST_transcripts_str = ''  # Build up a string that can be part of the dump'd query .tsv filename
        for transcript in ENST_transcripts:
            if not first_one:
                transcript_in_str += ','
            transcript_in_str += "'%s'"%transcript.id
            if len(ENST_transcripts_str) > 1:
                ENST_transcripts_str += '_'
            ENST_transcripts_str += transcript.id
            first_one = False
        transcript_in_str += ')'

        # import pdb; pdb.set_trace()
        
        query_str = ("SELECT distinct GC.protein_pos,GC.ref_amino_acid,GC.alt_amino_acid,\n"
                     "GC.chrom as chrom,GC.pos as pos,GC.end as end,GC.id as id,\n"
                     "GC.transcript,GC.allele,GC.ref_codon,GC.alt_codon\n"
                     " FROM GenomicConsequence AS GC\n")
        if label == 'clinvar': 
            query_str += "INNER JOIN clinvar on GC.chrom = clinvar.chrom and GC.pos = clinvar.pos\n"
            query_str += "  AND clinvar.clnsig %s ('Likely_pathogenic','Pathogenic','Pathogenic/Likely_pathogenic')" %('NOT IN' if variants_flavor == 'Neutral' else 'IN')
        elif label == 'clinvar38': 
            query_str += "INNER JOIN clinvar38 on GC.chrom = clinvar38.chrom and GC.pos = clinvar38.pos\n"
            query_str += "  AND clinvar38.clnsig %s ('Likely_pathogenic','Pathogenic','Pathogenic/Likely_pathogenic')" %('NOT IN' if variants_flavor == 'Neutral' else 'IN')
        elif label == 'drug':
            query_str += "INNER JOIN clinvar on GC.chrom = clinvar.chrom and GC.pos = clinvar.pos"
            query_str  += " AND clinvar.clnsig like '%drug%'\n"
            label = 'clinvar' # <- careful.  This is a clinvar query in the end: clinvar drug variants...
        elif label == 'cosmic': 
            query_str += "INNER JOIN cosmic on GC.chrom = cosmic.chrom and GC.pos = cosmic.pos"
            query_str  += " AND cosmic.cnt > 1\n"  # COSMIC queries only include count > 1
        elif label == 'cosmic38': 
            query_str += "INNER JOIN cosmic38 on GC.chrom = cosmic38.chrom and GC.pos = cosmic38.pos"
            query_str  += " AND cosmic38.cnt > 1\n"  # COSMIC38 queries only include count > 1
        elif label == 'tcga': 
            query_str += "INNER JOIN tcga on GC.chrom = tcga.chrom and GC.pos = tcga.pos"
        elif label == 'gnomad' or label == 'gnomad38': # or label == 'gnomad38': 
            query_str += "INNER JOIN GenomicData ON GC.label = GenomicData.label AND GC.chrom = GenomicData.chrom"
            query_str += " AND GC.pos = GenomicData.pos AND GC.end = GenomicData.end"
            query_str += " AND GenomicData.maf >= 1E-5 \n"  # GNOMAD only include maf of 1 in 10,000 or more to be neutral

        query_str += ( " WHERE GC.label=%s and " + transcript_in_str +
                       " and consequence LIKE '%%missense_variant%%' and length(ref_amino_acid)=1 and length(alt_amino_acid)=1 ")

        query_str += "ORDER BY GC.protein_pos,GC.ref_amino_acid,GC.alt_amino_acid"

        """ # With kwargs, a caller can append "and pathogenic=1' to the query
        for kwarg in kwargs:
            if isinstance(kwargs[kwarg],int):
                query_str += "and %s = %d"%(kwarg,kwargs[kwarg])
            else:
                query_str += "and %s = '%s'"%(kwarg,kwargs[kwarg])"""

        force_chain = None # <- For all these SQL queries there is no chain
        prior_variant_count = 0
        extended_variant_count = 0
        with PDBMapSQLdb() as db:
            db.execute(query_str,(label,))
            column_names = [x[0] for x in db.description]
            protein_changes_and_genomic_coordinates_df = pd.DataFrame(db.fetchall(),columns=column_names)
            # protein_changes_and_genomic_coordinates_df = pd.read_sql(query_str,db,params=(label,))
            LOGGER.info("%d rows returned from query",len(protein_changes_and_genomic_coordinates_df))
            csv_logfile = os.path.join(args.outdir,variants_flavor + '_' + label + ENST_transcripts_str + ".tsv")
            LOGGER.info("Writing all rows to tsv file %s:"%csv_logfile)
            protein_changes_and_genomic_coordinates_df.to_csv(csv_logfile,sep='\t')
            # if rows: If we have no rows, fine - then the variant_set is an emoty set
            # Add to variants for this flavor, creating a new set if needed
            _variant_set = id_to_variant_set.get(variant_set_id,PDBMapVariantSet(variants_flavor))

            prior_variant_count = len(_variant_set._variants)
            for index,row in protein_changes_and_genomic_coordinates_df.iterrows():
                _variant_set._variants.append([int(row['protein_pos']),row['ref_amino_acid'],row['alt_amino_acid'],force_chain])
            extended_variant_count = len(_variant_set._variants)

            if variant_set_id in id_to_variant_set:
                pass  # Because variant_set is the reference already in the dictionary from before
            else: # It's a new PDBMapVariantSet() instance, so add to dictionary
                id_to_variant_set[variant_set_id] = _variant_set

        return extended_variant_count - prior_variant_count

    def to_DataFrame(self,chain_id:str=None) -> pd.DataFrame:
        df = pd.DataFrame(self.variant_list,columns=["unp_pos","ref","alt","chain"]) ## ,"dcode","qt"])
        return df

    @staticmethod
    def parse_variants(variants_flavor: str,variants_or_filename): # -> Dict[str,PDBMapVariantSet]:
        """Parse comma,tab, or new-line delimited variants in form
           {chain|unp|enst}:[.AaaNNNNBbb

           parameters:
              flavor: string such as 'pathogenic' or 'benign'
              variants_or_filename: string or filename containing the variants"""

        if not variants_or_filename:
            return {}

        id_to_variant_set = {}

        if os.path.isfile(variants_or_filename):
            LOGGER.info("Parsing %s variants from file: %s",variants_flavor,variants_or_filename)

            raw_variants = [variant for line in open(variants_or_filename,'rt') for variant in line.strip().split(',')] 
        else:
             LOGGER.info("Parsing %s variants from command line argument:\n %s",variants_flavor,variants_or_filename)
             raw_variants = [variant.strip() for variant in variants_or_filename.split(',')]

        for raw_variant in raw_variants:
            dict_key = None

            # Step 1: Pull out the transcript or chain ID as a dict_key if you find one.  
            split_at_colon = raw_variant.split(':')
            variant_str_remaining = raw_variant.upper()
            if len(split_at_colon) > 2:
                msg = "%s variant %s has more than 1 colon.  Execution halting."%(flavor,raw_variant)
                LOGGER.critical(msg)
                sys_exit_failure(msg)
            if len(split_at_colon) == 2:
                dict_key = split_at_colon[0]   # Could be a chain, unp, or enst...
                variant_str_remaining = split_at_colon[1].upper()

            # Add to variants for this flavor, creating a new set if needed
            variant_set = id_to_variant_set.get(dict_key,PDBMapVariantSet(variants_flavor))

            if variant_str_remaining.startswith("P."):
                variant_str_remaining = variant_str_remaining[2:]

            force_chain = None
            # There _could_ be a chain ID at the very end - which is rarely used holdover
            chain_delimiter = variant_str_remaining.rfind('.')
            if chain_delimiter > -1:
                force_chain = (variant_str_remaining[chain_delimiter+1])
                variant_str_remaining = variant_str_remaining[0:chain_delimiter]
            else: # It could ALSO be the case that the chain was to the left before colon
                if dict_key and len(dict_key) < 6:
                    force_chain = dict_key

            result = PDBMapVariantSet._one_letter_matcher.match(variant_str_remaining)
            if result: # Then we parsed XnnY
                v_dict = result.groupdict()
            else:
                result = PDBMapVariantSet._three_letter_matcher.match(variant_str_remaining)
                if result:
                    v_dict = result.groupdict()
                    v_dict['refAA'] = seq1(v_dict['refAA'])
                    v_dict['altAA'] = seq1(v_dict['altAA'])
                else:
                    msg = "Unable to parse variant input %s, specifically fragment [%s]"%(raw_variant,variant_str_remaining)
                    LOGGER.critical(msg)
                    sys_exit_failure(msg)

            variant_set._variants.append([int(v_dict['position']),v_dict['refAA'],v_dict['altAA'],force_chain])
            if dict_key in id_to_variant_set:
                pass  # Because variant_set is the reference already in the dictionary from before
            else: # It's a new PDBMapVariantSet() instance, so add to dictionary
                id_to_variant_set[dict_key] = variant_set
            #else we were manipulating a reference already in the dictionary
       
        retained_variant_count =  sum([len(id_to_variant_set[dict_key].variant_list) for dict_key in id_to_variant_set])
        if len(raw_variants) == retained_variant_count:
            LOGGER.info("Raw variant parse: %d %s variants loaded",len(variant_set.variant_list),variants_flavor)
        else:
            LOGGER.warning("Blank/invalid variants dected in raw variant parse.  Started with %d %s, then reduced to %d",len(raw_variants),variants_flavor,retained_variant_count)

        return id_to_variant_set


def parse_qt(varset):
  if not varset:
    return []
  var,q = [],[]
  with open(varset,'r') as f:
    for line in f:
      var.append(line.split('\t')[0].strip())
      q.append(line.split('\t')[1].strip())
  var = ','.join(var)
  var = PDBMapVariantSet.parse_variants("Quantity weights",var)
  return list(zip(var,q))


@np.vectorize  # First argument is a vector, np.vectorize hands each element to function
def NeighborWeight(vus_to_known_distance_vectorized,lower_bound=8.,upper_bound=24):
    """For each distance between a variant of unknown significance (VUS) and a
    known (pathogenic or neutral) variant, return a spatial proximity metric that 
    weights by emphasizing closeness to the lower_bound of the interval of interest

    See page 2 of 10 in "Methods:Protein structural" analysis in... 
    "Sivley et al. BMC Bioinformatics (2018) 19:18" 
    which notes,
    "We quantified the spatial proximity of each VUS to each known pathogenic and entural variants
     using the NeighborWeight transformation of the 3D Euclidean distance between the centroid
     of each amino acid side chain"
    

    Parameters
    ----------
    d: Array-like set of distances between centroids of a VUS and known variants
       Code treats the variable like scalar, iteration controlled by @np.vectorize
    
    Returns
    -------
    Same shape'd np.array with weights in the interval [0.0,1.0]

    In practice, the lower_bound can equal the upper_bound, as when they
    calculated by a preceding the Ripley's K analysis"""

    
    if vus_to_known_distance_vectorized <= TOL: # Inter-residue distances this small are indicative of "same residue to same residue" or problems
        retval = 1+1e-10
        LOGGER.warning("NeighborWeight encountered tiny inter-residue distance %g Returning %g",vus_to_known_distance_vectorized,retval)
        return retval
    elif vus_to_known_distance_vectorized <= lower_bound:
        return 1.
    elif vus_to_known_distance_vectorized >= upper_bound:
        return 0.
    else:
        return 0.5*(np.cos(np.pi*(vus_to_known_distance_vectorized-lower_bound)/(upper_bound-lower_bound))+1)

def uniprox(cands,v,nwlb=8.,nwub=24.,weights=None):
    """ Predicts a univariate constraint score from weight vector 'v' """
    if weights:
      Cw,Vw = weights
      Cw = Cw.values.reshape((Cw.size,1))
      Vw = Vw.values.reshape((Vw.size,1))
    distance_matrix = cdist(cands,v)   # Calculate Distance matrix
    distance_matrix[distance_matrix==0] = np.inf     # Do not use variants to score their own residues
    NWv = NeighborWeight(distance_matrix,lower_bound=nwlb,upper_bound=nwub)
    if weights:
      # Substitution severity adjustment
      NWv = NWv * Vw.T# * np.dot(Cw,Vw.T)
    # np.fill_diagonal(NWv,0.) # NeighborWeight of self = 0.0. Inplace.
    cscores = np.sum(NWv,axis=1)/v.size
    return list(cscores)

def pathprox(cands,path,neut,nwlb=8.,nwub=24.,cross_validation_flag=None,weights=None):
  """ Predicts pathogenicty of variants of unknown significance (VUS)s
  given their coordinates, as well as coordinates of known pathogenic
  and neutral resides"""
  if weights:
    Cw,Aw,Bw = weights
    Cw = Cw.values.reshape((Cw.size,1))
    Aw = Aw.values.reshape((Aw.size,1))
    Bw = Bw.values.reshape((Bw.size,1))
  ccount = len(cands)
  pcount = len(path)
  ncount = len(neut)
  N = pcount + ncount
  # NeighborWeight of each candidate for each neutral/pathogenic variant
  NWn = NeighborWeight(cdist(cands,neut),lower_bound=nwlb,upper_bound=nwub)
  NWp = NeighborWeight(cdist(cands,path),lower_bound=nwlb,upper_bound=nwub)
  if weights:
    # Substitution severity of each candidate for each neutral/pathogenic variant
    NWn  = NWn * Bw.T# * np.dot(Cw,Bw.T)
    # NWp  = NWp * Aw.T# * np.dot(Cw,Aw.T) # DO NOT WEIGHT PATHOGENIC VARIANTS  
  psum = np.sum(NWn)
  nsum = np.sum(NWp)
  # Set self-weights to 0. for cross-validation
  if   cross_validation_flag == "N":
    assert(len(set(NWn.shape))==1) # Fail if not square matrix
    np.fill_diagonal(NWn,0.) # NeighborWeight of self = 0.0. Inplace.
  elif cross_validation_flag == "P":
    assert(len(set(NWp.shape))==1) # Fail if not square metrix
    np.fill_diagonal(NWp,0.) # NeighborWeight of self = 0.0. Inplace.
  NWs = np.hstack((NWn,NWp)) # NeighborWeight matrix stack
  nmask = np.zeros(N).astype(bool)
  nmask[:ncount] = True
  pmask = np.abs(nmask-1).astype(bool)
  pscores = np.sum(NWs[:,pmask],axis=1)/pcount
  nscores = np.sum(NWs[:,nmask],axis=1)/ncount
  cscores = pscores - nscores # subtraction binds c to -1..1
  return list(cscores)

def calc_auc(pathogenic_scores,neutral_scores):
  """Calculate the Area under the ROC curve (roc_auc)

  Parameters
  ----------
  Pathogenic scores
  Neutral scores

  Returns
  -------
  False Positive Rates[]
  True Positive Rates[]
  roc_auc 
  precision,recall pairs
  average precision score
  """
 
  # labels are the known 'true binary labels' or 'y_true' parameter
  LOGGER.info("calc_auc starting with %d pathogenic scores and %d neutral scores",len(pathogenic_scores),len(neutral_scores))
  known_binary_labels = [0]*len(neutral_scores)+[1]*len(pathogenic_scores)

  # catenate the lists of scores we computed
  all_scores  = neutral_scores+pathogenic_scores

  # Compute the increasing false positive, and true positive rates
  # The ith element of fpr[i] and tpr[i]  is the false(true) positive rate of scores > threshholds[i]
  fpr,tpr,thresholds  = sklearn.metrics.roc_curve(known_binary_labels,all_scores,drop_intermediate=False)

  roc_auc    =          sklearn.metrics.auc(fpr,tpr)
  LOGGER.info("roc_auc=%f  Compare to roc_auc_score()=%f",roc_auc,sklearn.metrics.roc_auc_score(known_binary_labels,all_scores))

  # Compute precision-recall pairs Precision is #true_positives/(true_positives+false_positives)
  precisions,recalls,_ = sklearn.metrics.precision_recall_curve(known_binary_labels,all_scores)
  precision_max_accumulate        = np.maximum.accumulate(precisions)
  pr_auc     =           sklearn.metrics.average_precision_score(known_binary_labels,all_scores,average="micro")
  
  LOGGER.info("Returning roc_auc=%g  pr_auc=%g",roc_auc,pr_auc)
  return fpr,tpr,roc_auc,precision_max_accumulate,recalls,pr_auc

def plot_roc(fpr,tpr,ax=None,save=True,label="PathProx",color='k'):
  ## Plot the ROC curve
  roc_auc = sklearn.metrics.auc(fpr,tpr)
  if not ax:
    fig,ax = plt.subplots(1,1,figsize=(5,5))
    # plt.title("%s ROC"%label)
    ax.plot([0,1],[0,1],'k--')
    ax.set_xlim([0.,1.])
    ax.set_ylim([0.,1.])
    ax.set_xlabel("False Positive Rate",fontsize=14)
    ax.set_ylabel("True Positive Rate",fontsize=14)
  l = "%s AUC: %5.2f"%(label,roc_auc)
  ax.plot(fpr,tpr,color=color,label=l,linewidth=4)
  if save:
    sns.despine(ax=ax)
    ax.legend(loc="lower right",fontsize=10)
    plt.sca(ax)
    plt.savefig("%s_pathprox_roc.pdf"%(args.label),dpi=300)
    plt.savefig("%s_pathprox_roc.png"%(args.label),dpi=300)
  return ax

def plot_pr(rec,prec,pr_auc,ax=None,save=True,label="PathProx",color='k'):
  ## Plot the PR curve
  if not ax:
    # Return one 5"x5" matplotlib.figure.Figure and matplotlib.axes.Axes
    fig,ax = plt.subplots(1,1,figsize=(5,5))
    # plt.title("%s PR"%label)
    ax.set_xlim([0.,1.])
    ax.set_ylim([0.,1.])
    ax.set_xlabel("Recall",fontsize=14)
    ax.set_ylabel("Precision",fontsize=14)
  l = "%s AUC: %5.2f"%(label,pr_auc)
  ax.plot(rec,prec,color=color,label=l,linewidth=4)
  if save:
    sns.despine(ax=ax)
    ax.legend(loc="lower left",fontsize=10)
    plt.sca(ax)
    plt.savefig("%s_pathprox_pr.pdf"%(args.label),dpi=300)
    plt.savefig("%s_pathprox_pr.png"%(args.label),dpi=300)
  return ax

def confidence(auc):
  if auc < 0.6:
    return "Low"
  elif auc < 0.75:
    return "Moderate"
  else:
    return "High"

def predict(cscores):
  return [["Neutral","Deleterious"][int(float(c>0))] for c in cscores]

"""
Seems to be dead code (obserbaton of Chris Moth Feb 5 2018
def write_results(cands,cscores,cpvalues,preds,conf,N,mwu_p,roc_auc,pr_auc,label):
  if not cands.size: return
  # If evaluating candidate mutations, continue...
  h = ["mut","score","pvalue","N","mwu_p","roc_auc","pr_auc","conf","prediction"]
  with open("%s_%s_results.txt"%(args.label.replace(' ','_'),'wb'),args.radius) as fout:
    writer = csv.writer(fout,delimiter='\t')
    writer.writerow(h)
    res = []
    for i in range(len(cscores)):
      row = [cands[i]]+["%.4f"%cscores[i],"%.3g"%cpvalues[i],"%d"%N,"%.3g"%mwu_p,"%.2f"%roc_auc,"%.2f"%pr_auc,conf,preds[i]]
      res.append(row)
    res = sorted(res,key=lambda x: float(x[1]),reverse=True)
    LOGGER.info("\t".join(h[:3]+h[-2:]))
    for row in res:
      writer.writerow(row)
      LOGGER.info("\t".join(row[:3]+row[-2:]))
"""

def var2coord(s,p,n,c,q=[]):
  # Construct a dataframe containing all variant sets
  # These new columns are initialized to np.nan and only change from nan
  # If they are loaded from known path/neutral/candidate/variants
  vdf = pd.DataFrame(columns=["unp_pos","ref","alt","chain","dcode","qt"])

  if p: # List of pathogenic variants
    pdf = pd.DataFrame(p,columns=["unp_pos","ref","alt","chain"])
    pdf["dcode"] = 1
    pdf["qt"]    = np.nan
    vdf = vdf.append(pdf)
  if n: # List of neutral variants
    ndf = pd.DataFrame(n,columns=["unp_pos","ref","alt","chain"])
    ndf["dcode"] = 0
    ndf["qt"]    = np.nan
    vdf = vdf.append(ndf)
  if c: # List of Candidate variants
    cdf = pd.DataFrame(c,columns=["unp_pos","ref","alt","chain"])
    cdf["dcode"] = -1
    cdf["qt"]    = np.nan
    vdf = vdf.append(cdf)
  if q: # List of quantified variants
    qdf = pd.DataFrame(q,columns=["unp_pos","ref","alt","chain","qt"])
    qdf["dcode"] = np.nan
    qdf["qt"]    = qdf["qt"].astype(np.float64)
    # Update column order to match  binary datasets
    qdf = qdf[["unp_pos","ref","alt","chain","dcode","qt"]]
    vdf = vdf.append(qdf)

  # drop_duplicates considers only the columns, not the IntIndex elements
  vdf.drop_duplicates(inplace=True)
 
  # Create a new MonotonicRangeIndex from 0 to len()-1.  Old index will be discarded
  vdf.reset_index(drop=True,inplace=True)

  msg = None
  if vdf.empty and not (
      args.add_exac or args.add_gnomad or args.add_gnomad38 or args.add_1kg or args.add_pathogenic or args.add_pathogenic38 or args.add_cosmic or args.cosmic38 or args.add_tcga):
    msg = "\nERROR: Must provide variants or request preloaded set with --add_<dataset>.\n"
  elif vdf.empty:
    msg = "\nERROR: No variants identified. Please manually provide pathogenic and neutral variant sets.\n"
  if msg:
    LOGGER.critical(msg); 
    sys_exit_failure(msg)

  # If a position has both pathogenic and neutral variants, THEN eliminate the neutral ones
  def defer_to_pathogenic(g):
    # If all sub-group elements are neutral or candidate VUSs, then 
    # simply return all of them.  Otherwise, return a group that has
    # excluded all neutrals, and is _only_ pathogenic or candidates 
    return g if all(g["dcode"]<=0) else g[g["dcode"]!=0]

  # vdf_neutrals_deferred = vdf.groupby(["unp_pos","ref","alt","chain"],as_index=False).apply(defer_to_pathogenic)
  vdf_neutrals_deferred = vdf.groupby(["unp_pos","ref","chain"],as_index=False).apply(defer_to_pathogenic)
  
  # Construct a dataframe from the coordinate file
  complex_df = pd.DataFrame([[r.get_parent().id,r.id] for r in s.get_residues()],columns=["chain","pdb_pos"])
  complex_df["unp_pos"] = [r.get_parent().alignment.pdb2seq[r.id] for r in s.get_residues()]

  # Calculate the coordinate center-of-mass for all residues
  def resicom(resi):
      if resi==None: return [None,None,None]
      return np.array([np.array(a.get_coord()) for a in resi.get_unpacked_list()]).mean(axis=0)

  coords   = [resicom(r) for r in s.get_residues()]
  coord_df = pd.DataFrame(coords,index=complex_df.index,columns=["x","y","z"])
  complex_df_coord = complex_df.merge(coord_df,left_index=True,right_index=True)
   
  # Merge the sequence-aligned structure dataframe with the sequence-based variant dataframe
  vdf_extract = vdf_neutrals_deferred[["unp_pos","ref","alt","chain","dcode","qt"]].copy()
  vdf_extract['unp_pos'] = vdf_extract['unp_pos'].apply(int)
  vdf_extract['chain'] = vdf_extract['chain'].apply(str)
  complex_df_coord['unp_pos'] = complex_df_coord['unp_pos'].apply(int)
  complex_df_coord['chain'] = complex_df_coord['chain'].apply(str)
  complex_df_merged = complex_df_coord.merge(vdf_extract,how="left",on=["chain","unp_pos"])
  # Ensure that no duplicate residues were introduced by the merge
  complex_df_merged.reset_index(drop=True,inplace=True)
  complex_df = complex_df_merged.drop_duplicates(["unp_pos","pdb_pos","chain","ref","alt","dcode"]).reset_index(drop=True)

  # Annotate each amino acid substitution with the blosum100 score
  complex_df["blosum100"] = complex_df.apply(lambda x: blosum100[(x["ref"],x["alt"])] if not np.isnan(x["dcode"]) else None,axis=1)

  # Fail if all variants outside structural coverage or mapped to missing residues
  msg = None
  vdf = complex_df[~complex_df["dcode"].isnull()]
  if vdf.empty:
    if not complex_df.empty:
      msg = "Structure successfully aligned to sequence, but no variants were mapped to non-missing residues."
    else:
      msg = "ERROR: Sequence-structure alignment failed."
    LOGGER.critical(msg); 
    sys_exit_failure(msg)


  # Check that both variant categories are populated (if input is categorical)
  if not args.quantitative:
    msg = None
    if (complex_df["dcode"]==0).sum() < 3:
      msg = "\nWARNING: Structure contains %d neutral variants (PathProx minimum 3).\n"%(complex_df["dcode"]==0).sum()
      if not args.add_exac and not args.add_gnomad and not args.add_1kg:
        msg += "Consider using --add_exac, --add_gnomad, or --add_1kg.\n"
      else:
        msg += "Please manually specify neutral variants.\n"
    if msg:
      LOGGER.warning(msg)
    if (complex_df["dcode"]==1).sum() < 3:
      msg = "\nWARNING: Structure contains %d pathogenic variants (PathProx minimum 3).\n"%(complex_df["dcode"]==1).sum()
      if not args.add_pathogenic and not args.add_cosmic and not args.add_tcga:
        msg += "Consider using --add_pathogenic, --add_cosmic, or --add_tcga.\n"
      else:
        msg += "Please manually specify pathogenic variants.\n"
      if msg:
        LOGGER.warning(msg)

  # Conversion from short-code to explicit description
  code2class = {-1 : "Candidate",
                 0 : "Neutral",
                 1 : "Pathogenic"}
  cnts   = complex_df.drop_duplicates(["unp_pos","ref","alt","dcode"]).groupby("dcode").apply(lambda x: (x["dcode"].values[0],len(x)))
  LOGGER.info("################################")
  LOGGER.info("Unique Variant Counts:")
  with open("%s_variant_counts.txt"%args.label,'w') as fout:
    writer = csv.writer(fout,delimiter='\t')    
    for dcode,cnt in cnts[::-1]:
      writer.writerow([code2class[dcode],cnt])
      LOGGER.info("%20s:  %d"%(code2class[dcode],cnt))
  LOGGER.info("")
  return complex_df

def TrangeCondensed(D):
  """From a compressed pdist() type list of inter-residue distances, return
     list of distances to consider"""

  min_observed_intervariant_distance = np.min(D)
  minT = max(np.ceil(min_observed_intervariant_distance),5)          # minimum observed inter-variant distance.  But no less than 5A
  maxT = min(np.ceil(np.max(D)),45)            # maximum observed inter-variant distance (bivariate) if <45A
  if maxT <= minT:
    maxT = np.ceil(np.max(D))
    msg = "Observations too close (min=%.0f, max=%.0f) to divide space; using full range.\n"%(minT,maxT)
    LOGGER.warning(msg)
  # Verify that the structure is large enough to analyze multiple distances
  if maxT == minT:
    LOGGER.warning("Skipped %s.%s: Structure is too small to analyze."%(sid,chain))
    return None
  # Retun a list of distances to consider, as consecutive Angstrom gradatons
  T = np.arange(minT,maxT+1,1) # inclusive range
  return T

def Trange(D,o=[]):
  """From the square symmetric matrix of inter-residue distances
  Calculate the minimum and maximum inter-variant distance Threshhlds.

  The optional parameter o is a specific list of variants
  that are of interest, without which all residues are considered

  The return value is a list of integers, the range of distance thresholds
  from minT to maxT"""
  
  min_observed_intervariant_distance = np.min(D[np.triu_indices(D.shape[0],k=1)])

  minT = max(np.ceil(min_observed_intervariant_distance),5)          # minimum observed inter-variant distance.  But no less than 5A
  # If a specific list of variants was provided, 
  # then consider only those variants in computing the max inter-variant distance
  if any(o):
    o = o.astype(bool)
    maxT = min(np.ceil(np.nanmax(D[o,:][:,o])),45)  # maximum observed inter-variant distance (univariate) if <45A
  else:
    maxT = min(np.ceil(np.nanmax(D)),45)            # maximum observed inter-variant distance (bivariate) if <45A
  if maxT <= minT:
    maxT = np.ceil(np.n)
    msg = "Observations too close (min=%.0f, max=%.0f) to divide space; using full range.\n"%(minT,maxT)
    LOGGER.warning(msg)
  # Verify that the structure is large enough to analyze multiple distances
  if maxT == minT:
    LOGGER.warning("Skipped %s.%s: Structure is too small to analyze."%(sid,chain))
    return None
  # Retun a list of distances to consider, as consecutive Angstrom gradatons
  T = np.arange(minT,maxT+1,1) # inclusive range
  return T

def permute(y,N):
  """ Generate permutations of the rows of vector/matrix y 
  Parameters
  ----------
  y: the vector or matrix that should be shuffled

  N: The number of shuffles that the generator must perform
  """
  for i in range(N):
    yield np.random.permutation(y)
    
def Kest(D,y,T=[],P=99999):
  """ Ripley's K-Function Estimator for Spatial Cluster Analysis (w/ Positional Constraints) (w/o Edge Correction)

  Parameters:
  -----------
  D: Distance matrix for all possible point pairs (observed and unobserved)

  y: Weight (or "use these residues") vector for all possible centroids (un-observed points must have NaN weight)

  T: Distance thresholds (integer range of distances of interest)

  P: Number of permutations to simulate the empiracal null distribution, and calculate a 95% confidence envelope 
    P=None implies a recursive call to calculate _ONE_ new permutation

  Caveat: Current implementation only handles positive weights
  
  Returns:
  --------
  K(t) is returned as in the paper and defined as
    K(t) = (SumOver(i, j != i)->(I(Dij < t))/N(N-1) for all thresholds t in T

  """

  if P: # Only log in the outer call, not all the recursive permutations
    LOGGER.info("Begin Kest() D.shape=%s (y > 0).sum=%s permutations=%s"%(str(D.shape),(y > 0).sum(),str(P)))
  assert(P!=1)

  # Ensure that all the weights are float values
  y = np.array(y,dtype=np.float64) # for NaN compatibility

  # The ys are deemed 'weighted' if they contain something
  # other than 0s, 1s, or np.Nans.
  weighted = (y==0).sum()+(y==1).sum() != y[~np.isnan(y)].size
  if not weighted: # convert 0/1 to nan/1, a cleanup of sorts
    y[y==0] = np.nan      # Convert all the 0s to Nans
    y[~np.isnan(y)] = 1.0 # Convert everything else to 1.0 floats

  o   = ~np.isnan(y) # o is a boolean np array
  Do  = D[o,:][:,o] # Distance matrix of only observed points, Nans on diagonal
  yo  = y[o]        # y values of only observed points
  R   = y.size      # Number of protein residues
  N   = yo.size     # Number of observed points
  if weighted:
    if P:
      Kest.DT = [Do>t for t in T] # Precompute distance masks
    Y  = np.outer(yo,yo) # NaN distance diagonal handles i=j pairs
    Y /= Y.sum() # normalize by all-pairs product-sum
    K  =  np.array([np.ma.masked_where(dt,Y).sum() for dt in Kest.DT])
  else:
    # K = np.array([((~np.isnan(Do)) & (Do<=t)).sum() for t in T],dtype=np.float64) / (N*(N-1))
    # The above original line emits many run-time warnings because of the nan in the diagnoal
    # Chris Moth changed Oct 12, 2018 to only use the upper triangle
    inter_residue_distances = Do[np.triu_indices(Do.shape[0],k=1)]
    K = np.array([(inter_residue_distances <= t).sum() for t in T],dtype=np.float64) / ((N * (N-1)) // 2)

  if P: # Then this is the "master" call to Kest, and about to return final results
    if weighted:
      # If weighted, shuffle values for observed points, via recursion
      K_perm = np.array([Kest(Do,yp,T,P=None) for yp in permute(yo,P)])
    else:
      # If unweighted, shuffle positions for all points, via recursion
      K_perm = np.array([Kest(D,yp,T,P=None) for yp in permute(y,P)])
    # Add the observed K vector to the permutation matrix
    K_perm = np.concatenate(([K],K_perm))
    # LOGGER.error("REMOVE THIS Writing to /tmp/temp.csv")
    # np.savetxt('/tmp/temp.csv',K_perm,delimiter='\t')
    # Calculate the simulated z-score 
    permutations_K_average,permutations_K_std,K_z = permutations_average_std_zscore(K_perm)
    if all(np.isnan(K_z)):
      # If all weights are equal, set K_z to 0. rather than NaN
      K_z = np.zeros(K_z.shape[0])
    # Calculate one-sided permutation p-value given K directionality
    K_p = []
    for i,z in enumerate(K_z):
      if z>0:
        K_p.append(min(1.,2.*(1.-percentileofscore(K_perm[:,i],K[i],'strict')/100.)))
      else:
        K_p.append(min(1.,2.*(1.-percentileofscore(-K_perm[:,i],-K[i],'strict')/100.)))

    # K_pz(t) = SurvivalFunction(o_z) = (1 - CumulativeDistributionFunction(o_z))*2.
    K_pz = norm.sf(abs(K_z))*2. # two-sided simulated p-value

    # Calculate the confidence envelope (the grey area on the K plot)
    # hce(t) = the values, below which 97.5% of the K_perm[t] fall 
    hce  = np.percentile(K_perm,97.5,axis=0)
    # lce(t) = the values, below which 2.5% of the K_perm[t] fall 
    lce  = np.percentile(K_perm, 2.5,axis=0)


    LOGGER.info(
      """End Kest() returning for t in (%d,%d):
           K(t)=%s
           K_p(t)=%s
           K_z(t)=%s
           K_pz=%s
           hce(t)=%s
           lce(t)=%s
           K_perm.shape=%s"""%(
             T[0],T[-1],
             K,K_p,K_z,K_pz,hce,lce,K_perm.shape))
    return K,permutations_K_average,permutations_K_std,K_p,K_z,K_pz,hce,lce,K_perm
  else: # P = None, and this is the recursion result for one permutation, so quickly return that permuted K
    return K

Kest.DT = []

def qtprox(cands,qt,nwlb=8.,nwub=24.):
  qt = qt[qt['qt'].notnull()]
  coordinates = qt[['x','y','z']].values
  qt = qt['qt'].values
  D = cdist(cands,coordinates)
  D[D==0] = np.inf # Do not scores residues using known values
  NWv = qt * NeighborWeight(D,lower_bound=nwlb,upper_bound=nwub)
  # np.fill_diagonal(NWv,0.)
  cscores = np.nansum(NWv,axis=1)/np.nansum(qt)
  return list(cscores)


def permutations_average_std_zscore(matrices):
  """Return the zscores of the 0th (observed)
     matrix in array of matrices.  Where standard
     deviations are 0, return 0 - rather than div by 0 exception"""

  # We do not use built in zscore is because
  # we have legit cases where the standard deviation is 0
  # (values are not distributed - as case of max_threshold
  # where 
  permutations_K_average = np.average(matrices,axis=0)
  permutations_K_std = np.std(matrices,axis=0)

  # This bool array will be true for std_matrix positions that are non-zero
  non_zero_stds = ~np.all(matrices == matrices[0],axis=0)

  zscores_or_0s = np.zeros(permutations_K_average.shape)

  # zscore = (observation-average)/sigma
  zscores_or_0s[non_zero_stds] = (matrices[0][non_zero_stds] - permutations_K_average[non_zero_stds]) / permutations_K_std[non_zero_stds]
  return permutations_K_average,permutations_K_std,zscores_or_0s


def pstats(Pmat):
  """ Calculates p-values and z-scores from the observation and
      (large) array of permutations of each 
      [Kpathogenic(t) - Kneutral(t)] matrix

      The First column of the matrix must contain the observation.
      Remaining (99999 etc) columns are the permutations
 """
  o = Pmat[0] # observed [Kpathogenic(t) - Kneutral(t)] matrix

  # Calculate the simulated z-score for the observed Pmat[0]
  _,_,o_z = permutations_average_std_zscore(Pmat)

  # 2018 October - Chris Moth removed below because 
  # the built in scipy zscore generates warnings on 
  # standard deviations of 0.
  # o_z = zscore(Pmat)[0]

  # Calculate one-sided permutation p-values

  o_p = []
  # Iterating over each Z(t)....
  # Directly calculate the p value from all the permutations
  for i,z in enumerate(o_z):
    if z>0: # i.e. if o[i] > meanOfPermutations[i], consoder samples to left
      # percentileofscore(series, observation) returns the percentage of 
      # values in the array 'series' which are 'strictly' below the observed
      # Kpathogenic(t) - Kneutral(t)]s
      # We double the 1-percentileofscore because normal distributions are two-sided
      # p values have max of 1, which should be rare given z test at outset
      o_p.append(min(1.,2*(1.-percentileofscore(Pmat[:,i],o[i],'strict')/100.)))
    else: # consider samples to right of observation
      # percentileofscore(-series, -observation) returns the percentage of 
      # values in the array 'series' which are 'strictly' above the observed
      # Kpathogenic(t) - Kneutral(t)]s
      # We double the 1-percentileofscore because normal distributions are two-sided
      # p values have max of 1, which should be rare given z test at outset
      o_p.append(min(1.,2*(1.-percentileofscore(-Pmat[:,i],-o[i],'strict')/100.)))

  # o_pz = SurvivalFunction(o_z) = 1 - CumulativeDistributionFunction(o_z)
  o_pz = norm.sf(abs(o_z))*2. # two-sided simulated p-value
  # Calculate the confidence envelope (High and Low Confidence Envelope)
  # hce(t) = the values, below which 97.5% of the Pmat[t] fall 
  hce = np.percentile(Pmat,97.5,axis=0)
  # lce(t) = the values, above which 97.5% of the Pmat[t] fall 
  lce = np.percentile(Pmat,2.5, axis=0)
  return o_p,o_z,o_pz,hce,lce

def k_plot(T,K,K_average,K_std,Kz,lce,hce,ax=None,weights=False):
  if not ax:
    fig,ax = plt.subplots(1,1,figsize=(15,5))
  # The 95% Confidence is the interval [lce(t),hce(t)] and colored grey
  ax.fill_between(T,lce,hce,alpha=0.2,
                      edgecolor='k',facecolor='k',
                      interpolate=True,antialiased=True)
  ax.plot(T,K_average,label='Average permuted K(t)',color='black',linestyle='dashed',linewidth=2.0) 
  ax.plot(T,K_average+K_std,label='+1 Std permuted K(t)',color='green',linestyle='dashed',linewidth=1.0) 
  ax.plot(T,K_average+2*K_std,label='+2 Std permuted K(t)',color='yellow',linestyle='dashed',linewidth=1.0) 
  ax.plot(T,K_average-K_std,label='-1 Std permuted K(t)',color='green',linestyle='dashed',linewidth=1.0) 
  ax.scatter(T,K,s=50,color='darkred',edgecolor='white',lw=1,label=["Un-Weighted Kobs(t)","Weighted Kobs(t)"][weights])
  # for t in range(K.shape[0]):
  #  text = ax.annotate("Kobs(t)=%0.2f, Kz(t)=%0.1f"%(K[t],Kz[t]), (T[t],K[t]),color='darkred',horizontalalignment='center',verticalalignment='top')
  #  text.set_rotation(+90.0)

    
  ax.set_xlabel("Distance Threshold t ($\\AA$)",fontsize=16)
  ax.set_ylabel("K(t) = Percentage of Inter-residue Distances Within t",fontsize=16,rotation=90)


  ax.set_xlim([T[0],T[-1]+0.5])
  # if any(K<0) or any(lce<0) or any(hce<0):
  #  ax.set_ylim([-1.05,1.05])
  # else:
  #  ax.set_ylim([-0.05,1.05])
  ax.set_ylim([-0.25,1.05])
  # Add a vertical line a the most extreme threshold
  dK = np.nanmax(np.abs(Kz),axis=0)    # maximum Kz
  t  = np.nanargmax(np.abs(Kz),axis=0) # t where Kz is maximized
  text = ax.annotate("Kobs(t)=%0.2f, Kz(t)=%0.1f"%(K[t],Kz[t]), (T[t],K[t]),color='darkred',horizontalalignment='right',verticalalignment='top')
  T,K = T[t],K[t]
  # ax.axhline(0,color='k',lw=1,ls="-")
  # This vertical dashed line on the grap shows 
  # the 't' on the X (i.e. T) axis where Kz is maximized
  ax.axvline(T,color='k',lw=2,ls="dashed",label="Most Significant K")
  return ax

def saveKplot(T,K,K_average,K_std,Kz,lce,hce,label="",weights=False):
  """Plot K(t) and its confidence envelope from

    Parameters
    ----------
    T:      Distance thresholds, the X axis
    K(t):   Ripley's K(t), fraction (probability) of residues with distance <= T(T)
    K_average(t): The average of the K(t) for all the permutations
    K_std(t): The standard deviation K(t) for all the permutations
    Kz:     The t where the zscore of K(t) is maximal, the 
    lce(t): Points above which 97.5% of the permuted K(t) were found
    hce:    Points below which 97.5% of the permuted K(t) were found
    label:  
  """
  global sid
  fig,ax = plt.subplots(1,1,figsize=(16,8))
  k_plot(T,K,K_average,K_std,Kz,lce,hce,ax)
  sns.despine()
  ax.set_title("Ripley's Kobs(t)",fontsize=16)
  ax.legend(loc="upper left",fontsize=14)
  plt.savefig("%s_%s_K_plot.pdf"%(args.label,label),dpi=300,bbox_inches='tight')
  plt.savefig("%s_%s_K_plot.png"%(args.label,label),dpi=300,bbox_inches='tight')
  plt.close(fig)

def d_plot(T,D,Dz,lce,hce,ax=None,weighted_flag=False):
  if not ax:
    fig,ax = plt.subplots(1,1,figsize=(15,5))
  # 95% Confidence
  ax.fill_between(T,lce,hce,alpha=0.2,
                      edgecolor='k',facecolor='k',
                      interpolate=True,antialiased=True)
  ax.scatter(T,D,s=50,color='darkred',edgecolor='white',lw=1,label=["Un-Weighted D","Weighted D"][weighted_flag])
  ax.set_xlabel("Distance Threshold (t)",fontsize=16)
  ax.set_ylabel("D",fontsize=16,rotation=90)
  ax.set_xlim([min(T),max(T)])
  if any(D<0) or any(lce<0) or any(hce<0):
    ax.set_ylim([-1.05,1.05])
  else:
    ax.set_ylim([-0.05,1.05])
  # Add a vertical line a the most extreme threshold
  dD = np.nanmax(np.abs(Dz),axis=0)    # maximum Kz
  t  = np.nanargmax(np.abs(Dz),axis=0) # t where Kz is maximized
  T,D = T[t],D[t]
  # ax.axhline(0,color='k',lw=1,ls="-")
  ax.axvline(T,color='k',lw=2,ls="dashed",label="Most Significant D")
  return ax

def saveDplot(T,D,Dz,lce,hce,label="",weights=False):
  global sid
  fig,ax = plt.subplots(1,1,figsize=(15,5))
  d_plot(T,D,Dz,lce,hce,ax)
  sns.despine()
  ax.set_title("Ripley's D",fontsize=16)
  ax.legend(loc="lower right",fontsize=14)
  plt.savefig("%s_D_plot.pdf"%args.label,dpi=300,bbox_inches='tight')
  plt.savefig("%s_D_plot.png"%args.label,dpi=300,bbox_inches='tight')
  plt.close(fig)

def uniK_calc_and_plot(D,y,P=99999,label="",weights=False):
  """ Univariate (i.e. neutral or pathogenic) K(t)
  calculation, with plot, and final return of K statistics
  of highest Z score.

  Parameters:
  -----------
  D:the symmetric inter-residue distance matrix
  for all pairs of residue centroids in the 3D structure,
  regardless of neutral/pathogenic/mutation/quantitative
  labels

  y:labels specific residues of interest for K(t) calculation
  If None, then all D residues are included in K(t) calculation

  P:count of requested permutations to create the empirical 
      null distribution
  label:string for output on the final result plot filenames

  Internally Computes and Plots:
  ------------------------------
     K(t):   Ripley's K statistic for residues

     Kp(t):  P values for each K(t) relative to Empirical Null permutations.

     max 


  Returns:
  --------
  K

  """
  LOGGER.info("Begin uniK_calc_and_plot(%s)  permutations=%d"%(label,P))
  # Get a list of Distance thresholds to test
  T = Trange(D,y.values)

  K,K_average,K_std,Kp,Kz,Kzp,hce,lce,_ = Kest(D,y,T,P=P)

  # Save the multi-distance K plot  
  saveKplot(T,K,K_average,K_std,Kz,lce,hce,label,weights=weights)

  # Determine the optimal K/T/p
  maxMagnitudeKzIndex = np.nanargmax(np.abs(Kz))
  K_ofMaxKz  = K[maxMagnitudeKzIndex]
  T_ofMaxKz  = T[maxMagnitudeKzIndex]
  P_ofMaxKz  = Kp[maxMagnitudeKzIndex]
  Kz_ofMaxKz = Kz[maxMagnitudeKzIndex]

  LOGGER.info("Ending uniK_calc_and_plot(%s) Returning for t of MaxKzScore: K=%.2f Kz=%.2f T=%.2f P=%.2f"%(label,K_ofMaxKz,Kz_ofMaxKz,T_ofMaxKz,P_ofMaxKz))

  return T,K_ofMaxKz,Kz_ofMaxKz,T_ofMaxKz,P_ofMaxKz

def biD(pathogenicCoordinates,neutralCoordinates,P=99999,label=""):
  """Compute Pathogenic less Neutral Bivariate D
  This is empirical null is calculated, via P permutations
  Arguments:
  array of (x,y,z) pathogenic variant residue centroids
  array of (x,y,z) neutral variant residue centroids
  P = number of permutations requested
  label= text for plot
  """

  LOGGER.info("Starting biD() pathogenicCoordinates.shape=%s  neutralCoordinates.shape=%s  Permutations=%d",str(pathogenicCoordinates.shape),str(neutralCoordinates.shape),P)

  # Counts of pathogenic and neutral variants
  nPathogenic = pathogenicCoordinates.shape[0] # Pathogenic variants
  nNeutral = neutralCoordinates.shape[0] # Neutral variants

  # combine neutral and pathogenic for a complete set of x,y,z spactial coordinates
  allCoordinates = np.concatenate((pathogenicCoordinates,neutralCoordinates))               

  # Initially, we have 1s (later unmasked) for the first (pathogenic) varian rows, 
  # 0s to later mask the neutrals
  maskPathogenic = np.array([True]*nPathogenic+[False]*nNeutral,dtype=bool)
  # Flip 0s to 1s, and vice versa so that 0=unmask

  # Neutrals are simply an inverse mask
  maskNeutral = ~maskPathogenic

  # pdist returns a condensed distance matrix as vectir
  allCoordinatesPdistVector = pdist(allCoordinates)  
  # squareform expands pdist vector to symmetric full matrix
  allCoordinatesDistanceMatrix  = squareform(allCoordinatesPdistVector)
  np.fill_diagonal(allCoordinatesDistanceMatrix,np.nan);

  # Importantly, the residue positions in the allCoordinates... will not change
  # in this calculation.  Only the flags (masks) of whether residues are pathogenic
  # or neutral are permuted

  # Distance thresholds to test.  See references for more info on 'T'
  T = TrangeCondensed(allCoordinatesPdistVector)

  ## Distance-dependent K derivatives
  # D is a compressed, masked, vector of inter-residue distances
  def biKest(distancesMask,N,T=[]):
    """ Abbreviated K function for bivariate analyses """
    # N = number of unmasked elements - was np.sqrt(D.count()) # qrt of Number of non-masked elements
    maskedDistanceMatrix = allCoordinatesDistanceMatrix[np.ix_(distancesMask,distancesMask)]
    return np.array([np.less(maskedDistanceMatrix[np.triu_indices(N,k=1)],t).sum() for t in T],dtype=np.float64) / (N*(N-1)//2)

  np.warnings.filterwarnings('error')

  KPathogenic  = biKest(maskPathogenic,nPathogenic,T)
  KNeutral  = biKest(maskNeutral,nNeutral,T)
  KPathogenicLessKNeutral = KPathogenic - KNeutral

  # Random label (neutral vs pathogenic) shuffling permutation test
  # .. which recalulates KPathogenic and KNeutral when entirely different residues
  # are selected as the Pathogenic and Neutral ones
  KPathogenicLessKNeutral_list = [KPathogenicLessKNeutral] # Initialize with observations

  for maskPathogenic_p in permute(maskPathogenic,P):
    maskNeutral_p = ~maskPathogenic_p

    # Calculate each multivariate K statistic
    # When mask entries are false that means that
    # _are_ to be considered in the calculation
    # Restated, False or 0 means "NOT masked out"
    KPathogenic_p  = biKest(maskPathogenic_p, nPathogenic, T)
    KNeutral_p  = biKest(maskNeutral_p, nNeutral, T)
    KPathogenicLessKNeutral_list.append(KPathogenic_p - KNeutral_p)

  # Permutation matrices
  KPathogenicLessKNeutral_matrices = np.array(KPathogenicLessKNeutral_list,dtype=np.float64)

  # Distance-dependent p-values and z-scores
  # The "two sided simulated p-value", DAB_zp, is not used further
  K_PlN_pvalues, K_PlN_zscores, DAB_zp,  K_PlN_high_confidence_envelope, K_PlN_low_confidence_envelope = pstats(KPathogenicLessKNeutral_matrices)

  # Plot and save the plot to graphics fle.
  saveDplot(T,KPathogenicLessKNeutral,K_PlN_zscores,K_PlN_low_confidence_envelope,K_PlN_high_confidence_envelope,label)

  # nanargmax returns the index where the maxMagnutdeKZ is found
  maxMagnitudeKzIndex = np.nanargmax(np.abs(K_PlN_zscores),axis=0)

  LOGGER.info("Ending biD() max (optimal) KPathogenicLessKNeutral max Zscore=%f found at distance threshold=%1.1f Angstroms",K_PlN_zscores[maxMagnitudeKzIndex],T[maxMagnitudeKzIndex])
  LOGGER.info("Also returning optimal KPathLessKNeutral[%1.1f A]=%.2f KPathLessKNeutral_pvalue[%1.1f A]=%1.2e",T[maxMagnitudeKzIndex],KPathogenicLessKNeutral[maxMagnitudeKzIndex],T[maxMagnitudeKzIndex],K_PlN_pvalues[maxMagnitudeKzIndex])
  return T,KPathogenicLessKNeutral[maxMagnitudeKzIndex],K_PlN_zscores[maxMagnitudeKzIndex],T[maxMagnitudeKzIndex],K_PlN_pvalues[maxMagnitudeKzIndex]

def ripley(complex_df,permutations=99999):
  """ Handles Ripley's univariate K and bivariate D
      analyses using variant pathogenic and neutral
      labels and a constructed inter-residue distance matrix

      Parameters
      ----------
      complex_df: Structural dataframe contining x/y/z values and a dcode for each residue
      """
  LOGGER.info("Begin ripley()  permutations=%d"%permutations)
  # Distance matrix for all residues, regardless of labeling
  D = squareform(pdist(complex_df[['x','y','z']]))
  # Set the diagnoal of the distance matrix to nans
  np.fill_diagonal(D,np.nan)
  # old way to do above lineD[np.identity(D.shape[0],dtype=bool)] = np.nan
  # Label vectors for pathogenic and neutral variants
  p = (complex_df["dcode"]==1).astype(int)
  n = (complex_df["dcode"]==0).astype(int)
  pT,pK,pKz,pKt,pKp = uniK_calc_and_plot(D,p,permutations,"pathogenic")
  nT,nK,nKz,nKt,nKp = uniK_calc_and_plot(D,n,permutations,"neutral")

  LOGGER.info("Ending ripley().  Univariate K:")
  LOGGER.info(" Neutral K:")
  LOGGER.info("  Most significant T = %d Angstroms from threshold range[%2d .. %2d]",nKt,nT[0],nT[-1])
  LOGGER.info("    where K = %.2f"%nK)
  LOGGER.info("    where Z = %.2f"%nKz)
  LOGGER.info("    where p = %g"%nKp)
  LOGGER.info(" Pathogenic K:")
  LOGGER.info("  Most significant T = %d Angstroms from threshold range[%2d .. %2d]",pKt,pT[0],pT[-1])
  LOGGER.info("    where K = %.2f"%pK)
  LOGGER.info("    where Z = %.2f"%pKz)
  LOGGER.info("    where p = %g"%pKp)

  ## Bivariate D
  # Coordinate matrices for pathogenic and neutral variants
  pathogenic_coordinates = complex_df.loc[complex_df["dcode"]==1,['x','y','z']]
  neutral_coordinates = complex_df.loc[complex_df["dcode"]==0,['x','y','z']]

  T,KPathogenicLessKNeutral_atMaxKzIndex,\
  K_PathogenicLessKNeutral_zscore_atMaxKzIndex,\
  T_atMaxKzIndex,\
  K_PathogenicLessKNeutral_pvalue_atMaxKzIndex = biD(pathogenic_coordinates,neutral_coordinates,permutations,"pathogenic-neutral")

  LOGGER.info("Bivariate D:")
  LOGGER.info("  Most significant distance threshold, t, is: %d A from threshold range[%2d .. %2d]",T_atMaxKzIndex,T[0],T[-1])
  LOGGER.info("    where PathProxScore[%d A] = %.2f",T_atMaxKzIndex,KPathogenicLessKNeutral_atMaxKzIndex)
  LOGGER.info("    where Zscore[%d A] of PathProxScore[%d A] = %.2f",T_atMaxKzIndex,T_atMaxKzIndex,K_PathogenicLessKNeutral_zscore_atMaxKzIndex)
  LOGGER.info("    and p[%d A] of PathProxScore = %1.2e",T_atMaxKzIndex,K_PathogenicLessKNeutral_pvalue_atMaxKzIndex)

  return  pK,pKz,pKt,\
          nK,nKz,nKt,\
          KPathogenicLessKNeutral_atMaxKzIndex,K_PathogenicLessKNeutral_zscore_atMaxKzIndex,T_atMaxKzIndex

def move_to_outdir(outdir):
  # Create the output directory and determine output labels
  timestamp = strftime("%Y-%m-%d")
  if not outdir:
    outdir = "results/PathProx_%s"%(args.label)
    # Only timestamp if no output directory explicitly specified
    if not args.no_timestamp and timestamp not in outdir:
      outdir += "_%s"%timestamp
  LOGGER.info("Output directory has been updated to %s"%outdir)
  cwd = os.getcwd()
  makedirs_capra_lab(outdir,'move_to_outdir')

  # Write current version of this script to output directory
  shutil.copy(__file__.rstrip('c'),outdir)
  os.chdir(outdir)
  return outdir

# 2019 October 14 Pathprox3: These attr files (which dump the input variant sets) must not be output if no dataframe rows are applicable.
def write_input_variants_to_chimera_attributes(dcode_filtered_df: pd.DataFrame,
                    variant_type: str,
                    attribute_text: str = None,
                    third_column:str = None):

    if not attribute_text:
        attribute_text = variant_type

    LOGGER.info("%d variants of type %s to write as chimera attributes %s*%s.attr",len(dcode_filtered_df),variant_type,args.label,variant_type)
    if len(dcode_filtered_df) < 1: # Don't write anything if no rows in the dataframe
        return

    with open("%s_%s.attr"%(args.label,variant_type),'w') as fout:
        fout.write("attribute: %s\n"%attribute_text)
        fout.write("match mode: 1-to-1\n")
        fout.write("recipient: residues\n")
        for i,row in dcode_filtered_df.iterrows():
            fout.write("\t:%s.%s\t%.3f\n"%(format_resid(row["resid"]),row["chain"],1.0 if not third_column else row[third_column]))

    with open("%s_renum_%s.attr"%(args.label,variant_type),'w') as fout:
        fout.write("attribute: %s\n"%attribute_text)
        fout.write("match mode: 1-to-1\n")
        fout.write("recipient: residues\n")
        for i,row in dcode_filtered_df.iterrows():
            fout.write("\t:%s.%s\t%.3f\n"%(row["unp_pos"],row["chain"],1.0 if not third_column else row[third_column]))

#=============================================================================#
if __name__ == "__main__":
    # Determine the script directory
    # Which is the directory that pathprox.py is running in plus "/bin"
    # A number of ancillary scripts are in that directory
    # This script will use pathvis.py to assist with chimera output
    script_dir = "%s/bin"%os.path.dirname(os.path.abspath(sys.argv[0]))

    # Parse remaining Command Line Options #
    import os,sys,argparse,configparser
    # from pdbmap import PDBMapModel

    curdir = os.getcwd() # Note the current directory


    # cmdline_parser.add_argument("refseq",type=str,
    #                  help="NM_.... refseq transcript ID")

    cmdline_parser.add_argument("--fasta",type=str,
                      help="UniProt fasta file for the reference sequence")
    cmdline_parser.add_argument("--pathogenic",type=str,
                      help="User-defined set of pathogenic variants")
    cmdline_parser.add_argument("--pathogenic_label",type=str,
                      help="Label for User-defined pathogenic variants")
    cmdline_parser.add_argument("--neutral",type=str,
                      help="User-defined set of neutral variants")
    cmdline_parser.add_argument("--neutral_label",type=str,
                      help="Label for User-defined neutral variants")
    cmdline_parser.add_argument("--quantitative",type=str,
                      help="User-defined set of residue quantitative traits")
    cmdline_parser.add_argument("--qlabel",type=str,
                      help="Label to use for the quantitative trait")
    cmdline_parser.add_argument("--no-blosum",action='store_true',default=False,
                      help="Do not weight training data by substitution severity (BLOSUM100)")

  # Variant database flags
    cmdline_parser.add_argument("--add_1kg",action="store_true",default=False,
                      help="Supplement neutral variant set with 1000 Genomes missense variants")
    cmdline_parser.add_argument("--add_exac",action="store_true",default=False,
                      help="Supplement neutral variant set with ExAC missense variants")
    cmdline_parser.add_argument("--add_gnomad",action="store_true",default=False,
                      help="Supplement neutral variant set with gnomAD missense variants")
    cmdline_parser.add_argument("--add_gnomad38",action="store_true",default=False,
                      help="Supplement neutral variant set with gnomAD 2.1.1 GRCh38 missense variants")
    cmdline_parser.add_argument("--add_benign",action="store_true",default=False,
                      help="Supplement neutral variant set with ClinVar benign missense variants")
    cmdline_parser.add_argument("--add_pathogenic",action="store_true",default=False,
                      help="Supplement pathogenic variant set with ClinVar pathogenic missense variants")
    cmdline_parser.add_argument("--add_pathogenic38",action="store_true",default=False,
                      help="Supplement pathogenic variant set with ClinVar GRCh38 pathogenic missense variants")
    cmdline_parser.add_argument("--add_drug",action="store_true",default=False,
                      help="Supplement pathogenic variant set with ClinVar drug response")
    cmdline_parser.add_argument("--add_cosmic",action="store_true",default=False,
                      help="Supplement pathogenic variant set with COSMIC somatic missense variants")
    cmdline_parser.add_argument("--add_cosmic38",action="store_true",default=False,
                      help="Supplement pathogenic variant set with COSMIC38 somatic missense variants")
    cmdline_parser.add_argument("--add_tcga",action="store_true",default=False,
                      help="Supplement pathogenic variant set with TCGA somatic missense variants")

  # Filter parameters
    cmdline_parser.add_argument("--chain",type=str,
                      help="Limit the analysis to one particular chain")
    cmdline_parser.add_argument("--isoform",type=str,
                      help="Explicit declaration of the reference isoform (ENST)")
    cmdline_parser.add_argument("--unp",type=str,
                      help="Explicit declaration of the reference uniprot id")
    cmdline_parser.add_argument("--use-residues",type=str,
                      help="Specify a range of residues to analyze in the form: `-500` or `1-500` or `500-`")

  # Analysis parameters
    cmdline_parser.add_argument("--radius",type=str,default="NW",
                      help="PathProx radius options: {'K'=Use Univariate 'max Kz' distance threshold,'D' use Ripley's bivariate 'max Kz','NW',<static radius e.g. 10>}")
    cmdline_parser.add_argument("--nwlb",type=float,default=8.,
                      help="Lower bound on the NeighborWeight function.")
    cmdline_parser.add_argument("--nwub",type=float,default=24.,
                      help="Upper bound on the NeighborWeight function")
    cmdline_parser.add_argument("--ripley",action='store_true',default=False,
                      help="Perform univariate K and bivariate D analyses. True by default if '--radius=K' or '--radius=D'.")
    cmdline_parser.add_argument("--permutations",type=int,default=99999,
                      help="Number of permutations for Ripley's K/D analyses")

  # Output parameters
    cmdline_parser.add_argument("--label",type=str,default='',
                      help="Optional analysis label (overrides entity inference)")
    cmdline_parser.add_argument("--no-timestamp","-nt",action="store_true",default=False,
                      help="Disables output directory timestamping")
    cmdline_parser.add_argument("--overwrite",action="store_true",default=False,
                      help="Overwrite previous results. Otherwise, exit with error")

    cmdline_parser.add_argument("--sqlcache",type=str,required=False,
                      help="Directory in which to cache pickled sql query returns, for fast repeat queries")

    LOGGER.info("Command: %s"%' '.join(sys.argv))

    # cmdline_parser.add_argument("collaboration",type=str,help="Collaboration ID (ex. UDN)")
    cmdline_parser.add_argument("mutation",nargs='?',type=str,default='',help="HGVS mutation string (ex S540A)")

    args,args_remaining = cmdline_parser.parse_known_args() # args_remaining)
    # if not os.path.exists(cache_dir):
    #  os.makedirs(cache_dir)


    LOGGER.info("Command: %s"%' '.join(sys.argv))

    statusdir = os.path.join(args.outdir,"status")
    LOGGER.info("Job status directory: %s"%statusdir)
    makedirs_capra_lab(statusdir,"Main statusdir creation")

    # Remove any prior statusdir contents
    for the_file in os.listdir(statusdir):
      file_path = os.path.join(statusdir, the_file)
      try:
        if os.path.isfile(file_path):
          os.unlink(file_path)
      except Exception as e:
        msg = "Unable to delete file %s from status directory"%file_path
        LOGGER.exception(msg)
        sys_exit_failure(msg)

    def __info_update(info):
      new_info_filename = os.path.join(statusdir,"info.new")
      with open(new_info_filename,'w') as f:
        f.write(info + '\n')
      final_info_filename = os.path.join(statusdir,"info")
      os.rename(new_info_filename,final_info_filename)
      LOGGER.info("%s now contains: %s"%(final_info_filename,info))

    def sys_exit_failure(info):
      __info_update(info)
      new_progress_filename = os.path.join(statusdir,"progress.new")
      with open(new_progress_filename,'w') as f:
        f.write("%s: %s\n"%(__file__,inspect.currentframe().f_back.f_lineno))
      os.rename(new_progress_filename,'%s/progress'%statusdir)
      # Mark this job as failed
      fail_filename = os.path.join(statusdir,"FAILED")
      open(fail_filename,'w').close()
      LOGGER.critical("Creating FAILURE file %s"%fail_filename)
      sys.exit(info)

    def statusdir_info(info):
      __info_update(info)
      new_progress_filename = os.path.join(statusdir,"progress.new")
      with open(new_progress_filename,'w') as f:
        f.write("%s: %s\n"%(__file__,inspect.currentframe().f_back.f_lineno))
      os.rename(new_progress_filename,'%s/progress'%statusdir)
      return info

    statusdir_info('Begun')

    required_config_items = ["dbhost","dbname","dbuser","dbpass",
      "collaboration",
      "pdb_dir",
      "sec2prim",
      "chimera_headless",
      "collaboration",
      "idmapping",
      "interpro_dir",
      "modbase2020_dir",
      "modbase2020_summary",
      # "modbase2013_dir",
      # "modbase2013_summary",
      # "modbase2016_dir",
      # "modbase2016_summary",
      "gnomad_dir",
      "gnomad_filename_template",
      "output_rootdir",
      "sprot",
      "swiss_dir",
      "swiss_summary",
      "vep","vep_cache_dir","vep_assembly"]

    config,config_dict = psb_config.read_config_files(args,required_config_items)

    config_dict_reduced = {x:config_dict[x] for x in required_config_items}

    config_dict = config_dict_reduced

    PDBMapGlobals.config = config_dict

    config_dict_shroud_password = {x:config_dict[x] for x in required_config_items}
    dbpass = config_dict.get('dbpass','?')
    config_dict_shroud_password['dbpass'] = '*' * len(dbpass)

    # pathprox_config_dict = dict(config.items('PathProx'))


    # LOGGER.info("Command Line Arguments")
    # pprint.pprint(vars(args))
    LOGGER.info("Command Line Arguments:\n%s"%pprint.pformat(vars(args)))
    LOGGER.info("Configuration File parameters:\n%s"%pprint.pformat(config_dict_shroud_password))
    # LOGGER.info("PathProx parameters:\n%s"%pprint.pformat(pathprox_config_dict))

    ## if args.sqlcache:
    ##   cache_dir = args.sqlcache
    ## else:
    ##   cache_dir = os.path.join(config_dict['output_rootdir'],config_dict['collaboration'],"sqlcache")

    ## # We must take care with the sql_cache because several processes can simultaneously attempt to make it!
    ## if not os.path.exists(cache_dir):
    ##   try:
    ##     os.makedirs(cache_dir)
    ##   except:
    ##     pass
    ##   try:
    ##     assert os.path.exists(cache_dir)
    ##   except:
    ##     sys_exit_failure("Fatal: Can't create sqlcache at destination path %s" %cache_dir)
    ## set_capra_group_sticky(cache_dir)

    ## LOGGER.info("SQL cache directory: %s"%cache_dir)

    structure = None
    chain_to_alignment = {}
    chain_to_transcript = {}

    # Step 1 is to parse out the assignments of chains to ids
    PDBMapProtein.load_idmapping(config_dict['idmapping'])

    # PROCESS USER-SPECIFIED transcript<->chain mappings
    # if the user has any transcript->chain assignments, parse them from command line
    chain_to_transcript = {}
    transcript_parser = argparse.ArgumentParser()
    # Chain letters are usually just one letter - but can be two
    # So allow for up to three
    chain_X_unp_re = re.compile('^--chain(.{1,3})unp=(.*)$')
    chain_X_enst_re = re.compile('^--chain(.{1,3})enst=(.*)$')
    chain_X_fasta_re = re.compile('^--chain(.{1,3})fasta=(.*)$')
    transcript = None
    for arg in args_remaining:
        re_match =  chain_X_unp_re.match(arg)
        unp = None
        enst = None
        fasta = None
        chain_letter = None

        if re_match: # A uniprot ID transcript ID was explicitly assigned to the chain
            chain_letter = re_match.group(1)
            unp = re_match.group(2)
            transcript = PDBMapTranscriptUniprot(unp)
            # if chain_letter in chain_to_transcript and transcript.id == chain_to_transcript[chain_letter].id:
            #    LOGGER.warning("argument %s is redundant.  Pathprox automatically assigned/aligned")
            LOGGER.info("Successful load of Uniprot transcript %s"%unp)
        else:
            re_match =  chain_X_enst_re.match(arg)
            if re_match: # An ENSEMBL transcript ID was explicitly assigned to the chain
                chain_letter = re_match.group(1)
                enst = re_match.group(2)
                transcript = PDBMapTranscriptEnsembl(enst)
                LOGGER.info("Successful load of Ensembl transcript %s"%enst)
            else:
                re_match =  chain_X_fasta_re.match(arg)
                if re_match: # A fasta amino acid string was explicitly assigned to the chain
                    chain_letter = re_match.group(1)
                    fasta = re_match.group(2)
                    transcript = PDBMapTranscriptFasta(fasta)
                    LOGGER.info("Successful load of fasta transcript %s"%fasta)
                else:
                    exitmsg = "Command line argument %s is invalid"%arg
                    LOGGER.critical(exitmsg)
                sys_exit_failure(exitmsg)

        if transcript:
            """alignment = PDBMapAlignment()
            (success,message) = alignment.align_trivial(transcript,structure,chain_id=chain_letter)
            if success:
                LOGGER.info("%s to chain %s Trivial align successful for residues in [%d,%d]\n%s"%(transcript.id,chain_letter,next(iter(alignment.seq_to_resid)),next(reversed(alignment.seq_to_resid)),alignment.aln_str))
            else:
                LOGGER.warning("Unable to trivially align transcript %s to %s.%s"%(transcript.id,structure.id,chain_letter))
                if args.pdb or args.biounit and unp: # Then user is mapping a deposited chain to a unp
                   if PDBMapProtein.isCanonicalByUniparc(transcript.id):
                       success,errormsg = alignment.align_sifts_canonical(transcript,structure,chain_id=chain_letter)
                   else:
                       success,errormsg = alignment.align_sifts_isoform_specific(transcript,structure,chain_id=chain_letter)

                if success:
                    LOGGER.info("Sifts alignment identity %f.  Align str:\n%s"%(alignment.perc_identity,alignment.aln_str))
                else:
                    (success,message) = alignment.align_biopython(transcript,structure,chain_id=chain_letter)
                    assert success,statusdir_info(message)
                    LOGGER.info("Biopython identity %f.  Align str:\n%s"%(alignment.perc_identity,alignment.aln_str))"""

            # chain_to_alignment[chain_letter] = alignment
            chain_to_transcript[chain_letter] = transcript

    is_biounit = 0  # Flag only set to 1 _If_ it is a PDB we load and _IF_ found in biounits
    if args.pdb: # Then this is going to be a single chain calculation - so do try a biounit which may lack chain of interest
        structure,is_biounit,chain_to_transcript,chain_to_alignment = PDBMapComplex_load_pdb_align_chains(args.pdb,False,chain_to_transcript)
    if args.biounit:
        is_biounit = 1
        structure,is_biounit,chain_to_transcript,chain_to_alignment = PDBMapComplex_load_pdb_align_chains(args.biounit,is_biounit,chain_to_transcript)
    elif args.usermodel:
        if not args.label:
            exitmsg = "--label is required on the command line when loading a --usermodel"
            LOGGER.critical(exitmsg)
            sys_exit_failure(exitmsg)
        structure = load_structure(args.label,args.usermodel)
    elif args.swiss: # This is a lenghty-ish swiss model ID, NOT the file location
        PDBMapSwiss.load_swiss_INDEX_JSON(config_dict['swiss_dir'],config_dict['swiss_summary']);
        structure = load_structure(args.swiss,PDBMapSwiss.get_coord_file(args.swiss))

        assigned_chain_id = next(iter(chain_to_transcript))

        # Swiss _could_ be a homo-oliomer - and we want to capture that!
        # by pointing the additional chain IDs at the same alignment
        for chain in structure[0]:
            if chain.id not in chain_to_transcript:
                first_residue = next(iter(structure[0][chain.id]))
                # Swiss model seems to put some HETATMs in some of the chains.  DO NOT align to those!
                if first_residue.id[0] == ' ':
                    chain_to_transcript[chain.id] = chain_to_transcript[assigned_chain_id]
            
    elif args.modbase: # This is a ENSP.... modbase file.  Currently we look only in modbase20
        modbase20 = PDBMapModbase2020(config_dict)
        modbase20_structure_filename = modbase20.get_coord_file(args.modbase)
        if os.path.exists(modbase20_structure_filename):
           structure = load_structure(args.modbase,modbase20_structure_filename)
        else:
           msg = "Modbase model %s not found in %s"%(args.modbase,modbase20_structure_Filename)
           LOGGER.critical(msg)
           sys.exit(msg)
        # modbase20 = PDBMapModbase2020(config_dict)
        # modbase16_structure_filename = modbase16.get_coord_file(args.modbase)
        # if os.path.exists(modbase16_structure_filename):
        #    structure = load_structure(args.modbase,modbase16_structure_filename)
        # else:
        #    modbase13 = PDBMapModbase2013(config_dict)
        #    modbase13_structure_filename = modbase13.get_coord_file(args.modbase)
        #    if os.path.exists(modbase13_structure_filename):
        #        structure = load_structure(args.modbase,modbase13_structure_filename)
        #    else:
        #        LOGGER.critical("Modbase model id %s not found in either modbase2013/2016 config directories");

    # Finally - if the chain left has no id, set the id to A for sanity
    for chain in list(structure.get_chains()):
        if chain.id in ('',' '):
            LOGGER.info('Renaming blank/missing chain ID to A')
            chain.id = 'A'

    if args.usermodel or args.swiss or args.modbase:
        alignment = PDBMapAlignment()
        for chain_letter in chain_to_transcript:
            alignment = PDBMapAlignment()
            (success,message) = alignment.align_trivial(chain_to_transcript[chain_letter],structure,chain_id=chain_letter)
            if success:
                LOGGER.info("%s to chain %s Trivial align successful for residues in [%d,%d]\n%s"%(transcript.id,chain_letter,next(iter(alignment.seq_to_resid)),next(reversed(alignment.seq_to_resid)),alignment.aln_str))
            elif args.usermodel:
                (success,message) = alignment.align_biopython(chain_to_transcript[chain_letter],structure,chain_id=chain_letter)
                if success:
                    LOGGER.warning("TRIVIAL ALIGNMENT OF USER MODEL FAILED\n%s to chain %s Biopython successful for residues in [%d,%d]\n%s"%(transcript.id,chain_letter,next(iter(alignment.seq_to_resid)),next(reversed(alignment.seq_to_resid)),alignment.aln_str))

            assert success,statusdir_info(message)
            chain_to_alignment[chain_letter] = alignment


    assert structure,statusdir_info("A structure file must be specified via --pdb, --biounit, --swiss, --modbase, or --usermodel")
    # Preprocess structural properties
    # Pathprox cannot use models beyond the 0th one - so drop those right away
    models = list(structure.get_models())
    if len(models) > 1: 
        first_model = True;
        for model in models:
            if not first_model:
                LOGGER.info("Removing extraneous model %d from %s",model.id,structure.id)
                structure.detach_child(model.id)
            first_model = False

    model = models[0]
    if args.chain:    # User wants analysis restricted to one chain
        if not args.chain in model:
            msg = "You requested restriction to chain %s.  But, this is not in structure %s"%(args.chain,structure.id)
            LOGGER.critical(msg)
            sys_exit_failure(msg)
        for chain in list(model.get_chains()):
            if chain.id != args.chain:
                LOGGER.info('Dropping chain %s from structure %s (args.chain=%s)',chain.id,structure.id,args.chain)
                model.detach_child(chain.id)

    # Finally - if the chain left has no id, set the id to A for sanity
    chain_count = len(list(model.get_chains()))

    residue_coms = {}

    for chain in model:
        for residue in chain:
                    # Compute the center of mass for all residues
            natom = float(len([atom for atom in residue]))
            com = sum([atom.coord for atom in residue]) / natom
            residue_coms[(chain.id,residue.id)] = com

    for chain in structure[0]:
        if chain.id not in chain_to_alignment:
            LOGGER.warning("*** %s.%s has not been aligned to a transcript.",structure.id,chain.id)

    renumbered_structure = structure_renumber_per_alignments(structure)

    # In a homo-dimer, multiple chains will associate with the same chain.  Create that reverse structure
    transcript_to_chains = {}
    for chain_letter in chain_to_transcript:
        chains = transcript_to_chains.get(chain_to_transcript[chain_letter].id,[])
        chains.append(chain_letter)
        transcript_to_chains[chain_to_transcript[chain_letter].id] = chains
        


    LOGGER.info("Remaining arguments are %s"%args_remaining)


    ## Parameter validity checks
    if not args.label:
      # Generate label from input parameters
      # Strip file extensions
      args.label = structure.id

    # Append relevant information to the label
    # if chain_to_alignment:
    #    args.label += "_%s"%args.chain ''.join(chain_to_alignment.keys())
    if args.chain:
      args.label += "_%s"%args.chain
    elif args.variants and len(args.variants.split(':')) > 1:
      args.label += "_%s"%args.variants.split(':')[0]
    if args.add_pathogenic:
      args.label += "_clinvar"
    if args.add_pathogenic38:
      args.label += "_clinvar38"
    if args.add_cosmic:
      args.label += "_cosmic"
    if args.add_cosmic38:
      args.label += "_cosmic38"
    if args.add_tcga:
      args.label += "_tcga"
    if args.pathogenic:
      if args.pathogenic_label:
        args.label += '_%s'%args.pathogenic_label
      else:
        args.label += '_pathogenic'
    if args.add_exac:
      args.label += "_exac"
    if args.add_gnomad:
      args.label += "_gnomad"
    if args.add_gnomad38:
      args.label += "_gnomad38"
    if args.neutral:
      if args.neutral_label:
        args.label += '_%s'%args.neutral_label
      else:
        args.label += '_neutral'
    if args.add_1kg:
      args.label += "_1kg"
    # Add the radius type and bounds if NeighborWeight
    if args.radius == "NW":
      args.label += "_NW-%.0f-%.0f"%(args.nwlb,args.nwub)
    elif args.radius == "K":
      args.label += "_K"
    elif args.radius == "D":
      args.label += "_D"
    else:
      args.radius = float(args.radius)
      args.label += "_%.1f"%args.radius

    LOGGER.info("Using analysis label: %s"%args.label)

    if args.use_residues:
      # Limit analysis to user-selected protein residues
      try:
        args.use_residues = tuple(args.use_residues.split('-'))
      except:
        msg = "Incorrect formatting for --use-residues. See help for details.\n"
        LOGGER.critical(msg)
        sys_exit_failure(msg)
    if args.radius not in ("K","D","NW"):
      # Check that a valid numeric radius was provided
      try:
        args.radius = float(args.radius)
      except:
        msg = "Invalid option for --radius. See help for details.\n"
        LOGGER.critical(msg)
        sys_exit_failure(msg)
    elif args.radius in ("K","D"):
      # Ripley's K/D analyses required for parameterization
      LOGGER.info("Setting args.ripley=True because args.radius=%s",args.radius)
      args.ripley = True

    if args.verbose:
      LOGGER.info("Active options:")
      for arg in vars(args):
        try:
          LOGGER.info("  %15s:  %s"%(arg,getattr(args,arg).name))
        except:
          LOGGER.info("  %15s:  %s"%(arg,getattr(args,arg)))
      LOGGER.info("")


    AAseq = None
    # Let this one go for now
    if False and args.isoform and not args.fasta:
      # Query the Ensembl API for the transcript
      # pdbmap will be in the PATH from psb_prep.bash
      cmd = "transcript_to_AAseq.pl %s"%args.isoform
      LOGGER.info("Executing: %s"%cmd)
      status, stdout_stderr =  subprocess.getstatusoutput(cmd)
      if status > 0:
        LOGGER.critical("Exit Code of %d from perl %s stdout:\n%s"%(status,cmd,stdout_stderr))
        sys.exit(1)
      else:
        AAseq = stdout_stderr.rstrip()
        LOGGER.info("%s -> %s",args.isoform,AAseq);
 

    # Warn user to include specific EnsEMBL transcripts
    if args.fasta and not args.isoform and \
            (args.add_gnomad or args.add_gnomad38 or \
             args.add_pathogenic or args.add_pathogenic38 or \
             args.add_cosmic or args.add_cosmic38 or \
             args.add_exac or args.add_1kg or args.add_benign or \
             args.add_drug or args.add_tcga):
      msg  = "\n!!!!===========================================!!!!\n"
      msg += "                        WARNING\n"
      msg += "   Reference EnsEMBL transcript was not specified. \n"
      msg += "    Are there multiple isoforms for this protein?\n"
      msg += "  Explictly declare isoform to avoid mis-alignments\n\n"
      msg += "!!!!===========================================!!!!\n\n"
      LOGGER.warning(msg)

    statusdir_info('Configured')
    #=============================================================================#
    ## Begin Analysis ##

    # Init PDBMap by load the summary information dictionaries that it will need to use:
    LOGGER.info("Initializing PDBMAP by loading idmapping sec2prom, sprot, modbase, and swiss model meta dictionaries")
    PDBMapProtein.load_sec2prim(config_dict['sec2prim'])
    PDBMapProtein.load_sprot(config_dict['sprot'])
    statusdir_info('Initialized')

    # Reinitialize variant sets for each structure
    quantitative_variant_weights = parse_qt(args.quantitative)
    pathogenic_variant_sets = PDBMapVariantSet.parse_variants("Pathogenic",args.pathogenic)
    neutral_variant_sets = PDBMapVariantSet.parse_variants("Neutral",args.neutral)
    candidate_variant_sets = PDBMapVariantSet.parse_variants("Candidate",args.variants)

    def log_variants(variant_source: str):
        LOGGER.info(variant_source)
        def log_variants_info(variants_flavor: str,variant_sets: Dict[str, PDBMapVariantSet]):
            LOGGER.info(" %-13.13s %d in %d transcripts"%(variants_flavor+':',sum(len(variant_sets[id].variant_list) for id in variant_sets),len(variant_sets)))
        def log_variants_details(variants_flavor: str,variant_sets: Dict[str, PDBMapVariantSet]):
            all_variants_described = ""
            for id in variant_sets:
                all_variants_described += "   %s: %d %s variants:\n"%(id,len(variant_sets[id].variant_list),variants_flavor)
                for variant in variant_sets[id].variant_list:
                    all_variants_described += "%s\n"%str(variant)
            LOGGER.info(" %-13.13s %d in %d transcripts\n%s"%(variants_flavor+':',
                    sum(len(variant_sets[id].variant_list) for id in variant_sets),len(variant_sets),all_variants_described))
        log_variants_info('Neutral',neutral_variant_sets)
        log_variants_details('Neutral',neutral_variant_sets)
        log_variants_info('Pathogenic',pathogenic_variant_sets)
        log_variants_info('Candidates',candidate_variant_sets)
        log_variants_info('Quantitative',quantitative_variant_weights)

    log_variants("User supplied variant counts:")

    LOGGER.info("SQL-Querying additional variants if requested for %s"%[(chain_id,chain_to_transcript[chain_id].id) for chain_id in chain_to_transcript])

    # Supplement with any requested variant datasets from SQL
    # 2019 October - loop over chain_to_transcript
    # convert to lists of ENST
    # load variants in a per-transcript manner
    def chain_to_ENST_transcripts(chain_id):
        """Given a chain ID from the structure, return a list of ensembl transcripts which 
           can be used to query genomic databases for known variants.
           
           Requires dict chain_to_transcript to map chain_id to a transcript """

        if chain_id not in chain_to_transcript:
            LOGGER.critical("Chain %s is not associated with any transcript ID.  Variants cannot be loaded from SQL for chain %s"%(chain_id,chain_id))
            return []

        transcript = chain_to_transcript[chain_id]

        # Uniprot identifiers are needed to query the PDBMapProtein idmapping file
        unps = [] # Empty list of uniprot identifiers

        ENST_transcripts = []
        if type(transcript) == PDBMapTranscriptUniprot:
            unps = [transcript.id]
        elif type(transcript) == PDBMapTranscriptEnsembl:
            ENST_transcripts = [transcript]
            unps = PDBMapProtein.enst2unp(transcript.id)
        else:
            unps = []


        for unp in unps:
            ensemble_transcript_ids = PDBMapProtein.unp2enst(unp)
            for ensembl_transcript_id in ensemble_transcript_ids:
                ensembl_transcript = PDBMapTranscriptEnsembl(ensembl_transcript_id)
                if not ensembl_transcript.aa_seq:
                    LOGGER.warning("Ensembl transcript %s has no associated aa_seq.  Skipping"%ensembl_transcript_id)
                    continue
                assert transcript.aa_seq == ensembl_transcript.aa_seq, statusdir_info("%s and %s AA sequences differ:\n%s"%(
                    transcript.id,ensembl_transcript.id,PDBMapTranscriptBase.describe_transcript_differences(transcript,ensembl_transcript)))
                LOGGER.info("Ensembl transcript %s has same aa_seq as transcript %s",ensembl_transcript.id,transcript.id)
                if ensembl_transcript not in ENST_transcripts:
                    ENST_transcripts.append(ensembl_transcript)

        if not ENST_transcripts:
            LOGGER.critical("No ENSMBL transcripts were referenced to chain %s transcripts %s",chain_id,str(unps))
        return ENST_transcripts

    for chain in structure[0]:
        # import pdb; pdb.set_trace()
        if chain.id not in chain_to_alignment:
            LOGGER.critical("Chain %s has not been aligned to a transcript.  Skipping SQL variant load",chain.id)
            continue
        ENST_transcripts = chain_to_ENST_transcripts(chain.id)
        if not ENST_transcripts:
            LOGGER.critical("Skipping SQL variant load for chain %s as no ENST transcripts were mapped",chain.id)
            continue

        # Typically, we want the variants associated to a unp identifier
        variant_set_id = chain_to_transcript[chain.id].id 

        if args.add_1kg:
            PDBMapVariantSet.query_and_extend('Neutral',neutral_variant_sets,variant_set_id,ENST_transcripts,'1kg')
        if args.add_exac:
            PDBMapVariantSet.query_and_extend('Neutral',neutral_variant_sets,variant_set_id,ENST_transcripts,'exac')
        if args.add_gnomad:
            PDBMapVariantSet.query_and_extend('Neutral',neutral_variant_sets,variant_set_id,ENST_transcripts,'gnomad')
        if args.add_gnomad38:
            extended_count = PDBMapVariantSet.query_and_extend('Neutral',neutral_variant_sets,variant_set_id,ENST_transcripts,'gnomad38')


            if extended_count == 0: # This should be modularized - chill for now.  Test code
                _variant_set = neutral_variant_sets.get(variant_set_id,PDBMapVariantSet('Neutral'))
                extended_variant_count = 0
                prior_variant_count = len(_variant_set._variants)

                # See if we can get some variants from the source Gnomad vcf files
                pdbmap_gnomad = PDBMapGnomad(config_dict)
                for enst_transcript in ENST_transcripts:
                    vep_echo_filename = os.path.join(args.outdir,"%s_vep_results.vcf"%enst_transcript.id)
                    LOGGER.info("SQL unsuccessful.  Retrieving Gnomad variants from source .vcf files.  VEP output echoed to %s"%vep_echo_filename)
                    gnomad_variant_df = pdbmap_gnomad.retrieve_gnomad_missense_variants(enst_transcript,output_directory=args.outdir)

                    for index,row in gnomad_variant_df.iterrows():
                        if row['maf'] >= 1E-5: # Sorry to hard code - but this is a gnomad thing
                            _variant_set._variants.append([
                               int(row['Protein_position']),
                               row['Ref_AminoAcid'],
                               row['Alt_AminoAcid'],
                               None]) # <- force_chain None correct
                extended_variant_count = len(_variant_set._variants)

                if variant_set_id in neutral_variant_sets:
                    pass  # Because variant_set is the reference already in the dictionary from before
                else: # It's a new PDBMapVariantSet() instance, so add to dictionary
                    neutral_variant_sets[variant_set_id] = _variant_set
                LOGGER.info("%d variants obtained from Gnomad source .vcf files"%extended_variant_count)
        if args.add_benign:
            PDBMapVariantSet.query_and_extend('Neutral',neutral_variant_sets,variant_set_id,ENST_transcripts,'clinvar')

        if args.add_pathogenic:
            PDBMapVariantSet.query_and_extend('Pathogenic',pathogenic_variant_sets,variant_set_id,ENST_transcripts,'clinvar')
        if args.add_pathogenic38:
            PDBMapVariantSet.query_and_extend('Pathogenic',pathogenic_variant_sets,variant_set_id,ENST_transcripts,'clinvar38')
        if args.add_drug:
            PDBMapVariantSet.query_and_extend('Pathogenic',pathogenic_variant_sets,variant_set_id,ENST_transcripts,'drug')
        if args.add_cosmic:
            PDBMapVariantSet.query_and_extend('Pathogenic',pathogenic_variant_sets,variant_set_id,ENST_transcripts,'cosmic')
        if args.add_cosmic38:
            PDBMapVariantSet.query_and_extend('Pathogenic',pathogenic_variant_sets,variant_set_id,ENST_transcripts,'cosmic38')
        if args.add_tcga:
            PDBMapVariantSet.query_and_extend('Pathogenic',pathogenic_variant_sets,variant_set_id,ENST_transcripts,'tcga')

    log_variants("Final supplied variant counts:")

    # Create (and move to) output directory
    args.outdir = move_to_outdir(args.outdir)

    # We now gather all the variants, for all transcripts, for all chains, into a single dataframe with 6 columns
    variants_DataFrame = pd.DataFrame()
    for (dcode,variant_set) in [(0,neutral_variant_sets), (1,pathogenic_variant_sets), (-1,candidate_variant_sets), (np.NaN,quantitative_variant_weights)]:
        for transcript_id in variant_set: # Careful, transcript_id could be None 
            df = pd.DataFrame()
            if not transcript_id:
                df = variant_set[transcript_id].to_DataFrame()
                # df['chain'] _must_ come from a chain-override in the user input unless this is a single chain structure
                if chain_count == 1 and len(df) and (df.iloc[0]['chain'] == None or df.iloc[0]['chain'] == np.NaN):
                    df['chain'] = list(structure[0].get_chains())[0].id
            else: # Gather up _all_ the chains associated with the transcript ID and map variants to all fo them
                if transcript_id not in transcript_to_chains:
                    if dcode == -1 and len(transcript_id) < 6: # Then the transcript is really a specific chain_id override on the variant
                        chain_id = transcript_id
                        new_df = variant_set[transcript_id].to_DataFrame(chain_id=chain_id)
                        new_df['chain'] = chain_id
                        df = pd.concat([df,new_df],ignore_index=True)
                    else:
                        LOGGER.warning("Transcript %s found in variant inputs, but is not in structure.  Skipping: %s",transcript_id,str(variant_set[transcript_id]))
                        continue # Nothing to concat to our frame
                else:  # transcript_id is in transcript_to_chains
                    for chain_id in transcript_to_chains[transcript_id]:
                        new_df = variant_set[transcript_id].to_DataFrame(chain_id=chain_id)
                        new_df['chain'] = chain_id
                        df = pd.concat([df,new_df],ignore_index=True)
                        if dcode == -1: # For now, we do NOT want to map a "candidate variant" twice just because the transcript is twice in complex
                            break
            df['dcode'] = dcode
            # The quantitative weights are not right - and must be fixed
            df['qt'] = np.NaN

            variants_DataFrame = pd.concat([variants_DataFrame,df],ignore_index=True)

    # With the variants now gathered into a dataframe with columns=["unp_pos","ref","alt","chain","dcode","qt"])
    # Remove all duplicate rows (indices are ignored)
    variants_DataFrame.drop_duplicates(inplace=True)

    # Create a new MonotonicRangeIndex from 0 to len()-1.  Old index will be discarded
    variants_DataFrame.reset_index(drop=True,inplace=True)

    # If we have no variants, then we have to exit.  But, give the user
    # some ideas about how to find variants in SQL or via text files
    msg = None
    if variants_DataFrame.empty:
        msg = "No variants loaded."
        dash_dash_add_arguments = [] # What --add_* arguments were on command line?
        vars_args = vars(args)
        for key in vars_args:
            if vars_args[key] and str(vars_args[key]).startswith("add_"):
                dash_dash_add_arguments.append("--"+key)
        if dash_dash_add_arguments:
            msg += " SQL queries " + str(dash_dash_add_arguments) + " found no variants."
        else:
            msg += " Consider --add_<dataset> arguments, or files of variants for this complex"

        LOGGER.critical(msg)
        sys_exit_failure(msg)

    # If a position has both pathogenic and neutral variants, THEN eliminate the neutral ones
    def defer_to_pathogenic(g):
        # If all sub-group elements are neutral or candidate VUSs, then 
        # simply return all of them.  Otherwise, return a group that has
        # excluded all neutrals, and is _only_ pathogenic or candidates 
        return g if all(g["dcode"]<=0) else g[g["dcode"]!=0]

    vdf_neutrals_deferred = variants_DataFrame.groupby(["unp_pos","ref","alt","chain"],as_index=False).apply(defer_to_pathogenic)

    len_before_drop_duplicates = len(vdf_neutrals_deferred)
    vdf_neutrals_deferred = vdf_neutrals_deferred.drop_duplicates(["unp_pos","ref","chain","dcode"])
    if len(vdf_neutrals_deferred) != len_before_drop_duplicates:
        LOGGER.info("After elimination of duplicate rows same [unp_pos,ref,chain,dcode] dataframe lost %d rows"%(len_before_drop_duplicates-len(vdf_neutrals_deferred)),)

    def resicom(resi):
       if resi==None: return [None,None,None]
       return np.array([np.array(a.get_coord()) for a in resi.get_unpacked_list()]).mean(axis=0)
    # Construct a dataframe from the coordinate file

    def unp_pos_if_aligned(chain_id,res_id):
        if chain_id not in chain_to_alignment:
            return None
        alignment = chain_to_alignment[chain_id]
        if res_id not in alignment.resid_to_seq:
            return None
        return int(alignment.resid_to_seq[res_id])

    # In one iteration through the entire structure, we create a ne dataframe consisting of
    # Each aligned residue's chain ID (r.get_parent()), the pdb residue ID tuple, the transcript position (unp_pos)
    # and the x/y/z center of mass of the residue
    # The "for r..." list comprehension includes only non-hetero residues with (blank,#,insert code) i.e. (not r.id(0).strip())
    # and only then of they are aligned to transcript positions i.e. unp_pos_if_aligned... != None
    structure_centroids_df = pd.DataFrame([[r.get_parent().id,r.id,unp_pos_if_aligned(r.get_parent().id,r.id)]+ list(resicom(r)) 
        for r in structure.get_residues() if (
           (not r.id[0].strip()) and (unp_pos_if_aligned(r.get_parent().id,r.id))
           ) 
        ],columns=["chain","resid","unp_pos","x","y","z"])

   

    # Merge the sequence-aligned structure dataframe with the sequence-based variant dataframe
    vdf_extract = vdf_neutrals_deferred[["unp_pos","ref","alt","chain","dcode","qt"]].copy()
    vdf_extract['unp_pos'] = vdf_extract['unp_pos'].apply(int)
    vdf_extract['chain'] = vdf_extract['chain'].apply(str)
    structure_centroids_df['unp_pos'] = structure_centroids_df['unp_pos'].apply(int)
    structure_centroids_df['chain'] = structure_centroids_df['chain'].apply(str)
    complex_df_merged = structure_centroids_df.merge(vdf_extract,how="left",on=["chain","unp_pos"])
    # Ensure that no duplicate residues were introduced by the merge
    complex_df_merged.reset_index(drop=True,inplace=True)
    complex_df = complex_df_merged.drop_duplicates(["unp_pos","resid","chain","ref","alt","dcode"]).reset_index(drop=True)
  
    # Annotate each amino acid substitution with the blosum100 score
    complex_df["blosum100"] = complex_df.apply(lambda x: blosum100[(x["ref"],x["alt"])] if not np.isnan(x["dcode"]) else None,axis=1)
  
    # Fail if all variants outside structural coverage or mapped to missing residues
    msg = None
    vdf = complex_df[~complex_df["dcode"].isnull()]
    if vdf.empty:
      if not complex_df.empty:
        msg = "Structure successfully aligned to sequence, but no variants were mapped to non-missing residues."
      else:
        msg = "ERROR: Sequence-structure alignment failed."
      LOGGER.critical(msg); 
      sys_exit_failure(msg)

    # Check if there are enough variants to calculate neutral/pathogenic constraint
    sufficient_neutral_variants = (complex_df["dcode"]==0).sum() > 2
    sufficient_pathogenic_variants = (complex_df["dcode"]==1).sum() > 2

    pClosestPathogenicMutationPDB = None
    pNextPathogenicMutationPDB = None
    pClosestPathogenicMutationUNP = None
    pNextPathogenicMutationUNP = None
    pClosestPathogenicDistance = None
    pNextPathogenicSeqDelta = None

    mapped_variant_count =  (complex_df["dcode"] == -1).sum()

    def format_resid(resid):
         # If that residue type first element of the biopython has anything, then
         # something weird is going on, so make some noise
         if len(resid[0].strip()):
           LOGGER.warning("Malformed residue tuple: %s   Reformatting for readability/chimera.",str(resid))

         # If there is an insertion code then return 123F for residue id
         if len(resid[2].strip()):
            resid_formatted = "%d%s"%(resid[1],resid[2])
            LOGGER.warning("Residue %s has insertion code and will be reformatted as %s"%(resid,resid_formatted))
            return resid_formatted

         # The usual case, a pdb resdidue number with out a type, without an insertion code
         return str(resid[1])
       
    # Check that candidate variants did not map to missing residues
    if args.variants and mapped_variant_count < 1:
        LOGGER.info("All candidate variants were aligned to missing residues.")
        # Unset the variants argument for downstream processing
        args.variants = []
    elif mapped_variant_count != 1: # Can't easily talk about "nearest" residues if more than 1 unknwn
        LOGGER.info("%d unknown variants were mapped to the structure.  Skipping calculation of nearest pathogenic residues."%mapped_variant_count)
    else: # mapped_variant_count == 1: # One variant is normal case for UDN
        df_variant = complex_df[complex_df["dcode"] == -1]      
        df_pathogenic = complex_df[complex_df["dcode"] == 1]      
        if len(df_pathogenic) < 1:
            LOGGER.warning("No pathogenic mutations for calculation of nearest-to-variant.  Skipping.")
        else:
            #First, get the distance in 3-space
            variantCOM = df_variant[['x','y','z']].values
            pathogenicCOMs = df_pathogenic[['x','y','z']].values
            distances_to_pathogenic = cdist(variantCOM,pathogenicCOMs) # Calculate distance matrix 
            # Get nearest pathogenic position
            min_subscript = np.argmin(distances_to_pathogenic[0])
            pClosestPathogenicDistance = distances_to_pathogenic[0][min_subscript]
  
            df_pathogenic_nearest = df_pathogenic.iloc[min_subscript]
            pClosestPathogenicMutationPDB = "%s%s%s"%(
              df_pathogenic_nearest['ref'],format_resid(df_pathogenic_nearest['resid']),df_pathogenic_nearest['alt'])
            pClosestPathogenicMutationUNP = "%s%d%s"%(
              df_pathogenic_nearest['ref'],int(df_pathogenic_nearest['unp_pos']),df_pathogenic_nearest['alt'])
   
            LOGGER.info("3-space Nearest to is %d %s %s"%(min_subscript,pClosestPathogenicMutationPDB,pClosestPathogenicMutationUNP))

            #Second, get the distance in sequence space
            variant_unp_pos = df_variant[['unp_pos']].values
            pathogenic_unp_poss = df_pathogenic[['unp_pos']].values
            distances_to_pathogenic = cdist(variant_unp_pos,pathogenic_unp_poss) # Distance Matrix
            pNextPathogenicSeqDelta = int(distances_to_pathogenic[0][min_subscript])
            # Get nearest pathogenic position
            min_subscript = np.argmin(distances_to_pathogenic[0])
            df_pathogenic_nearest = df_pathogenic.iloc[min_subscript]
            pNextPathogenicMutationPDB = "%s%s%s"%(
              df_pathogenic_nearest['ref'],format_resid(df_pathogenic_nearest['resid']),df_pathogenic_nearest['alt'])
            pNextPathogenicMutationUNP = "%s%d%s"%(
              df_pathogenic_nearest['ref'],int(df_pathogenic_nearest['unp_pos']),df_pathogenic_nearest['alt'])
            print("Sequence Nearest is %d %s %s"%(min_subscript,pNextPathogenicMutationPDB,pNextPathogenicMutationUNP))


    basePDBfilename = "%s_%s_%s"%(args.label,structure.id,is_biounit)
    renumberedPDBfilename = "%s_renum.pdb"%basePDBfilename         
    renumberedCIFfilename = "%s_renum.cif"%basePDBfilename         
    renumberedPDBfilenameSS = "%s_renumSS.pdb"%basePDBfilename # Renumbered + chimera HELIX/SHEET notations
    renumberedCIFfilenameSS = "%s_renumSS.cif"%basePDBfilename         
    originalPDBfilename = "%s.pdb"%basePDBfilename
    originalCIFfilename = "%s.cif"%basePDBfilename

    # 
    # Write the renumbered PDB to the output directory
    # 
    LOGGER.info("Writing renumbered CIF to file %s"%renumberedPDBfilename)
    cifio = MMCIFIO()
    cifio.set_structure(renumbered_structure)
    cifio.save(renumberedCIFfilename)
    LOGGER.info("Writing renumbered PDB to file %s"%renumberedPDBfilename)
    pdbio = PDBIO()
    pdbio.set_structure(renumbered_structure)
    pdbio.save(renumberedPDBfilename)

    # We now have a "clean" renumbered PDB - but for improved ngl visualization, helps to
    # let chimera add some secondary structure annotation
    script = '"%s/writePdbWithSecondary.py %s %s"'%(script_dir,renumberedPDBfilename,renumberedPDBfilenameSS)
    cmd    = "TEMP=$PYTHONPATH; unset PYTHONPATH; %s --nogui --silent --script %s; export PYTHONPATH=$TEMP"%(config_dict['chimera_headless'],script)
    # Allow Mac OSX to use the GUI window
    if platform.system() == 'Darwin':
        cmd  = "TEMP=$PYTHONPATH; unset PYTHONPATH; chimera --silent --script %s; export PYTHONPATH=$TEMP"%script
    try:
        LOGGER.info("Running Chimera script: %s"%cmd)
        # status  = os.system(cmd)
        p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
        (stdout, stderr) = p.communicate()

        if stdout and len(str(stdout)) > 0:
          LOGGER.info("chimera stdout: %s"%str(stdout))
        if stderr and len(str(stderr)) > 0:
          LOGGER.info("chimera stderr: %s"%str(stderr))
        if p.returncode != 0:
          raise Exception("Chimera process returned non-zero exit status.")
    except Exception as e:
        LOGGER.exception("Chimera failed")
        raise

    #
    # Write the original structure to the output directory
    #
    LOGGER.info("Writing original structure to file %s"%originalCIFfilename)
    cifio.set_structure(structure)
    cifio.save(originalCIFfilename)

    ## Write pathogenic, neutral, candidate, and/or quantiative trait Chimera attribute files

    if args.variants and args.quantitative:
        sorted_complex_df = complex_df[complex_df["qt"].notnull()].sort_values(by=["chain","unp_pos"])
    elif sufficient_neutral_variants or sufficient_pathogenic_variants or args.variants:
        sorted_complex_df = complex_df[complex_df["dcode"].notnull()].sort_values(by=["chain","unp_pos"])
    else:
        sorted_complex_df = pd.DataFrame([])


    if not sorted_complex_df.empty:
        # Dump the neutral variants which have been loaded
        write_input_variants_to_chimera_attributes(sorted_complex_df.loc[sorted_complex_df["dcode"]==0],"neutral")
        # Dump the pathological variants which have been loaded
        write_input_variants_to_chimera_attributes(sorted_complex_df.loc[sorted_complex_df["dcode"]==1],"pathogenic")

        write_input_variants_to_chimera_attributes(sorted_complex_df.loc[sorted_complex_df["qt"].notnull()],"quantitative",third_column="qt")

        # WAIT - these double neutral+pathogenic character ones should no long be in there
        # Label as 0.5 if both neutral and pathogenic variants are mapped to a residue
        filterwarnings('ignore')
        vt          = sorted_complex_df.drop_duplicates(["chain","resid"])
        vt["dcode"] = sorted_complex_df.groupby(["chain","resid"]).apply(lambda x: np.mean(x["dcode"])).values
        resetwarnings()
    else:
        vt = pd.DataFrame([])

    # Write attribute file with dcode for all the variants, with the "attribute:" set to pathogenicty
    write_input_variants_to_chimera_attributes(vt, "variants", attribute_text="pathogenicity", third_column="dcode")

    # If requested, run the univariate K and bivariate D analyses
    if sufficient_neutral_variants and sufficient_pathogenic_variants and args.ripley:
        LOGGER.info("Neutral and Pathogenic residue counts are sufficient to calculate both Ripley's univariate K and bivariate D...")
        pK,pKz,pKt,nK,nKz,nKt,D,Dz,Dt = ripley(complex_df,args.permutations)
        qK = qKz = qKt = qKp = np.nan
    # If only one dataset present, process the univariate K only
    elif args.ripley and (sufficient_neutral_variants or sufficient_pathogenic_variants):
        # Distance matrix for all residues
        D = squareform(pdist(complex_df[['x','y','z']]))
        D[np.identity(D.shape[0],dtype=bool)] = np.nan
        if sufficient_neutral_variants:
            n = (complex_df["dcode"]==0).astype(int)
            nT,nK,nKz,nKt,nKp = uniK_calc_and_plot(D,n,args.permutations,"neutral")
            pK = pKz = pKt = np.nan
            qK = qKz = qKt = qKp = np.nan
        if sufficient_pathogenic_variants:
            p = (complex_df["dcode"]==1).astype(int)
            pT,pK,pKz,pKt,pKp = uniK_calc_and_plot(D,p,args.permutations,"pathogenic")
            nK = nKz = nKt = np.nan
            qK = qKz = qKt = qKp = np.nan
        LOGGER.info("Insufficient variant counts to perform Ripley's D analysis.")
        LOGGER.info("Using default neighbor-weight parameters for PathProx.")
        args.radius = "NW"
        D = Dz = Dt = np.nan
    # We input some quantitative annotations and have a ripley...
    elif args.quantitative and args.ripley:
        LOGGER.info("Calculating Ripley's univariate weighted K...")
        D = squareform(pdist(complex_df[['x','y','z']]))
        D[np.identity(D.shape[0],dtype=bool)] = np.nan
        qT,qK,qKz,qKt,qKp = uniK_calc_and_plot(D,complex_df['qt'],weights=True)
        pK = pKz = pKt = np.nan
        nK = nKz = nKt = np.nan
        D  = Dz  = Dt  = np.nan
    # If none of the datasets are populated, skip the Ripley analyses
    else:
        if args.ripley: # < We only issue the warning if they wanted ripley analysis
            LOGGER.info("Insufficient variant counts to perform Ripley's K/D analyses.")
            LOGGER.info("Using default neighbor-weight parameters for PathProx.")
        # Set all Ripley's results to NaN
        pK = pKz = pKt = nK = nKz = nKt = qK = qKz = qKt = D = Dz = Dt = np.nan
        args.radius = "NW"

    LOGGER.info("Phase 3 - Perform Pathprox cross-validation and report performance");
    # Run the PathProx cross-validation and report performance
    pathogenic_coordinates  = complex_df.loc[complex_df["dcode"]==1,['x','y','z']]     # Pathogenic
    neutral_coordinates  = complex_df.loc[complex_df["dcode"]==0,['x','y','z']]     # Neutral
    if not args.no_blosum:
        LOGGER.info("Integrating blosum-100 weights for all pathogenic and neutral residue");
        pathogenic_blosum100_weights = complex_df.loc[complex_df["dcode"]==1,"blosum100"]
        neutral_blosum100_weights = complex_df.loc[complex_df["dcode"]==0,"blosum100"]

    # Determine the radius or NeighborWeight parameters
    if args.radius == "NW":
        nwlb,nwub = args.nwlb,args.nwub
        LOGGER.info("args.radius='NW'.  NeighborWeight bounds set to args.nwlb=%f,args.nwub=%f",args.nwlb,args.nwub);
    elif args.radius == "K":
        nwlb = nwub = pKt
        LOGGER.info("args.radius='K'. NeighborWeight() will return 1.0 for distances <= pKt=%f, 0 for distances > %f",nwlb,nwub);
    elif args.radius == "D":
        nwlb = nwub = Dt
        LOGGER.info("args.radius='D'. NeighborWeight() will return 1.0 for distances <= Dt=%f, 0 for distances > %f",nwlb,nwub);
    else:
        nwlb = nwub = args.radius
        LOGGER.info("args.radius=%s. NeighborWeight() will return 1.0 for distances <= %f, 0 for distances > %f",nwlb,nwub);

    # Measure the predictive performance of PathProx
    if sufficient_neutral_variants and sufficient_pathogenic_variants:
        if args.no_blosum:
            LOGGER.info("args.no_blosum=%s Measuring unweighted PathProx cross-validation performance...",args.no_blosum)
            pathogenic_cross_validation_scores = pathprox(pathogenic_coordinates,pathogenic_coordinates,neutral_coordinates,nwlb=nwlb,nwub=nwub,cross_validation_flag="P")
            neutral_cross_validation_scores = pathprox(neutral_coordinates,pathogenic_coordinates,neutral_coordinates,nwlb=nwlb,nwub=nwub,cross_validation_flag="N")
        else:
            LOGGER.info("args.no_blosum=%s Measuring Blosum100-weighted PathProx cross-validation performance...",args.no_blosum)
            pathogenic_cross_validation_scores = pathprox(pathogenic_coordinates,pathogenic_coordinates,neutral_coordinates,
                nwlb=nwlb,nwub=nwub,cross_validation_flag="P",weights=(pathogenic_blosum100_weights,pathogenic_blosum100_weights,neutral_blosum100_weights))
            neutral_cross_validation_scores = pathprox(neutral_coordinates,pathogenic_coordinates,neutral_coordinates,
                nwlb=nwlb,nwub=nwub,cross_validation_flag="N",weights=(neutral_blosum100_weights,pathogenic_blosum100_weights,neutral_blosum100_weights))
            fpr,tpr,roc_auc,prec,rec,pr_auc = calc_auc(pathogenic_cross_validation_scores,neutral_cross_validation_scores)
        LOGGER.info("PathProx ROC AUC: %.2f"%roc_auc)
        LOGGER.info("PathProx PR  AUC: %.2f"%pr_auc)
        # Plot the ROC and PR curves
        fig_roc = plot_roc(fpr,tpr,save=False)
        fig_pr  = plot_pr(rec,prec,pr_auc=pr_auc,save=False)
        res = np.c_[fpr,tpr]
        np.savetxt("%s_pathprox_roc.txt.gz"%args.label,res,"%.4g",'\t')
        res = np.c_[prec,rec]
        np.savetxt("%s_pathprox_pr.txt.gz"%args.label,res,"%.4g",'\t')
        if args.radius in ("K","D"):
            nwlb,nwub = nKt,nKt
        LOGGER.info("Measuring Neutral Constraint cross-validation performance...")
        if args.no_blosum:
            ascores = uniprox(pathogenic_coordinates,neutral_coordinates,nwlb=nwlb,nwub=nwub)
            bscores = uniprox(neutral_coordinates,neutral_coordinates,nwlb=nwlb,nwub=nwub)
        else:
            ascores = uniprox(pathogenic_coordinates,neutral_coordinates,nwlb=nwlb,nwub=nwub)# ,w=(Aw,Bw)) DO NOT WEIGHT PATHOGENIC VARIANTS
            bscores = uniprox(neutral_coordinates,neutral_coordinates,nwlb=nwlb,nwub=nwub,weights=(neutral_blosum100_weights,neutral_blosum100_weights))
        fpr,tpr,roc_auc,prec,rec,pr_auc = calc_auc(list(1-np.array(ascores)),list(1-np.array(bscores)))
        LOGGER.info("Neutral Constraint ROC AUC: %.2f"%roc_auc)
        LOGGER.info("Neutral Constraint PR  AUC: %.2f"%pr_auc)
        # Plot the ROC and PR curves
        fig_roc = plot_roc(fpr,tpr,ax=fig_roc,save=False,label="Neutral Constraint",color="blue")
        fig_pr  = plot_pr(rec,prec,ax=fig_pr,pr_auc=pr_auc,save=False,label="Neutral Constraint",color="blue")
        if args.radius in ("K","D"):
            nwlb,nwub = pKt,pKt
        LOGGER.info("Measuring Pathogenic Constraint cross-validation performance...")
        if args.no_blosum:
            ascores = uniprox(pathogenic_coordinates,pathogenic_coordinates,nwlb=nwlb,nwub=nwub)
            bscores = uniprox(neutral_coordinates,pathogenic_coordinates,nwlb=nwlb,nwub=nwub)
        else:
            ascores = uniprox(pathogenic_coordinates,pathogenic_coordinates,nwlb=nwlb,nwub=nwub)# ,w=(Aw,Aw)) DO NOT WEIGHT PATHOGENIC VARIANTS
            bscores = uniprox(neutral_coordinates,pathogenic_coordinates,nwlb=nwlb,nwub=nwub,weights=(neutral_blosum100_weights,pathogenic_blosum100_weights))
        fpr,tpr,roc_auc,prec,rec,pr_auc = calc_auc(ascores,bscores)
        LOGGER.info("Pathogenic Constraint ROC AUC: %.2f"%roc_auc)
        LOGGER.info("Pathogenic Constraint PR  AUC: %.2f"%pr_auc)
        # Plot and save the ROC and PR curves
        fig_roc = plot_roc(fpr,tpr,ax=fig_roc,label="Pathogenic Constraint",color="red")
        fig_pr  = plot_pr(rec,prec,ax=fig_pr,pr_auc=pr_auc,label="Pathogenic Constraint",color="red")
    else:
        roc_auc = pr_auc = np.nan

    # Calculate PathProx scores for all residues
    all_coordinates  = complex_df[['x','y','z']]
    if sufficient_neutral_variants and sufficient_pathogenic_variants:
        if args.no_blosum:
            LOGGER.info("Calculating unweighted PathProx scores (and z-scores)...")
            complex_df["pathprox"] = pathprox(C,A,B,nwlb=nwlb,nwub=nwub)
        else:
            LOGGER.info("Calculating Blosum-100-weighted PathProx scores (and z-scores)...")
            all_blosum100_weights= complex_df["blosum100"]
            complex_df["pathprox"] = pathprox(all_coordinates,pathogenic_coordinates,neutral_coordinates,nwlb=nwlb,nwub=nwub,weights=(neutral_blosum100_weights,pathogenic_blosum100_weights,neutral_blosum100_weights))
    # Calculate neutral constraint scores for all residues
    if sufficient_neutral_variants:
        LOGGER.info("Calculating neutral constraint scores (and z-scores)...")
        # Now calculate for all residues
        if args.no_blosum:
            complex_df["neutcon"] = 1-np.array(uniprox(all_coordinates,neutral_coordinates,nwlb=nwlb,nwub=nwub))
        else:
            all_blosum100_weights= complex_df["blosum100"]
            complex_df["neutcon"] = 1-np.array(uniprox(all_coordinates,neutral_coordinates,nwlb=nwlb,nwub=nwub,weights=(all_blosum100_weights,neutral_blosum100_weights)))
        # nneut = (complex_df["dcode"]==0).sum()
        # pneutcon = [uniprox(complex_df[['x','y','z']],
        #                     complex_df.ix[np.random.choice(complex_df.index,nneut),['x','y','z']],
        #                     nwlb=nwlb,nwub=nwub) for i in range(args.permutations)]
        # pneutcon = np.concatenate(([complex_df['neutcon']],pneutcon))
        # complex_df["neutcon_z"] = np.nan_to_num(zscore(pneutcon)[0])

    # Calculate pathogenic constraint scores for all residues
    if sufficient_pathogenic_variants:
        LOGGER.info("Calculating pathogenic constraint scores (and z-scores)...")
        # Now calculate for all residues
        if args.no_blosum:
            complex_df["pathcon"] = uniprox(all_coordinates,pathogenic_coordinates,nwlb=nwlb,nwub=nwub)
        else:
            complex_df["pathcon"] = uniprox(all_coordinates,pathogenic_coordinates,nwlb=nwlb,nwub=nwub)# ,w=(Cw,Aw)) DO NOT WEIGHT PATHOGENIC VARIANTS
        # npath = (complex_df["dcode"]==1).sum()
        # ppathcon = [uniprox(complex_df[['x','y','z']],
        #                     complex_df.ix[np.random.choice(complex_df.index,npath),['x','y','z']],
        #                     nwlb=nwlb,nwub=nwub) for i in range(args.permutations)]
        # ppathcon = np.concatenate(([complex_df['pathcon']],ppathcon))
        # complex_df["pathcon_z"] = np.nan_to_num(zscore(ppathcon)[0])
    # Calculate quantitative trait constraint scores for all residues
    if args.quantitative:
        LOGGER.info("Calculating quantitative trait constraint scores (and z-scores)...")
        if args.radius in ("K","D"):
            nwlb,nwub = qKt,qKt
        complex_df["qtprox"] = qtprox(complex_df[['x','y','z']],complex_df[['x','y','z','qt']],nwlb=nwlb,nwub=nwub)

    complex_df.to_csv(args.label + "_complex_df.csv",sep='\t',header=True,index=False)

    # Report PathProx/constraint scores for candidate variants
    if args.variants and args.quantitative:
        LOGGER.info("Quantitative constraint scores for candidate missense variants:")
        LOGGER.info(complex_df.loc[complex_df["dcode"]<0,["unp_pos","ref","alt","qtprox"]].sort_values( \
              by=["qtprox","unp_pos"],ascending=[False,True]).groupby(["unp_pos"]).apply(np.mean).to_string(index=False))
    elif args.variants:
        if sufficient_neutral_variants:
            LOGGER.info("Neutral constraint scores for candidate missense variants:")
            LOGGER.info(complex_df.loc[complex_df["dcode"]<0,["unp_pos","ref","alt","neutcon"]].sort_values( \
              by=["neutcon","unp_pos"],ascending=[False,True]).groupby(["unp_pos"]).apply(np.mean).to_string(index=False))
        if sufficient_pathogenic_variants:
            LOGGER.info("Pathogenic constraint scores for candidate missense variants:")
            LOGGER.info(complex_df.loc[complex_df["dcode"]<0,["unp_pos","ref","alt","pathcon"]].sort_values( \
              by=["pathcon","unp_pos"],ascending=[False,True]).groupby(["unp_pos"]).apply(np.mean).to_string(index=False))
        if sufficient_neutral_variants and sufficient_pathogenic_variants: # meaning, we had BOTH enough pathogenic and enough neutrals
            LOGGER.info("PathProx scores for candidate missense variants:")
            LOGGER.info(complex_df.loc[complex_df["dcode"]<0,["unp_pos","ref","alt","pathprox"]].sort_values( \
                by=["pathprox","unp_pos"],ascending=[False,True]).groupby(["unp_pos"]).apply(np.mean).to_string(index=False))
        else:
            LOGGER.warning("Pathprox scores not output due to insufficient neutral AND pathogenic variants")


    # Write a particular score column to Chimera attribute file with obvious file name
    def write_score_column_to_chimera_attributes(complex_df: pd.DataFrame,score_column: str):
        # score_column: column headers of interest for per-residue scores/values
        #                ....  neutcon, pathcon, pathprox, or qtprox
        with open("%s_%s.attr"%(args.label,score_column),'w') as attribute_file, \
             open("%s_renum_%s.attr"%(args.label,score_column),'w') as attribute_renum_file:
            attribute_file.write("attribute: %s\n"%score_column)
            attribute_file.write("match mode: 1-to-1\n")
            attribute_file.write("recipient: residues\n")
            attribute_renum_file.write("attribute: %s\n"%score_column)
            attribute_renum_file.write("match mode: 1-to-1\n")
            attribute_renum_file.write("recipient: residues\n")

            for i,row in complex_df.iterrows():
                attribute_file.write("\t:%s.%s\t%.3f\n"%(format_resid(row["resid"]),row["chain"],row[score_column]))
                attribute_renum_file.write("\t:%d.%s\t%.3f\n"%(row["unp_pos"],row["chain"],row[score_column]))

    # Neutral constraint
    if sufficient_neutral_variants:
        write_score_column_to_chimera_attributes(complex_df,"neutcon")

    # Pathogenic constraint
    if sufficient_pathogenic_variants:
        write_score_column_to_chimera_attributes(complex_df,"pathcon")

    # PathProx
    if sufficient_neutral_variants and sufficient_pathogenic_variants:
        write_score_column_to_chimera_attributes(complex_df,"pathprox")

    # Quantitative constraint
    if args.quantitative:
        write_score_column_to_chimera_attributes(complex_df,"qtprox")

    # Write summary results to file
    # Ripley's K/D results
    head = ["Kz_path","Kt_path","Kz_neut","Kt_neut","Dz","Dt","Kz_quant","Kt_quant"]
    vals = [pKz,pKt,nKz,nKt,Dz,Dt,qKz,qKt]
    # PathProx Cross Validation, ROC and PR
    head.extend(["roc_auc","pr_auc"])
    vals.extend([roc_auc,pr_auc])
    # Extract and report on Variants of Unknown Significance, and their analysis
    if (complex_df["dcode"]<0).sum()>0:
        vus_df = complex_df.loc[complex_df["dcode"]<0,["ref","unp_pos","alt"]]
        vus_strs = vus_df["ref"].str.cat([vus_df["unp_pos"].astype(str),vus_df["alt"]])
        # PathProx scores
        if sufficient_pathogenic_variants and sufficient_neutral_variants:
            pp_vus = complex_df.loc[complex_df["dcode"]<0,["unp_pos","ref","alt","pathprox"]].sort_values( \
                  by=["pathprox","unp_pos"],ascending=[False,True]).groupby(["unp_pos"]).apply(np.mean)
            head.extend(["%s_pathprox"%m for m in vus_strs])
            vals.extend(list(pp_vus["pathprox"].values))
        # Pathogenic constraint scores
        if sufficient_pathogenic_variants:
            pc_vus = complex_df.loc[complex_df["dcode"]<0,["unp_pos","ref","alt","pathcon"]].sort_values( \
                  by=["pathcon","unp_pos"],ascending=[False,True]).groupby(["unp_pos"]).apply(np.mean)
            head.extend(["%s_pathcon"%m for m in vus_strs])
            vals.extend(list(pc_vus["pathcon"].values))
        # Neutral constraint scores
        if sufficient_neutral_variants:
            nc_vus = complex_df.loc[complex_df["dcode"]<0,["unp_pos","ref","alt","neutcon"]].sort_values( \
                  by=["neutcon","unp_pos"],ascending=[False,True]).groupby(["unp_pos"]).apply(np.mean)
            head.extend(["%s_neutcon"%m for m in vus_strs])
            vals.extend(list(nc_vus["neutcon"].values))

    # Logically if any of these "nearest pathogenic residue" calculations are set, then they all should be
    # So output these in a way that will crash hard if they are not all set
    if (pClosestPathogenicMutationPDB or pNextPathogenicMutationPDB or 
          pClosestPathogenicMutationUNP or pNextPathogenicMutationUNP or
          pClosestPathogenicDistance or pNextPathogenicSeqDelta):
        head.extend(["ClosestPDB","ClosestUNP","ClosestDistance"])
        vals.extend([pClosestPathogenicMutationPDB,pClosestPathogenicMutationUNP,pClosestPathogenicDistance])
        head.extend(["NextPDB","NextUNP","SeqDelta"])
        vals.extend([pNextPathogenicMutationPDB,pNextPathogenicMutationUNP,pNextPathogenicSeqDelta])

    # Write all results to file
    summary_filename = "%s_summary.csv"%args.label
    summary_strings = []
    for summary_tuple in zip(head,vals):
        summary_strings.append("  %-15.15s %s"%(summary_tuple[0],summary_tuple[1]))
    LOGGER.info("Writing summary of pathprox results to %s:\n%s",summary_filename,'\n'.join(summary_strings))
    with open(summary_filename,'w') as fout:
        writer = csv.writer(fout,delimiter='\t')
        writer.writerow(head)
        writer.writerow(vals)

    # Generating structural images
    LOGGER.info("Visualizing with Chimera (this may take a while)...")
    params = [args.label,structure.id,str(is_biounit)]
    # Run the Chimera visualization script
    script = '"%s/pathvis.py %s"'%(script_dir,' '.join(params))
    cmd    = "TEMP=$PYTHONPATH; unset PYTHONPATH; %s --nogui --silent --script %s; export PYTHONPATH=$TEMP"%(config_dict['chimera_headless'],script)
    # Allow Mac OSX to use the GUI window
    if platform.system() == 'Darwin':
        cmd  = "TEMP=$PYTHONPATH; unset PYTHONPATH; chimera --silent --script %s; export PYTHONPATH=$TEMP"%script
    try:
        LOGGER.info("Running Chimera script: %s"%cmd)
        # status  = os.system(cmd)
        p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
        (stdout, stderr) = p.communicate()

        if stdout and len(str(stdout)) > 0:
            LOGGER.info("chimera stdout: %s"%str(stdout))
        if stderr and len(str(stderr)) > 0:
            LOGGER.info("chimera stderr: %s"%str(stderr))
        if p.returncode != 0:
            raise Exception("Chimera process returned non-zero exit status.")
    except Exception as e:
        LOGGER.exception("Chimera failed")
        raise
    # Generate a dictionary of information that psb_rep.py can use to create a robust ngl viewer 
    # experience for users.  Javascript format is human readable, and opens possibility for integration with javascript
    json_filename = "%s_ResiduesOfInterest.json"%args.label
    
    residuesOfInterest = {}  
     
    # Strip out all punctuation to make a reasonablly unique token if psb_rep.py needs to generate global variables in .html
    residuesOfInterest['unique_div_name'] = os.path.basename(args.uniquekey).translate(str.maketrans('','',string.punctuation))

    # pdb_rep.py will point ngl to the psb with HELIX and SHEET information (Secondary Structure)
    residuesOfInterest['pdbSSfilename'] =  os.path.join(args.outdir,renumberedPDBfilenameSS)
    residuesOfInterest['cifSSfilename'] =  os.path.join(args.outdir,renumberedCIFfilenameSS)

    # write out the various residues in per-chain dictionaries.  psb_rep.py will use these later 
    from collections import defaultdict 
    variant_residues = defaultdict(list)
    neutral_residues = defaultdict(list)
    pathogenic_residues = defaultdict(list)
    attribute_residues =defaultdict(list) 
    for i,r in complex_df.iterrows():
      residue = r["unp_pos"]
      chain =  r["chain"]
      if r['dcode'] == 0:
        neutral_residues[chain].append(residue)
      elif r['dcode'] == 1:
        pathogenic_residues[chain].append(residue)
      elif r['dcode'] == -1:
        variant_residues[chain].append(residue)
      # If chimera-specific per-residue attributes are in memory, dump them too
      if (r['qt'] is not None) and (not np.isnan(r['qt'])):
        attribute_residues[chain].append({'residue': residue, 'attribute': r['qt']})

    residuesOfInterest['variants'] = variant_residues
    residuesOfInterest['neutrals'] = neutral_residues
    residuesOfInterest['pathogenics'] = pathogenic_residues
    residuesOfInterest['user_attributes'] = attribute_residues

    with open(json_filename,'w') as f:
         json.dump(residuesOfInterest, f)
    LOGGER.info("%d variants, %d neutral, %d pathogenic, and %d user_attributes recorded to %s"%(
      len(residuesOfInterest['variants']),
      len(residuesOfInterest['neutrals']),
      len(residuesOfInterest['pathogenics']),
      len(residuesOfInterest['user_attributes']),
      json_filename) )
    
    # Return back to the right directory!
    os.chdir(curdir)
    statusdir_info('Completed')
    # Mark this analysis as complete
    open(os.path.join(statusdir,"complete"),'w').close()

    LOGGER.info("All analyses complete.")
    sys.exit(0)



    ### IMPORTANT - The legacy code below can process MULTIPLE PDB files.  And we need to think about this!

    """unp_flag,olabel = False,args.label
    args.quantitative  = args.quantitative
    curdir = os.getcwd() # Note the current directory
    for sid,bio,cf in get_coord_files(args.entity,io):

      # Begin each analysis in the original directory
      if os.getcwd() != curdir:
        os.chdir(curdir)

      # Reinitialize variant sets for each structure
      q = parse_qt(args.quantitative)
      p = parse_variants("Pathogenic",args.pathogenic)
      n = parse_variants("Neutral",args.neutral)
      c = parse_variants("Candidate",args.variants)

      if unp_flag:
        tlabel = olabel.split('_')
        tlabel.insert(1,sid)
        args.label = '_'.join(tlabel)
        LOGGER.info("Sub-analysis label set to %s"%args.label)

      # Check for existing results
      if os.path.exists("%s/.%s.complete"%(args.outdir,args.label)):
        # Do not overwrite unless specified by the user
        if not args.overwrite:
          msg = "\nStructure %s[%s] has been processed. Use --overwrite to overwrite.\n"%(sid,bio)
          LOGGER.warning(msg); continue
        else:
          # Remove the complete flag and reanalyze
          os.remove("%s/.%s.complete"%(args.outdir,args.label))

      LOGGER.info("Processing %s[%s]..."%(sid,bio))

      # Read and renumber the coordinate file:
      LOGGER.info("Reading coordinates from %s..."%cf)
      s_renum,_,_,_    = read_coord_file(cf,sid,bio,chain=args.chain,
                                      fasta=args.fasta,residues=args.use_residues,
                                      renumber=True,AAseq=AAseq)

      unp = None
      # Read the coordinate file, align, etc. Do not renumber
      s,indb,unp,chains = read_coord_file(cf,sid,bio,chain=args.chain,
                                      fasta=args.fasta,residues=args.use_residues,
                                      renumber=False,AAseq=AAseq)
      if args.fasta:
        LOGGER.info("UniProt AC derived from FASTA: %s"%unp)
      else:
        LOGGER.info("UniProt AC: %s"%unp)

      # Check that any user-specified chain is present in the structure
      if args.chain and args.chain not in chains:
        msg = "Biological assembly %s does not contain chain %s. Skipping\n"%(bio,args.chain)
        LOGGER.warning(msg); continue

      # Reduce the considered chains to the user-specified chain if present
      chains = chains if not args.chain else [args.chain]

      # If user-supplied variants did not include chain specifiers, add the
      # user-specified or inferred chain IDs to each variant as necessary.
      p = [v if len(v)==4 else v+[ch] for v in p for ch in chains]
      n = [v if len(v)==4 else v+[ch] for v in n for ch in chains]
      c = [v if len(v)==4 else v+[ch] for v in c for ch in chains]
      if args.quantitative: #FIXME: Necessary condition?
        q = [v+[qt] if len(v)==4 else v+[ch]+[qt] for v,qt in q for ch in chains]

      LOGGER.info("User supplied variant counts:")
      LOGGER.info(" Neutral:      %d"%len(n))
      LOGGER.info(" Pathogenic:   %d"%len(p))
      LOGGER.info(" Candidates:   %d"%len(c))
      LOGGER.info(" Quantitative: %d"%len(q))

      LOGGER.info("Querying additional variants for %s (%s)..."%(sid,unp))

      # Supplement with any requested variant datasets
      if args.add_1kg:
        n.extend(query_1kg(io,sid,unp,chains,indb))
      if args.add_exac:
        LOGGER.info("Adding ExAC missense variants...")
        n.extend(query_exac(io,sid,unp,chains,indb))
      if args.add_gnomad:
        LOGGER.info("Adding gnomAD missense variants...")
        n.extend(query_gnomad(io,sid,unp,chains,indb))
      if args.add_benign:
        LOGGER.info("Adding ClinVar benign variants...")
        n.extend(query_benign(io,sid,unp,chains,indb))
      if args.add_pathogenic:
        LOGGER.info("Adding ClinVar pathogenic variants...")
        p.extend(query_pathogenic(io,sid,unp,chains,indb))
      if args.add_drug:
        LOGGER.info("Adding ClinVar drug-related variants...")
        p.extend(query_drug(io,sid,unp,chains,indb))
      if args.add_cosmic:
        LOGGER.info("Adding COSMIC recurrent variants...")
        p.extend(query_cosmic(io,sid,unp,chains,indb))
      if args.add_tcga:
        p.extend(query_tcga(io,sid,unp,chains,indb))

      LOGGER.info("Final supplied variant counts:")
      LOGGER.info(" Neutral:      %d"%len(n))
      LOGGER.info(" Pathogenic:   %d"%len(p))
      LOGGER.info(" Candidates:   %d"%len(c))
      LOGGER.info(" Quantitative: %d"%len(q))

      # Create (and move to) output directory
      args.outdir = move_to_outdir(args.outdir)

      # Annotate coordinate file with pathogenic, neutral, and candidate labels
      complex_df   = var2coord(s,p,n,c,q)

      # On exit from above
      # complex_df["dcode"] is 0 for neutral variants, 1 for pathogenic, -1 for our target variant(s) of interest
      # and NaN for the structural entries with no assignments
    """
 
      



    """ Chris has note into Mike Sivley as this seems to do nothing for us
      # Write scores to tab-delimited file
      # Neutral constraint
      if sufficient_neutral_variants:
        c = complex_df.sort_values(by="neutcon").drop_duplicates(["chain","resid"])
        c = complex_df.sort_values(by=["chain","resid"])
        c.to_csv("%s_neutcon.txt"%args.label,sep='\t',header=True,index=False)
      # Pathogenic constraint
      if sufficient_pathogenic_variants:
        c = complex_df.sort_values(by="pathcon").drop_duplicates(["chain","resid"])
        c = complex_df.sort_values(by=["chain","resid"])
        c.to_csv("%s_pathcon.txt"%args.label,sep='\t',header=True,index=False)
      # PathProx
      if sufficient_neutral_variants and sufficient_pathogenic_variants:
        c = complex_df.sort_values(by="pathprox").drop_duplicates(["chain","resid"])
        c = complex_df.sort_values(by=["chain","resid"])
        c.to_csv("%s_pathprox.txt"%args.label,sep='\t',header=True,index=False)
      # Quantative trait constraint
      if args.quantitative:
        c = complex_df.sort_values(by="qtprox").drop_duplicates(["chain","resid"])
        c = complex_df.sort_values(by=["chain","resid"])
        c.to_csv("%s_qtprox_%s.txt"%(args.label,args.qlabel),sep='\t',header=True,index=False)
    """
