#!/usr/bin/env python
#
# Project        : PathProx
# Filename       : pathprox.py  (pathprox2.py reflects new logging, configuration, output structure of pipelineV2 architecture)
# Authors        : R. Michael Sivley
#                : Xiaoyi Dou
# Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2016-11-04
# Description    : Predicts the pathogenicity of missense variants by their
#                : spatial position relative to known pathogenic and neutral
#                : variation in protein structure.
#=============================================================================#
"""Pathprox predicts (scores) the pathogenicity of a missense variant by its
spatial position relative to known pathogenic and neutral variation in its 
3D protein structure."""

#logging
import logging
import json
import re

from logging.handlers import RotatingFileHandler
from logging import handlers
PATHPROX_LOG_FORMATTER = logging.Formatter('%(asctime)s %(levelname)-4s [%(filename)16s:%(lineno)4d] %(message)s',datefmt='%H:%M:%S')
ch = logging.StreamHandler()
ch.setFormatter(PATHPROX_LOG_FORMATTER)
LOGGER = logging.getLogger()
LOGGER.addHandler(ch)
LOGGER.setLevel(logging.INFO)

import inspect # For status updates

## Package Dependenecies ##
# Standardction with a given name return the same logger instance. This means that logger instances never need to be passed between d
import sys,os,shutil,gzip,csv,platform,grp,stat,subprocess
import string
import subprocess as sp
from time import strftime
from itertools import combinations
from collections import OrderedDict
import pprint
from subprocess import Popen, PIPE

cache_dir = '/tmp/sqlcache'  # Must be overridden early

# Warnings
from warnings import filterwarnings,resetwarnings

# Numerical
import pandas as pd
import numpy as np
import math
TOL = 1e-5 # zero-tolerancection with a given name return the same logger instance. This means that logger instances never need to be passed between d

# Stats
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.stats import norm, percentileofscore, mannwhitneyu
from scipy.spatial import KDTree

# Removed October - because built-in zscore chokes on sigma=0 from scipy.stats.mstats import zscore

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

# Biopython
filterwarnings('ignore',category=ImportWarning)
from Bio.SubsMat.MatrixInfo import blosum100
minv,maxv = np.min(list(blosum100.values())),np.max(list(blosum100.values()))
for key,val in list(blosum100.items()):
  # Transform to 0..1, then invert so that severe mutations get high scores
  val = float(val)
  blosum100[key]             = 1-((val-minv)/(maxv-minv))
  # Store in reverse for easy queries
  blosum100[(key[1],key[0])] = 1-((val-minv)/(maxv-minv))

from Bio.PDB.PDBParser import PDBParser
resetwarnings()
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# PDBMap
from pdbmap import PDBMapIO,PDBMapParser,PDBMapStructure,PDBMapSwiss
from pdbmap import PDBMapProtein,PDBMapTranscript,PDBMapAlignment
from lib.amino_acids import longer_names


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

def read_fasta(fasta):
  """ Aligns the observed sequence with a user-provided reference """
  with open(fasta,'r') as fin:
    try:
      unp,refid = fin.readline().split('|')[1:3]
    except:
      msg = "Malformed fasta header. Please provide unaltered UniProt fasta files.\n"
      LOGGER.critical(msg); 
      sys_exit_failure(msg)
    if len(unp.split('-')[0])!=6 or not unp[0].isalpha():
      msg = "Fasta header |%s| lacks uniprot ID in position 2. Please provide unaltered UniProt fasta files."%unp
      LOGGER.critical(msg); 
      sys_exit_failure(msg)
    refid     = refid.split()[0]
    return unp,refid,''.join([l.strip() for l in fin.readlines() if l[0]!=">"])

def query_alignment_new(sid):
  """ Query the reference->observed alignment from PDBMap """
  s   = "SELECT unp,trans_seqid,b.chain,actr.chain_res_num as chain_seqid,b.rescode as pdb_res,d.rescode as trans_res "
  s  += "FROM AlignChainTranscript act "
  s  += "INNER JOIN Residue b ON act.label=b.label AND act.structid=b.structid "
  s  += "LEFT OUTER JOIN AlignChainTranscriptResidue actr on actr.al_id = act.al_id "
  s  += "AND act.chain=b.chain AND actr.chain_res_num=b.seqid AND actr.chain_res_icode = b.icode "
  s  += "INNER JOIN Chain c ON b.label=c.label AND b.structid=c.structid "
  s  += "AND b.label=c.label AND b.biounit=c.biounit "
  s  += "AND b.model=c.model AND b.chain=c.chain "
  s  += "INNER JOIN Transcript d ON act.label=d.label AND act.transcript=d.transcript "
  s  += "AND actr.trans_seqid=d.seqid "
  w   = "WHERE (act.label='pdb' or act.label='modbase' or act.label='swiss') AND act.structid=%s AND b.biounit=0 AND b.model=0"
  q   = s+w
  res = list(io.secure_cached_query(cache_dir,q,(sid,)))
  if (not res) or (len(res) == 0): # Then we are having a very bad day indeed
    return None,None

  unp = res[0]["unp"]
  # chain,chain_seqid -> ref_seqid,rescode
  chains = set([r["chain"] for r in res])
  aln = dict(((r["chain"],r["chain_seqid"]),
              (r["trans_seqid"],r["pdb_res"])) for r in res)
  return unp,aln

def query_alignment(sid): #Deprecate this one after pdbs in new structures
  """ Query the reference->observed alignment from PDBMap """
  s   = "SELECT unp,trans_seqid,b.chain,chain_seqid,b.rescode as pdb_res,d.rescode as trans_res from Alignment a "
  s  += "INNER JOIN Residue b ON a.label=b.label AND a.structid=b.structid "
  s  += "AND a.chain=b.chain AND a.chain_seqid=b.seqid "
  s  += "INNER JOIN Chain c ON b.label=c.label AND b.structid=c.structid "
  s  += "AND b.label=c.label AND b.biounit=c.biounit "
  s  += "AND b.model=c.model AND b.chain=c.chain "
  s  += "INNER JOIN Transcript d ON a.label=d.label AND a.transcript=d.transcript "
  s  += "AND a.trans_seqid=d.seqid "
  w   = "WHERE (a.label='pdb' or a.label='modbase' or a.label='swiss') AND b.structid=%s AND b.biounit=0 AND b.model=0"
  q   = s+w
  res = list(io.secure_cached_query(cache_dir,q,(sid,)))

  if (not res) or (len(res) == 0): # Then we need to query from the new AlignmentChain Setup
    return query_alignment_new(sid)
  unp = res[0]["unp"]
  # chain,chain_seqid -> ref_seqid,rescode
  chains = set([r["chain"] for r in res])
  aln = dict(((r["chain"],r["chain_seqid"]),
              (r["trans_seqid"],r["pdb_res"])) for r in res)
  return unp,aln

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

def model_lookup(io,mid):
  """ Returns coordinate files for a ModBase ID """
  PDBMapModel.load_modbase(config_dict['modbase2016_dir'],config_dict['modbase2016_summary'])
  PDBMapModel.load_modbase(config_dict['modbase2013_dir'],config_dict['modbase2013_summary'])
  f = PDBMapModel.get_coord_file(mid.upper())
  LOGGER.info("File location for %s: %s"%(mid.upper(),f))
  #f = "%s/Homo_sapiens_2016/model/%s.pdb.gz"%(args.modbase_dir,mid.upper())
  if not os.path.exists(f):
    try:
      cmd = ["xz","-d",'.'.join(f.split('.')[:-1])+'.xz']
      sp.check_call(cmd)
    except:
      msg  = "Coordinate file missing for %s\n"%mid
      msg += "Expected: %s\n"%f
      raise Exception(msg)
    cmd = ["gzip",'.'.join(f.split('.')[:-1])]
    sp.check_call(cmd)
  return [(mid,0,f)]

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

def uniprot_lookup(io,ac):
  """ Returns coordinate files for a UniProt AC """
  entities = io.load_unp(ac)
  if not entities:
    msg = "No structures or models associated with %s\n"%ac
    LOGGER.warning(msg)
    return None
  flist = []
  for etype,ac in entities:
    if etype == "structure":
      flist.extend(structure_lookup(io,ac))
    elif etype == "model":
      flist.extend(model_lookup(io,ac))
    elif etype == "swiss":
      flist.extend(swiss_lookup(io,ac))
  return flist

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

def query_cosmic(io,sid,refid=None,chains=None,indb=False):
  """ Query somatic variants (COSMIC) from PDBMap """
  if indb:
    s,w = default_var_query()
    if chains:
      w += "AND chain in (%s) "%','.join(["'%s'"%c for c in chains])
    f   = ("cosmic",sid)
    c   = ["unp_pos","ref","alt","chain"]
  else:
    s,w = sequence_var_query()
    f   = ["cosmic"]
    c   = ["unp_pos","ref","alt"]
  m   = "INNER JOIN cosmic d "
  m  += "ON b.chr=d.chr and b.start=d.start "
  m  += "AND d.cnt>1 "
  q   = s+m+w
  res = [list(r) for r in io.secure_cached_query(cache_dir,q,f,cursorclass="Cursor")]
  # If user-specified model...
  if not indb:
    # Add chains IDs
    res = [r+[c] for r in res for c in chains]
    # # Add null pdb_pos (not yet aligned)
    # res = [[None]+r for r in res]
  return res

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

def get_coord_files(entity,io):
  """ Returns the relevant coordinate file or None if no file found """
  # Check if the entity is a filename. If so, assume PDB/ENT format
  global unp_flag
  if os.path.isfile(entity):
    exts = 1 + int(os.path.basename(entity).split('.')[-1]=='gz')
    sid  = '.'.join(os.path.basename(entity).split('.')[:-exts]) 
    return [(sid,-1,entity)]
  else:
    LOGGER.info("Attempting to identify entity: %s..."%entity)
    etype = io.detect_entity_type(entity)
    if   etype == "structure":
      LOGGER.info("Input is a PDB ID. Locating coordinate file...")
      return structure_lookup(io,entity,args.chain is None)
    elif etype == "model":
      LOGGER.info("Input is a ModBase model ID. Locating coordinate file...")
      return model_lookup(io,entity)
    elif etype == "swiss":
      LOGGER.info("Input is a Swissmodel ID. Locating coordinate file...")
      return swiss_lookup(io,entity)
    elif etype == "unp":
      unp_flag = True
      LOGGER.info("Input is a UniProt AC. Identifying all relevant PDB structures and ModBase models...")
      return uniprot_lookup(io,entity)
    elif etype: # HGNC returned a UniProt ID
      unp_flag = True
      LOGGER.info("Input is an HGNC gene name. Identifying all relevant PDB structures and ModBase models...")
      LOGGER.info("HGNC: %s => UniProt: %s..."%(entity,etype))
      return uniprot_lookup(io,etype)
    else:
      msg = "Could not identify %s."%entity
      LOGGER.critical(msg);
      sys_exit_failure(msg)

def read_coord_file(coord_filename,sid,bio,chain,fasta=None,residues=None,renumber=False,AAseq=None):
  """ Reads the coordinate file into a PDBMapStructure object """
  # If no fasta provided, query Uniprot AC and alignment
  indb = False

  unp = None
  seq = AAseq
  if fasta:
    aln = {}
    unp,refid,seq = read_fasta(fasta)
  else:
    if bio < 0 and not seq:
        msg = "FASTA files or ENST... transcript id required for user-defined protein models\n"
        LOGGER.critical(msg)
        sys_exit_failure(msg)
    # if not seq:
    unp,aln = query_alignment(sid)
    indb = (unp != None)
    refid = None

  with gzip.open(coord_filename,'rt') if coord_filename.split('.')[-1]=="gz" else open(coord_filename,'r') as fin:
    filterwarnings('ignore',category=PDBConstructionWarning)
    p = PDBParser()
    s = p.get_structure(sid,fin)
    resetwarnings()
  p = PDBMapParser()
  s = p.process_structure(s,force=True)
  # Preprocess structural properties
  for m in s:
    for c in list(m.get_chains()):
      if c.id in ('',' '):
        # In Biopython 1.73, this below calls @setter.id, and sorts out the parent relationships beautifully
        c.id = 'A'

  # Reduce to the specified chain if given
  if chain:
    # Drop chains that do not match the specified chain
    for m in s:
      for c in list(m.get_chains()):
        if c.id != chain:
          LOGGER.info("Ignoring chain %s (user specified %s)."%(c.id,chain))
          c.get_parent().detach_child(c.id)

  s = PDBMapStructure(s,refseq=seq,pdb2pose={})
  # Dirty hack: Reset pdb2pose
  s._pdb2pose = dict((key,key) for key,val in s._pdb2pose.items())

  # If we in fact got an alignment from the database, then pathprox should use that
  if aln and len(aln) > 2:
    # Manually create a PDBMapAlignment from the queried alignment
    LOGGER.info("Chains in %s: %s"%
       (s.id,','.join(sorted([c.id for c in s.get_chains()]))))
    LOGGER.info("Chains in %s with pre-calculated alignments: %s"%
       (s.id,','.join(sorted([c for c in set([k[0] for k in list(aln.keys())])]))))
    for c in s.get_chains():
      LOGGER.info("Using pre-calculated alignment for chain %s"%c.id)
      refdict = dict((val[0],(val[1],"NA",0,0,0)) for key,val in aln.items() if key[0]==c.id)
      c.transcript = PDBMapTranscript("ref","ref","ref",refdict)
      c.alignment  = PDBMapAlignment(c,c.transcript)
      s.transcripts.append(c.transcript)
      s.alignments.append(c.alignment)
    # Remove any residues that do not map to any reference residues
    for m in s:
      for c in m:
        for r in list(c.get_residues()):
          if r.id not in c.alignment.pdb2seq:
            c.detach_child(r.id)

  if renumber:
    from Bio.PDB.Chain import Chain
    # Update all residue numbers to match the reference sequence
    for m in s:
      for c in list(m.get_chains()):
        # Align pdb residue numbers to the reference sequence
        if not c.alignment.trivial_alignment: # Then we need to renumber things a bit
          residues_list = list(c.get_residues())
          renumbered_chain = Chain(c.id)

          for r in residues_list:
            LOGGER.info("%s",str(r))
            LOGGER.info("%s",str(c.alignment.pdb2seq[r.id]))
            new_id = (' ',c.alignment.pdb2seq[r.id],' ')
            r.id = new_id 
            renumbered_chain.add(r)

          m.detach_child(c.id)

          # Chris Moth does NOT like the monkey patching
          # But next bit of code needs the c.alignment.perc_identity
          renumbered_chain.alignment = c.alignment

          m.add(renumbered_chain)

        
          

  # If user-provided structure, QA the chain alignment
  if bio < 0:
    # Drop chains that do not match the reference sequence
    for m in s:
      for c in list(m.get_chains()):
        if c.alignment.perc_identity<0.9:
          LOGGER.info("Chain %s does not match reference sequence. Ignoring."%c.id)
          c.get_parent().detach_child(c.id)

  # Reduce to specified residue range if given
  if residues:
    for c in list(s.get_chains()):
      for r in list(c.get_residues()):
        if (args.use_residues[0] and r.id[1] < args.use_residues[0]) or \
            (args.use_residues[1] and r.id[1] > args.use_residues[1]):
          c.detach_child(r.id)
  return s,indb,unp,[c.id for c in s.get_chains()]

def chain_match(io,unp,sid,bio):
  """ Identifies which chains are associated with a UniProt AC """
  q  = "SELECT distinct chain FROM Chain "
  q += "WHERE label in ('pdb','modbase','swiss') AND structid=%s AND biounit=%s AND unp=%s"
  return [r[0] for r in io.secure_cached_query(cache_dir,q,(sid,bio,unp),
                                          cursorclass="Cursor")]

def check_coverage(s,chain,pos,refpos=True):
  """ Test if a protein position is covered by the protein structure """
  return bool(s.get_residue(chain,pos,refpos=refpos))

# The parser has to cope with AnnnB, AaaNNNBbb, and p.AaaNNNBbb (official) input formats
# See http://varnomen.hgvs.org/
# Strings ending with .X are consired to terminate with chain indicators
from Bio.SeqUtils import seq1
def parse_variants(variantsFlavor,varset):
  if not varset:
    return []
  if os.path.isfile(varset):
    LOGGER.info("Parsing %s variants from file: %s",variantsFlavor,varset)
    raw_variants = [line.strip() for line in open(varset,'r')]
  else:
    LOGGER.info("Parsing %s variants from command line argument:\n %s",variantsFlavor,varset)
    raw_variants = [variant.strip() for variant in varset.split(',')]

  # Strip out opening "p." if it is there
  variants_without_p = []
  for v in raw_variants:
    if v.startswith("p."):
      variants_without_p.append(v[2:])
    else :
      variants_without_p.append(v)

  if len(variants_without_p) == len(raw_variants):
    LOGGER.info("Raw variant parse: %d %s variants loaded",len(variants_without_p),variantsFlavor)
  else:
    LOGGER.warning("Blank/invalid variants dected in raw variant parse.  Started with %d %s, then reduced to %d",len(raw_variants),variantsFlavor,len(variants_without_p))

  parsed_variants = []
  chains = []

  one_letter_matcher = re.compile(r"^(?P<refAA>[A-Z])(?P<position>[0-9]+)(?P<altAA>[A-Z])$")
  assert one_letter_matcher

  three_letter_matcher = re.compile(r"^(?P<refAA>[A-Z]{3})(?P<position>[0-9]+)(?P<altAA>[A-Z]{3})$")
  assert three_letter_matcher

  for v in variants_without_p:
    chain_delimiter = v.rfind('.')
    if chain_delimiter > -1:
      chains.append(v[chain_delimiter])
      v_without_chain = v[0:chain_delimiter].upper()
    else:
      chains.append(None)
      v_without_chain = v.upper()
    result = one_letter_matcher.match(v_without_chain)
    if result:
      v_dict = result.groupdict()
    else:
      result = three_letter_matcher.match(v_without_chain)
      if result:
        v_dict = result.groupdict()
        v_dict['refAA'] = seq1(v_dict['refAA'])
        v_dict['altAA'] = seq1(v_dict['altAA'])
      else:
        msg = "Unable to parse variant fragment %s"%v_without_chain
        LOGGER.critical(msg)
        sys_exit_failure(msg)
    parsed_variants.append([int(v_dict['position']),v_dict['refAA'],v_dict['altAA']])

  # pdb_pos, unp_pos, ref, alt, chain
  final_variants = [[None]+v+[chains[i]] if chains[i] else v for i,v in enumerate(parsed_variants)]
  LOGGER.info("%d variants retained by parser",len(final_variants))
  if len(final_variants) > 1:
    LOGGER.info("Variants: %s .... %s",final_variants[0],final_variants[-1])
  elif len(final_variants) == 1:
    LOGGER.info("Variant: %s",final_variants[0])
  return final_variants


def parse_qt(varset):
  if not varset:
    return []
  var,q = [],[]
  with open(varset,'r') as f:
    for line in f:
      var.append(line.split('\t')[0].strip())
      q.append(line.split('\t')[1].strip())
  var = ','.join(var)
  var = parse_variants("Quantity weights",var)
  return list(zip(var,q))

def resicom(resi):
  if resi==None: return [None,None,None]
  return np.array([np.array(a.get_coord()) for a in resi.get_unpacked_list()]).mean(axis=0)

@np.vectorize
def NeighborWeight(d,lower_bound=8.,upper_bound=24):
  """For each distance between a variant of unknown significance (VUS) and a
  known (pathogenic or neutral) variant, return a spatial proximity metric that 
  weights by emphasizing closeness to the lower_bound of the interval of interest

  See page 2 of 10 in "Methods:Protein structural" analysis in... 
  "Sivley et al. BMC Bioinformatics (2018) 19:18" 

  Parameters
  ----------
  d: Array-like set of distances between centroids of a VUS and known variants
  
  Returns
  -------
  Same shape'd np.array with weights in the interval [0.0,1.0]

  In practice, the lower_bound can equal the upper_bound, as when they
  calculated by a preceding the Ripley's K analysis"""

  if d <= TOL: # Inter-residue distances this small are indicative of "same residue to same residue" or problems
    retval = 1+1e-10
    LOGGER.warning("NeighborWeight encountered tiny inter-residue distance %g Returning %g",d,retval)
    return retval
  elif d <= lower_bound:
    return 1.
  elif d >= upper_bound:
    return 0.
  else:
    return 0.5*(np.cos(np.pi*(d-lower_bound)/(upper_bound-lower_bound))+1)

def uniprox(cands,v,nwlb=8.,nwub=24.,weights=None):
  """ Predicts a univariate constraint score from weight vector 'v' """
  if weights:
    Cw,Vw = weights
    Cw = Cw.values.reshape((Cw.size,1))
    Vw = Vw.values.reshape((Vw.size,1))
  D       = cdist(cands,v)
  D[D==0] = np.inf     # Do not use variants to score their own residues
  NWv = NeighborWeight(D,lower_bound=nwlb,upper_bound=nwub)
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
  if vdf.empty and not (args.add_exac or args.add_gnomad or args.add_1kg or args.add_pathogenic or args.add_cosmic or args.add_tcga):
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

  vdf_neutrals_deferred = vdf.groupby(["unp_pos","ref","alt","chain"],as_index=False).apply(defer_to_pathogenic)
  
  # Construct a dataframe from the coordinate file
  sdf = pd.DataFrame([[r.get_parent().id,r.id] for r in s.get_residues()],columns=["chain","pdb_pos"])
  sdf["unp_pos"] = [r.get_parent().alignment.pdb2seq[r.id] for r in s.get_residues()]

  # Calculate the coordinate center-of-mass for all residues
  coords   = [resicom(r) for r in s.get_residues()]
  coord_df = pd.DataFrame(coords,index=sdf.index,columns=["x","y","z"])
  sdf_coord = sdf.merge(coord_df,left_index=True,right_index=True)
   
  # Merge the sequence-aligned structure dataframe with the sequence-based variant dataframe
  vdf_extract = vdf_neutrals_deferred[["unp_pos","ref","alt","chain","dcode","qt"]].copy()
  vdf_extract['unp_pos'] = vdf_extract['unp_pos'].apply(int)
  vdf_extract['chain'] = vdf_extract['chain'].apply(str)
  sdf_coord['unp_pos'] = sdf_coord['unp_pos'].apply(int)
  sdf_coord['chain'] = sdf_coord['chain'].apply(str)
  sdf_merged = sdf_coord.merge(vdf_extract,how="left",on=["chain","unp_pos"])
  # Ensure that no duplicate residues were introduced by the merge
  sdf_merged.reset_index(drop=True,inplace=True)
  sdf = sdf_merged.drop_duplicates(["unp_pos","pdb_pos","chain","ref","alt","dcode"]).reset_index(drop=True)

  # Annotate each amino acid substitution with the blosum100 score
  sdf["blosum100"] = sdf.apply(lambda x: blosum100[(x["ref"],x["alt"])] if not np.isnan(x["dcode"]) else None,axis=1)

  # Fail if all variants outside structural coverage or mapped to missing residues
  msg = None
  vdf = sdf[~sdf["dcode"].isnull()]
  if vdf.empty:
    if not sdf.empty:
      msg = "Structure successfully aligned to sequence, but no variants were mapped to non-missing residues."
    else:
      msg = "ERROR: Sequence-structure alignment failed."
    LOGGER.critical(msg); 
    sys_exit_failure(msg)


  # Check that both variant categories are populated (if input is categorical)
  if not args.quantitative:
    msg = None
    if (sdf["dcode"]==0).sum() < 3:
      msg = "\nWARNING: Structure contains %d neutral variants (PathProx minimum 3).\n"%(sdf["dcode"]==0).sum()
      if not args.add_exac and not args.add_gnomad and not args.add_1kg:
        msg += "Consider using --add_exac, --add_gnomad, or --add_1kg.\n"
      else:
        msg += "Please manually specify neutral variants.\n"
    if msg:
      LOGGER.warning(msg)
    if (sdf["dcode"]==1).sum() < 3:
      msg = "\nWARNING: Structure contains %d pathogenic variants (PathProx minimum 3).\n"%(sdf["dcode"]==1).sum()
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
  cnts   = sdf.drop_duplicates(["unp_pos","ref","alt","dcode"]).groupby("dcode").apply(lambda x: (x["dcode"].values[0],len(x)))
  LOGGER.info("################################")
  LOGGER.info("Unique Variant Counts:")
  with open("%s_variant_counts.txt"%args.label,'w') as fout:
    writer = csv.writer(fout,delimiter='\t')    
    for dcode,cnt in cnts[::-1]:
      writer.writerow([code2class[dcode],cnt])
      LOGGER.info("%20s:  %d"%(code2class[dcode],cnt))
  LOGGER.info("")
  return sdf

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
    
def Kest(D,y,T=[],P=9999):
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
      Remaining (9999 etc) columns are the permutations
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

def uniK_calc_and_plot(D,y,P=9999,label="",weights=False):
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

def biD(pathogenicCoordinates,neutralCoordinates,P=9999,label=""):
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

def ripley(sdf,permutations=9999):
  """ Handles Ripley's univariate K and bivariate D
      analyses using variant pathogenic and neutral
      labels and a constructed inter-residue distance matrix

      Parameters
      ----------
      sdf: Structural dataframe contining x/y/z values and a dcode for each residue
      """
  LOGGER.info("Begin ripley()  permutations=%d"%permutations)
  # Distance matrix for all residues, regardless of labeling
  D = squareform(pdist(sdf[['x','y','z']]))
  # Set the diagnoal of the distance matrix to nans
  np.fill_diagonal(D,np.nan)
  # old way to do above lineD[np.identity(D.shape[0],dtype=bool)] = np.nan
  # Label vectors for pathogenic and neutral variants
  p = (sdf["dcode"]==1).astype(int)
  n = (sdf["dcode"]==0).astype(int)
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
  pathogenic_coordinates = sdf.loc[sdf["dcode"]==1,['x','y','z']]
  neutral_coordinates = sdf.loc[sdf["dcode"]==0,['x','y','z']]

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

#=============================================================================#
if __name__ == "__main__":
  # Determine the script directory
  # Which is the directory that pathprox.py is running in plus "/bin"
  # A number of ancillary scripts are in that directory
  # This script will use pathvis.py to assist with chimera output
  script_dir = "%s/bin"%os.path.dirname(os.path.abspath(sys.argv[0]))
  ## Parse Command Line Options ##
  import os,sys,argparse,configparser
  from pdbmap import PDBMapModel
  from psb_shared import psb_config,psb_progress

  # Setup the Command Line Argument Parser
  cmdline_parser = psb_config.create_default_argument_parser(__doc__,".")

  # Input parameters
  cmdline_parser.add_argument("entity",type=str,
                      help="Gene ID, UniProt AC, PDB ID, or PDB filename")
  cmdline_parser.add_argument("refseq",type=str,
                      help="NM_.... refseq transcript ID")
  cmdline_parser.add_argument("variants",type=str,nargs='?',
                      help="Comma-separated list of protein HGVS or \
                            a filename containing one identifier per line")
  cmdline_parser.add_argument("--fasta",type=str,
                      help="UniProt fasta file for the reference sequence")
  cmdline_parser.add_argument("--pathogenic",type=str,
                      help="User-defined set of pathogenic variants")
  cmdline_parser.add_argument("--pathogenic_label",type=str,
                      help="Label for User-defined pathogenic variants")
  cmdline_parser.add_argument("--neutral",type=str,
                      help="User-defined set of neutral variants")
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
  cmdline_parser.add_argument("--add_benign",action="store_true",default=False,
                      help="Supplement neutral variant set with ClinVar benign missense variants")
  cmdline_parser.add_argument("--add_pathogenic",action="store_true",default=False,
                      help="Supplement pathogenic variant set with ClinVar pathogenic missense variants")
  cmdline_parser.add_argument("--add_drug",action="store_true",default=False,
                      help="Supplement pathogenic variant set with ClinVar drug response")
  cmdline_parser.add_argument("--add_cosmic",action="store_true",default=False,
                      help="Supplement pathogenic variant set with COSMIC somatic missense variants")
  cmdline_parser.add_argument("--add_tcga",action="store_true",default=False,
                      help="Supplement pathogenic variant set with TCGA somatic missense variants")

  # Filter parameters
  cmdline_parser.add_argument("--chain",type=str,
                      help="Limit the analysis to a particular chain")
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
  cmdline_parser.add_argument("--permutations",type=int,default=9999,
                      help="Number of permutations for Ripley's K/D analyses")

  # Output parameters
  cmdline_parser.add_argument("-o","--outdir",type=str,
                      help="Directory to use for output and results")
  cmdline_parser.add_argument("--uniquekey",type=str,required=True,
                      help="A gene/refseq/mutation/structure/chain/flavor unique identifer for this particular job")
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

  args = cmdline_parser.parse_args()
  if not args.outdir:
    if args.uniquekey:
       args.outdir = args.uniquekey
    else:
       args.outdir = 'results'
    timestamp = strftime("%Y-%m-%d")
    if not args.no_timestamp and timestamp not in args.outdir:
      args.outdir += "_%s"%timestamp
    LOGGER.info("The --outdir parameter is missing.  Set to %s"%args.outdir)

  makedirs_capra_lab(args.outdir,'Main outdir creation')
       
  # Write current version of this script to output directory
  # pw_name = pwd.getpwuid( os.getuid() ).pw_name # example jsheehaj or mothcw
  log_filename = "%s/%s.log"%(args.outdir,args.uniquekey)
  # cache_dir = "/tmp/cache"
  # if not os.path.exists(cache_dir):
  #  os.makedirs(cache_dir)

  LOGGER.info("Log file is %s"%log_filename)
  needRoll = os.path.isfile(log_filename)

  fh = RotatingFileHandler(log_filename, maxBytes=(1048576*5), backupCount=7)
  # formatter = logging.Formatter('%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s')
  fh.setFormatter(PATHPROX_LOG_FORMATTER)
  fh.setLevel(logging.INFO)
  LOGGER.addHandler(fh)

  if needRoll:
    fh.doRollover()

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
    sys_exit_failure(info)

  def statusdir_info(info):
    __info_update(info)
    new_progress_filename = os.path.join(statusdir,"progress.new")
    with open(new_progress_filename,'w') as f:
      f.write("%s: %s\n"%(__file__,inspect.currentframe().f_back.f_lineno))
    os.rename(new_progress_filename,'%s/progress'%statusdir)

  statusdir_info('Begun')

  required_config_items = ["dbhost","dbname","dbuser","dbpass",
    "collaboration",
    "pdb_dir",
    "sec2prim",
    "chimera_headless",
    "collaboration",
    "idmapping",
    "interpro_dir",
    "modbase2013_dir",
    "modbase2013_summary",
    "modbase2016_dir",
    "modbase2016_summary",
    "output_rootdir",
    "sprot",
    "swiss_dir",
    "swiss_summary"]

  config,config_dict = psb_config.read_config_files(args,required_config_items)

  config_dict_reduced = {x:config_dict[x] for x in required_config_items}

  config_dict = config_dict_reduced

  config_dict_shroud_password = {x:config_dict[x] for x in required_config_items}
  dbpass = config_dict.get('dbpass','?')
  config_dict_shroud_password['dbpass'] = '*' * len(dbpass)

  # pathprox_config_dict = dict(config.items('PathProx'))


  # LOGGER.info("Command Line Arguments")
  # pprint.pprint(vars(args))
  LOGGER.info("Command Line Arguments:\n%s"%pprint.pformat(vars(args)))
  LOGGER.info("Configuration File parameters:\n%s"%pprint.pformat(config_dict_shroud_password))
  # LOGGER.info("PathProx parameters:\n%s"%pprint.pformat(pathprox_config_dict))

  if args.sqlcache:
    cache_dir = args.sqlcache
  else:
    cache_dir = os.path.join(config_dict['output_rootdir'],config_dict['collaboration'],"sqlcache")

  # We must take care with the sql_cache because several processes can simultaneously attempt to make it!
  if not os.path.exists(cache_dir):
    try:
      os.makedirs(cache_dir)
    except:
      pass
    try:
      assert os.path.exists(cache_dir)
    except:
      sys_exit_failure("Fatal: Can't create sqlcache at destination path %s" %cache_dir)
  set_capra_group_sticky(cache_dir)

  LOGGER.info("SQL cache directory: %s"%cache_dir)

  ## Parameter validity checks
  if not args.label:
    # Generate label from input parameters
    # Strip file extensions
    if os.path.exists(args.entity):
      idx = -2 if args.entity[-2:]=='gz' else -1
      args.label = '_'.join(args.entity.split('/')[-1].split('.')[:idx])
    else:
      args.label = args.entity
  # Append relevant information to the label
  if args.chain:
    args.label += "_%s"%args.chain
  if args.add_pathogenic:
    args.label += "_clinvar"
  if args.add_cosmic:
    args.label += "_cosmic"
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
  if args.isoform and not args.fasta:
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
          (args.add_gnomad or args.add_pathogenic or args.add_cosmic or \
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

  # Prepare a database connection to load database variant sets, in case requested later
  io = PDBMapIO(config_dict['dbhost'],config_dict['dbuser'],config_dict['dbpass'],config_dict['dbname'])

  # Init PDBMap by load the summary information dictionaries that it will need to use:
  LOGGER.info("Initializing PDBMAP by loading idmapping sec2prom, sprot, modbase, and swiss model meta dictionaries")
  PDBMapProtein.load_idmapping(config_dict['idmapping'])
  PDBMapProtein.load_sec2prim(config_dict['sec2prim'])
  PDBMapProtein.load_sprot(config_dict['sprot'])
  statusdir_info('Initialized')

  unp_flag,olabel = False,args.label
  qflag  = args.quantitative
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
    # import pdb; pdb.set_trace()
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
    sdf   = var2coord(s,p,n,c,q)

    # On exit from above
    # sdf["dcode"] is 0 for neutral variants, 1 for pathogenic, -1 for our target variant(s) of interest
    # and NaN for the structural entries with no assignments

    # Check if there are enough variants to calculate neutral/pathogenic constraint
    nflag = (sdf["dcode"]==0).sum() > 2
    pflag = (sdf["dcode"]==1).sum() > 2

    pClosestPathogenicMutationPDB = None
    pNextPathogenicMutationPDB = None
    pClosestPathogenicMutationUNP = None
    pNextPathogenicMutationUNP = None
    pClosestPathogenicDistance = None
    pNextPathogenicSeqDelta = None

    mapped_variant_count =  (sdf["dcode"] == -1).sum()
    def format_resid(resid):
       # If that residue type first element of the biopython has anything, then
       # something weird is going on, so make some noise
       if len(resid[0].strip()):
         logger.warning("Malformed residue tuple: %s   Reformatting for readability/chimera.",str(resid))

       # If there is an insertion code then return 123F for residue id
       if len(resid[2].strip()):
          resid_formatted = "%d%s"%(resid[1],resid[2])
          logger.warning("Residue %s has insertion code and will be reformatted as %s"%(resid,resid_formatted))
          return resid_formatted

       # The usual case, a pdb resdidue number with out a type, without an insertion code
       return str(resid[1])
     
    # Check that candidate variants did not map to missing residues
    if args.variants and mapped_variant_count < 1:
      LOGGER.info("All candidate variants were aligned to missing residues.")
      # Unset the variants argument for downstream processing
      args.variants = []
    elif mapped_variant_count == 1: # One variant is normal case for UDN
      df_variant = sdf[sdf["dcode"] == -1]      
      df_pathogenic = sdf[sdf["dcode"] == 1]      
      if len(df_pathogenic) < 1:
        LOGGER.info("No pathogenic mutation points to calculate nearest-to-variant.")
      else:
        #First, get the distance in 3-space
        variantCOM = df_variant[['x','y','z']].values
        pathogenicCOMs = df_pathogenic[['x','y','z']].values
        distances_to_pathogenic = cdist(variantCOM,pathogenicCOMs)
        # Get nearest pathogenic position
        min_subscript = np.argmin(distances_to_pathogenic[0])
        pClosestPathogenicDistance = distances_to_pathogenic[0][min_subscript]

        df_pathogenic_nearest = df_pathogenic.iloc[min_subscript]
        pClosestPathogenicMutationPDB = "%s%s%s"%(
          df_pathogenic_nearest['ref'],format_resid(df_pathogenic_nearest['pdb_pos']),df_pathogenic_nearest['alt'])
        pClosestPathogenicMutationUNP = "%s%d%s"%(
          df_pathogenic_nearest['ref'],int(df_pathogenic_nearest['unp_pos']),df_pathogenic_nearest['alt'])
 
        print("3-space Nearest is %d %s %s"%(min_subscript,pClosestPathogenicMutationPDB,pClosestPathogenicMutationUNP))
        #Second, get the distance in sequence space
        variant_unp_pos = df_variant[['unp_pos']].values
        pathogenic_unp_poss = df_pathogenic[['unp_pos']].values
        distances_to_pathogenic = cdist(variant_unp_pos,pathogenic_unp_poss)
        pNextPathogenicSeqDelta = int(distances_to_pathogenic[0][min_subscript])
        # Get nearest pathogenic position
        min_subscript = np.argmin(distances_to_pathogenic[0])
        df_pathogenic_nearest = df_pathogenic.iloc[min_subscript]
        pNextPathogenicMutationPDB = "%s%s%s"%(
          df_pathogenic_nearest['ref'],format_resid(df_pathogenic_nearest['pdb_pos']),df_pathogenic_nearest['alt'])
        pNextPathogenicMutationUNP = "%s%d%s"%(
          df_pathogenic_nearest['ref'],int(df_pathogenic_nearest['unp_pos']),df_pathogenic_nearest['alt'])
        print("Sequence Nearest is %d %s %s"%(min_subscript,pNextPathogenicMutationPDB,pNextPathogenicMutationUNP))
    else:
      LOGGER.info("%d variants were mapped.  Impossible to calculate nearest pathogenic residues."%mapped_variant_count)

    


    # Write the renumbered structure to the output directory

    basePDBfilename = "%s_%s_%s"%(args.label,sid,bio)
    renumberedPDBfilename = "%s_renum.pdb"%basePDBfilename
    renumberedPDBfilenameSS = "%s_renumSS.pdb"%basePDBfilename
    originalPDBfilename = "%s.pdb"%basePDBfilename

    LOGGER.info("Writing renumbered PDB to file %s"%renumberedPDBfilename)
    io.set_structure(s_renum)
    io.save(renumberedPDBfilename)

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


    # Write the original structure to the output directory
    LOGGER.info("Writing original PDB to file %s"%originalPDBfilename)
    io.set_structure(s)
    io.save(originalPDBfilename)

    ## Write pathogenic, neutral, candidate, and/or quantiative trait Chimera attribute files
    if args.variants and qflag:
      v = sdf[sdf["qt"].notnull()].sort_values(by=["chain","unp_pos"])
    elif nflag or pflag or args.variants:
      v = sdf[sdf["dcode"].notnull()].sort_values(by=["chain","unp_pos"])
    else:
      v = pd.DataFrame([])
    # Write attribute file with neutral variants
    if nflag:
      # Original PDB numbering
      with open("%s_neutral.attr"%args.label,'w') as fout:
        fout.write("attribute: neutral\n")
        fout.write("match mode: 1-to-1\n")
        fout.write("recipient: residues\n")
        for i,r in v.loc[v["dcode"]==0].iterrows():
          fout.write("\t:%s.%s\t%.3f\n"%(format_resid(r["pdb_pos"]),r["chain"],1.0))
      # Renumbered PDB
      with open("%s_renum_neutral.attr"%args.label,'w') as fout:
        fout.write("attribute: neutral\n")
        fout.write("match mode: 1-to-1\n")
        fout.write("recipient: residues\n")
        for i,r in v.loc[v["dcode"]==0].iterrows():
          fout.write("\t:%d.%s\t%.3f\n"%(r["unp_pos"],r["chain"],1.0))
    # Write attribute file with pathogenic variants
    if pflag:
      # Original PDB numbering
      with open("%s_pathogenic.attr"%args.label,'w') as fout:
        fout.write("attribute: pathogenic\n")
        fout.write("match mode: 1-to-1\n")
        fout.write("recipient: residues\n")
        for i,r in v.loc[v["dcode"]==1].iterrows():
          fout.write("\t:%s.%s\t%.3f\n"%(format_resid(r["pdb_pos"]),r["chain"],1.0))
      # Renumbered PDB
      with open("%s_renum_pathogenic.attr"%args.label,'w') as fout:
        fout.write("attribute: pathogenic\n")
        fout.write("match mode: 1-to-1\n")
        fout.write("recipient: residues\n")
        for i,r in v.loc[v["dcode"]==1].iterrows():
          fout.write("\t:%d.%s\t%.3f\n"%(r["unp_pos"],r["chain"],1.0))
    # Write attribute file with quantitative trait
    if qflag:
      with open("%s_quantitative.attr"%args.label,'w') as fout:
        fout.write("attribute: quantitative\n")
        fout.write("match mode: 1-to-1\n")
        fout.write("recipient: residues\n")
        for i,r in v.loc[v["qt"].notnull()].iterrows():
          fout.write("\t:%d.%s\t%.3f\n"%(r["unp_pos"],r["chain"],r["qt"]))
    # Label as 0.5 if both neutral and pathogenic variants are mapped to a residue
    filterwarnings('ignore')
    if not v.empty:
      vt          = v.drop_duplicates(["chain","pdb_pos"])
      vt["dcode"] = v.groupby(["chain","pdb_pos"]).apply(lambda x: np.mean(x["dcode"])).values
    else:
      vt = pd.DataFrame([])
    resetwarnings()
    # Write attribute file with all variants
    # Original PDB numbering
    with open("%s_variants.attr"%args.label,'w') as fout:
      fout.write("attribute: pathogenicity\n")
      fout.write("match mode: 1-to-1\n")
      fout.write("recipient: residues\n")
      for i,r in vt.iterrows():
        fout.write("\t:%s.%s\t%.3f\n"%(format_resid(r["pdb_pos"]),r["chain"],r["dcode"]))
    # Renumbered PDB
    with open("%s_renum_variants.attr"%args.label,'w') as fout:
      fout.write("attribute: pathogenicity\n")
      fout.write("match mode: 1-to-1\n")
      fout.write("recipient: residues\n")
      for i,r in vt.iterrows():
        fout.write("\t:%d.%s\t%.3f\n"%(r["unp_pos"],r["chain"],r["dcode"]))

    # If requested, run the univariate K and bivariate D analyses
    if nflag and pflag and args.ripley:
      LOGGER.info("Neutral and Pathogenic residue counts are sufficient to calculate both Ripley's univariate K and bivariate D...")
      pK,pKz,pKt,nK,nKz,nKt,D,Dz,Dt = ripley(sdf,args.permutations)
      qK = qKz = qKt = qKp = np.nan
    # If only one dataset present, process the univariate K only
    elif args.ripley and (nflag or pflag):
      # Distance matrix for all residues
      D = squareform(pdist(sdf[['x','y','z']]))
      D[np.identity(D.shape[0],dtype=bool)] = np.nan
      if nflag:
        n = (sdf["dcode"]==0).astype(int)
        nT,nK,nKz,nKt,nKp = uniK_calc_and_plot(D,n,args.permutations,"neutral")
        pK = pKz = pKt = np.nan
        qK = qKz = qKt = qKp = np.nan
      if pflag:
        p = (sdf["dcode"]==1).astype(int)
        pT,pK,pKz,pKt,pKp = uniK_calc_and_plot(D,p,args.permutations,"pathogenic")
        nK = nKz = nKt = np.nan
        qK = qKz = qKt = qKp = np.nan
      LOGGER.info("Insufficient variant counts to perform Ripley's D analysis.")
      LOGGER.info("Using default neighbor-weight parameters for PathProx.")
      args.radius = "NW"
      D = Dz = Dt = np.nan
    elif qflag and args.ripley:
      LOGGER.info("Calculating Ripley's univariate weighted K...")
      D = squareform(pdist(sdf[['x','y','z']]))
      D[np.identity(D.shape[0],dtype=bool)] = np.nan
      qT,qK,qKz,qKt,qKp = uniK_calc_and_plot(D,sdf['qt'],weights=True)
      pK = pKz = pKt = np.nan
      nK = nKz = nKt = np.nan
      D  = Dz  = Dt  = np.nan
    # If none of the datasets are populated, skip the Ripley analyses
    else:
      if args.ripley:
        LOGGER.info("Insufficient variant counts to perform Ripley's K/D analyses.")
        LOGGER.info("Using default neighbor-weight parameters for PathProx.")
      # Set all Ripley's results to NaN
      pK = pKz = pKt = nK = nKz = nKt = qK = qKz = qKt = D = Dz = Dt = np.nan
      args.radius = "NW"

    LOGGER.info("Phase 3 - Perform Pathprox cross-validation and report performance");
    # Run the PathProx cross-validation and report performance
    pathogenic_coordinates  = sdf.loc[sdf["dcode"]==1,['x','y','z']]     # Pathogenic
    neutral_coordinates  = sdf.loc[sdf["dcode"]==0,['x','y','z']]     # Neutral
    if not args.no_blosum:
      LOGGER.info("Integrating blosum-100 weights for all pathogenic and neutral residue");
      pathogenic_blosum100_weights = sdf.loc[sdf["dcode"]==1,"blosum100"]
      neutral_blosum100_weights = sdf.loc[sdf["dcode"]==0,"blosum100"]

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
    if nflag and pflag:
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
    all_coordinates  = sdf[['x','y','z']]
    if nflag and pflag:
      if args.no_blosum:
        LOGGER.info("Calculating unweighted PathProx scores (and z-scores)...")
        sdf["pathprox"] = pathprox(C,A,B,nwlb=nwlb,nwub=nwub)
      else:
        LOGGER.info("Calculating Blosum-100-weighted PathProx scores (and z-scores)...")
        all_blosum100_weights= sdf["blosum100"]
        sdf["pathprox"] = pathprox(all_coordinates,pathogenic_coordinates,neutral_coordinates,nwlb=nwlb,nwub=nwub,weights=(neutral_blosum100_weights,pathogenic_blosum100_weights,neutral_blosum100_weights))
    # Calculate neutral constraint scores for all residues
    if nflag:
      LOGGER.info("Calculating neutral constraint scores (and z-scores)...")
      # Now calculate for all residues
      if args.no_blosum:
        sdf["neutcon"] = 1-np.array(uniprox(all_coordinates,neutral_coordinates,nwlb=nwlb,nwub=nwub))
      else:
        all_blosum100_weights= sdf["blosum100"]
        sdf["neutcon"] = 1-np.array(uniprox(all_coordinates,neutral_coordinates,nwlb=nwlb,nwub=nwub,weights=(all_blosum100_weights,neutral_blosum100_weights)))
      # nneut = (sdf["dcode"]==0).sum()
      # pneutcon = [uniprox(sdf[['x','y','z']],
      #                     sdf.ix[np.random.choice(sdf.index,nneut),['x','y','z']],
      #                     nwlb=nwlb,nwub=nwub) for i in range(args.permutations)]
      # pneutcon = np.concatenate(([sdf['neutcon']],pneutcon))
      # sdf["neutcon_z"] = np.nan_to_num(zscore(pneutcon)[0])
    # Calculate pathogenic constraint scores for all residues
    if pflag:
      LOGGER.info("Calculating pathogenic constraint scores (and z-scores)...")
      # Now calculate for all residues
      if args.no_blosum:
        sdf["pathcon"] = uniprox(all_coordinates,pathogenic_coordinates,nwlb=nwlb,nwub=nwub)
      else:
        sdf["pathcon"] = uniprox(all_coordinates,pathogenic_coordinates,nwlb=nwlb,nwub=nwub)# ,w=(Cw,Aw)) DO NOT WEIGHT PATHOGENIC VARIANTS
      # npath = (sdf["dcode"]==1).sum()
      # ppathcon = [uniprox(sdf[['x','y','z']],
      #                     sdf.ix[np.random.choice(sdf.index,npath),['x','y','z']],
      #                     nwlb=nwlb,nwub=nwub) for i in range(args.permutations)]
      # ppathcon = np.concatenate(([sdf['pathcon']],ppathcon))
      # sdf["pathcon_z"] = np.nan_to_num(zscore(ppathcon)[0])
    # Calculate quantitative trait constraint scores for all residues
    if qflag:
      LOGGER.info("Calculating quantitative trait constraint scores (and z-scores)...")
      if args.radius in ("K","D"):
        nwlb,nwub = qKt,qKt
      sdf["qtprox"] = qtprox(sdf[['x','y','z']],sdf[['x','y','z','qt']],nwlb=nwlb,nwub=nwub)

    # Write scores to tab-delimited file
    # Neutral constraint
    if nflag:
      c = sdf.sort_values(by="neutcon").drop_duplicates(["chain","pdb_pos"])
      c = sdf.sort_values(by=["chain","pdb_pos"])
      c.to_csv("%s_neutcon.txt"%args.label,sep='\t',header=True,index=False)
    # Pathogenic constraint
    if pflag:
      c = sdf.sort_values(by="pathcon").drop_duplicates(["chain","pdb_pos"])
      c = sdf.sort_values(by=["chain","pdb_pos"])
      c.to_csv("%s_pathcon.txt"%args.label,sep='\t',header=True,index=False)
    # PathProx
    if nflag and pflag:
      c = sdf.sort_values(by="pathprox").drop_duplicates(["chain","pdb_pos"])
      c = sdf.sort_values(by=["chain","pdb_pos"])
      c.to_csv("%s_pathprox.txt"%args.label,sep='\t',header=True,index=False)
    # Quantative trait constraint
    if qflag:
      c = sdf.sort_values(by="qtprox").drop_duplicates(["chain","pdb_pos"])
      c = sdf.sort_values(by=["chain","pdb_pos"])
      c.to_csv("%s_qtprox_%s.txt"%(args.label,args.qlabel),sep='\t',header=True,index=False)

    # Write scores to Chimera attribute file
    # Neutral constraint
    if nflag:
      # Original PDB numbering
      with open("%s_neutcon.attr"%args.label,'w') as fout:
        fout.write("attribute: neutcon\n")
        fout.write("match mode: 1-to-1\n")
        fout.write("recipient: residues\n")
        for i,r in c.iterrows():
          fout.write("\t:%s.%s\t%.3f\n"%(format_resid(r["pdb_pos"]),r["chain"],r["neutcon"]))
      # Renumbered PDB
      with open("%s_renum_neutcon.attr"%args.label,'w') as fout:
        fout.write("attribute: neutcon\n")
        fout.write("match mode: 1-to-1\n")
        fout.write("recipient: residues\n")
        for i,r in c.iterrows():
          fout.write("\t:%d.%s\t%.3f\n"%(r["unp_pos"],r["chain"],r["neutcon"]))
    # Pathogenic constraint
    if pflag:
      # Original PDB numbering
      with open("%s_pathcon.attr"%args.label,'w') as fout:
        fout.write("attribute: pathcon\n")
        fout.write("match mode: 1-to-1\n")
        fout.write("recipient: residues\n")
        for i,r in c.iterrows():
          fout.write("\t:%s.%s\t%.3f\n"%(format_resid(r["pdb_pos"]),r["chain"],r["pathcon"]))
      # Renumbered PDB
      with open("%s_renum_pathcon.attr"%args.label,'w') as fout:
        fout.write("attribute: pathcon\n")
        fout.write("match mode: 1-to-1\n")
        fout.write("recipient: residues\n")
        for i,r in c.iterrows():
          fout.write("\t:%d.%s\t%.3f\n"%(r["unp_pos"],r["chain"],r["pathcon"]))
    # PathProx
    if nflag and pflag:
      # Original PDB numbering
      with open("%s_pathprox.attr"%args.label,'w') as fout:
        fout.write("attribute: pathprox\n")
        fout.write("match mode: 1-to-1\n")
        fout.write("recipient: residues\n")
        for i,r in c.iterrows():
          fout.write("\t:%s.%s\t%.3f\n"%(format_resid(r["pdb_pos"]),r["chain"],r["pathprox"]))
      # Renumbered PDB
      with open("%s_renum_pathprox.attr"%args.label,'w') as fout:
        fout.write("attribute: pathprox\n")
        fout.write("match mode: 1-to-1\n")
        fout.write("recipient: residues\n")
        for i,r in c.iterrows():
          fout.write("\t:%d.%s\t%.3f\n"%(r["unp_pos"],r["chain"],r["pathprox"]))
    # Quantitative constraint
    if qflag:
      # Original PDB numbering
      with open("%s_qtprox.attr"%args.label,'w') as fout:
        fout.write("attribute: qtprox\n")
        fout.write("match mode: 1-to-1\n")
        fout.write("recipient: residues\n")
        for i,r in c.iterrows():
          fout.write("\t:%d.%s\t%.3f\n"%(format_resid(r["pdb_pos"]),r["chain"],r["qtprox"]))
      # Renumbered PDB
      with open("%s_renum_qtprox.attr"%args.label,'w') as fout:
        fout.write("attribute: qtprox\n")
        fout.write("match mode: 1-to-1\n")
        fout.write("recipient: residues\n")
        for i,r in c.iterrows():
          fout.write("\t:%d.%s\t%.3f\n"%(r["unp_pos"],r["chain"],r["qtprox"]))  

    # Report PathProx/constraint scores for candidate variants
    if args.variants and qflag:
      LOGGER.info("Quantitative constraint scores for candidate missense variants:")
      LOGGER.info(sdf.loc[sdf["dcode"]<0,["unp_pos","ref","alt","qtprox"]].sort_values( \
            by=["qtprox","unp_pos"],ascending=[False,True]).groupby(["unp_pos"]).apply(np.mean).to_string(index=False))
    elif args.variants:
      if nflag:
        LOGGER.info("Neutral constraint scores for candidate missense variants:")
        LOGGER.info(sdf.loc[sdf["dcode"]<0,["unp_pos","ref","alt","neutcon"]].sort_values( \
            by=["neutcon","unp_pos"],ascending=[False,True]).groupby(["unp_pos"]).apply(np.mean).to_string(index=False))
      if pflag:
        LOGGER.info("Pathogenic constraint scores for candidate missense variants:")
        LOGGER.info(sdf.loc[sdf["dcode"]<0,["unp_pos","ref","alt","pathcon"]].sort_values( \
            by=["pathcon","unp_pos"],ascending=[False,True]).groupby(["unp_pos"]).apply(np.mean).to_string(index=False))
        if nflag: # meaning, we had BOTH enough pathogenic and enough neutrals
          LOGGER.info("PathProx scores for candidate missense variants:")
          LOGGER.info(sdf.loc[sdf["dcode"]<0,["unp_pos","ref","alt","pathprox"]].sort_values( \
            by=["pathprox","unp_pos"],ascending=[False,True]).groupby(["unp_pos"]).apply(np.mean).to_string(index=False))

    # Write summary results to file
    # Ripley's K/D results
    head = ["Kz_path","Kt_path","Kz_neut","Kt_neut","Dz","Dt","Kz_quant","Kt_quant"]
    vals = [pKz,pKt,nKz,nKt,Dz,Dt,qKz,qKt]
    # PathProx Cross Validation, ROC and PR
    head.extend(["roc_auc","pr_auc"])
    vals.extend([roc_auc,pr_auc])
    # Extract VUS names
    if (sdf["dcode"]<0).sum()>0:
      muts = sdf.loc[sdf["dcode"]<0,["ref","unp_pos","alt"]]
      muts = muts["ref"].str.cat([muts["unp_pos"].astype(str),muts["alt"]])
      # PathProx scores
      if pflag and nflag:
          pp_vus = sdf.loc[sdf["dcode"]<0,["unp_pos","ref","alt","pathprox"]].sort_values( \
                by=["pathprox","unp_pos"],ascending=[False,True]).groupby(["unp_pos"]).apply(np.mean)
          head.extend(["%s_pathprox"%m for m in muts])
          vals.extend(list(pp_vus["pathprox"].values))
      # Pathogenic constraint scores
      if pflag:
          pc_vus = sdf.loc[sdf["dcode"]<0,["unp_pos","ref","alt","pathcon"]].sort_values( \
                by=["pathcon","unp_pos"],ascending=[False,True]).groupby(["unp_pos"]).apply(np.mean)
          head.extend(["%s_pathcon"%m for m in muts])
          vals.extend(list(pc_vus["pathcon"].values))
      # Neutral constraint scores
      if nflag:
          nc_vus = sdf.loc[sdf["dcode"]<0,["unp_pos","ref","alt","neutcon"]].sort_values( \
                by=["neutcon","unp_pos"],ascending=[False,True]).groupby(["unp_pos"]).apply(np.mean)
          head.extend(["%s_neutcon"%m for m in muts])
          vals.extend(list(nc_vus["neutcon"].values))

    # Logically if any of these "nearest pathogenic residue" calculations are set, then the all should be
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
    params = [args.label,sid,str(bio)]
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

    variant_residues = []
    neutral_residues = []
    pathogenic_residues = []
    attribute_residues = []
    for i,r in sdf.iterrows():
      residue_tuple = (r["unp_pos"],r["chain"])
      if r['dcode'] == 0:
        neutral_residues.append(residue_tuple)
      elif r['dcode'] == 1:
        pathogenic_residues.append(residue_tuple)
      elif r['dcode'] == -1:
        variant_residues.append(residue_tuple)
      # If chimera-specific per-residue attributes are in memory, dump them too
      if (r['qt'] is not None) and (not np.isnan(r['qt'])):
        attribute_residues.append({'residue': residue_tuple, 'attribute': r['qt']})

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
