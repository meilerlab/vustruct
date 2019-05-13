#!/usr/bin/env python
#
# Project        : PDBMap
# Filename       : get_pdbs.py
# Authors        : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2017-01-22
# Description    : Identifies all PDB structures of a given protein, and  
#                : selects those with coverage of the given mutation.
#=============================================================================#
import pandas as pd
import numpy as np
from pdbmap import PDBMapProtein

def get_pdbs(io,unp,mut_pos=None,maybe_threshold=20):
    """Return a dataframe of all structures and models of a uniprot identifier, optionally labeling
    whether a supplied mutation falls within each structure or model

    :param connection: The pre-opened MySQL connection
    :type connection: MySQL connection
    :param mut_pos: The optional transcript mutation position, for calculating whether a structure covers the mutation
    :type mut_pos: int
    :param maybe_threshold: the distance outside the ends of a model or structure, within which the "Maybe" label should be assigned
    :type maybe_threshold: int
    :return: dataframe of found structures, with column "dist_to_boundary" containing "Yes", "No", or "Maybe"
    :rtype: pandas Dataframe
    """
    
    
    # an isoform, as in ABC123-*"""
    sql = """
    SELECT a.label,
       d.method as method,a.structid,a.chain,c.transcript,
       min(c.trans_seqid) as seq_start,
       max(c.trans_seqid) as seq_end,
       COUNT(distinct b.seqid) as nresidues,
       d.resolution, NULL as identity,
       NULL as template
    FROM Chain a
    STRAIGHT_JOIN Residue b
    ON a.label=b.label AND a.structid=b.structid AND a.chain=b.chain AND a.biounit=b.biounit AND a.model=b.model
    INNER JOIN Alignment c
    ON b.label=c.label and b.structid=c.structid and b.chain=c.chain and b.seqid=c.chain_seqid
    LEFT JOIN Structure d
    ON a.label=d.label AND a.structid=d.pdbid
   WHERE a.biounit=0 AND a.model=0 AND (a.unp=%(base_unp)s OR a.unp LIKE %(base_unp_dash_pct)s) AND a.label in ('pdb')
    GROUP BY a.label,a.structid,a.chain
    ORDER BY d.pdbid IS NOT NULL,nresidues DESC,d.resolution ASC
    """

    # Read the query results into pandas DataFrame
    base_unp = unp.upper().split('-')[0]
   # appears that mut_pos is not used by this sql query above
    df  = pd.read_sql(sql,io._con,params={'mut_pos':mut_pos,'base_unp':base_unp,'base_unp_dash_pct':base_unp+"-%%"})
    if df.empty:        # ChrisMoth
        return df    # ChrisMoth
    df["resolution"] = df["resolution"].astype(np.float64)
    df["identity"]   = df["identity"].astype(np.float64)


    # Measure the distance from the mutated residue (if specified) to the sequence edge
    def dist_to_boundary(x,pos):
        seq_start = x['seq_start']
        seq_end   = x['seq_end']
        if seq_start < pos < seq_end:
            return 0
        elif pos < seq_start:
            return seq_start - pos
        else:
            return pos - seq_end
    if mut_pos:
        df["dist_to_boundary"]  = df[['seq_start','seq_end']].apply(dist_to_boundary,axis=1,pos=mut_pos)
        # Can we perform an analyses for this mutation?
        # Yes:   Mutation is covered by the protein structure
        # Maybe: Mutation is within maybe_threshold amino acids of the sequence edge
        # No:    Mutation is distal from the structural coverage
        def analysis_status(d):
            if d == 0:
                return 'Yes'
            elif d < maybe_threshold:
                return 'Maybe'
            else:
                return 'No'
        df['analysis_possible'] = df['dist_to_boundary'].apply(analysis_status)
    else:
        df["dist_to_boundary"]  = 0
        df["analysis_possible"] = "Maybe"

    # Rename the columns for reporting
    df.columns = ["Label","Method","Structure ID","Chain","Transcript",
                                "Sequence Start","Sequence End","Residues","Resolution (PDB)",
                                "Sequence Identity","PDB Template",
                                "Distance to Boundary","Analysis Possible?"]

    # Eliminate all structures that do not match transcript
    enstlist = PDBMapProtein.unp2enst(unp);
    return df[df.Transcript.isin(enstlist)]

def get_modbase_swiss(io,unp,mut_pos=None,maybe_threshold=20):
    """Return a dataframe of all structures and models of a uniprot identifier, optionally labeling
    whether a supplied mutation falls within each structure or model

    :param connection: The pre-opened MySQL connection
    :type connection: MySQL connection
    :param mut_pos: The optional transcript mutation position, for calculating whether a structure covers the mutation
    :type mut_pos: int
    :param maybe_threshold: the distance outside the ends of a model or structure, within which the "Maybe" label should be assigned
    :type maybe_threshold: int
    :return: dataframe of found structures, with column "dist_to_boundary" containing "Yes", "No", or "Maybe"
    :rtype: pandas Dataframe
    """
    
    
    # an isoform, as in ABC123-*"""
    sql = """
    SELECT a.label,
                 IF(e.method IS NOT NULL,e.method,'Swiss') as method,a.structid,a.chain,act.transcript,
                    min(actr.trans_seqid) as seq_start,
                    max(actr.trans_seqid) as seq_end,
                    COUNT(distinct b.seqid) as nresidues,
                    NULL as resolution,IF(e.identity IS NOT NULL,e.identity,Swiss.str_id) AS identity,
                    IF(e.pdbid IS NOT NULL,UPPER(e.pdbid),Swiss.template) as template
    FROM Chain a
    STRAIGHT_JOIN Residue b
    ON a.label=b.label AND a.structid=b.structid AND a.chain=b.chain AND a.biounit=b.biounit AND a.model=b.model
    INNER JOIN AlignChainTranscript act
    ON b.label=act.label and b.structid=act.structid and b.chain=act.chain 
    LEFT JOIN AlignChainTranscriptResidue actr
    ON actr.al_id = act.al_id and b.seqid=actr.chain_res_num and b.icode=actr.chain_res_icode
    LEFT JOIN Model e
    ON a.label=e.label AND a.structid=e.modelid
    LEFT JOIN Swiss
    ON a.label=Swiss.label AND a.structid=Swiss.modelid
   WHERE a.biounit=0 AND a.model=0 AND (a.unp=%(base_unp)s OR a.unp LIKE %(base_unp_dash_pct)s) AND a.label in ('modbase','swiss')
    GROUP BY a.label,a.structid,a.chain
    ORDER BY nresidues DESC,identity DESC
    """

    # Read the query results into pandas DataFrame
    base_unp = unp.upper().split('-')[0]
   # appears that mut_pos is not used by this sql query above
    df  = pd.read_sql(sql,io._con,params={'mut_pos':mut_pos,'base_unp':base_unp,'base_unp_dash_pct':base_unp+"-%%"})
    if df.empty:        # ChrisMoth
        return df    # ChrisMoth
    df["resolution"] = df["resolution"].astype(np.float64)
    df["identity"]   = df["identity"].astype(np.float64)


    # Measure the distance from the mutated residue (if specified) to the sequence edge
    def dist_to_boundary(x,pos):
        seq_start = x['seq_start']
        seq_end   = x['seq_end']
        if seq_start < pos < seq_end:
            return 0
        elif pos < seq_start:
            return seq_start - pos
        else:
            return pos - seq_end
    if mut_pos:
        df["dist_to_boundary"]  = df[['seq_start','seq_end']].apply(dist_to_boundary,axis=1,pos=mut_pos)
        # Can we perform an analyses for this mutation?
        # Yes:   Mutation is covered by the protein structure
        # Maybe: Mutation is within maybe_threshold amino acids of the sequence edge
        # No:    Mutation is distal from the structural coverage
        def analysis_status(d):
            if d == 0:
                return 'Yes'
            elif d < maybe_threshold:
                return 'Maybe'
            else:
                return 'No'
        df['analysis_possible'] = df['dist_to_boundary'].apply(analysis_status)
    else:
        df["dist_to_boundary"]  = 0
        df["analysis_possible"] = "Maybe"

    # Rename the columns for reporting
    df.columns = ["Label","Method","Structure ID","Chain","Transcript",
                                "Sequence Start","Sequence End","Residues","Resolution (PDB)",
                                "Sequence Identity","PDB Template",
                                "Distance to Boundary","Analysis Possible?"]

    # Eliminate all structures that do not match transcript
    enstlist = PDBMapProtein.unp2enst(unp);
    return df[df.Transcript.isin(enstlist)]
