#!/usr/bin/env python3
# coding: utf-8
# Project        : Integrate PeSto into vustruct pipeline
# Filename       : run_pesto.py
# Authors        : Chris Moth
# Organization   : Meiler Lab
#                : Vanderbilt University
# Email          : chris.moth@vanderbilt.edu
# Date           : 2025-
#
# PeSTo is published at
# Krapp, L.F., Abriata, L.A., Cortés Rodriguez, F. et al. 
# PeSTo: parameter-free geometric deep learning for accurate 
# prediction of protein binding interfaces. 
# Nat Commun 14, 2175 (2023). https://doi.org/10.1038/s41467-023-37701-8
#
# The code below is lifted from the apply_model.ipynb Juyper Notebook at:
# https://github.com/LBM-EPFL/PeSTo
# 
# =============================================================================#

"""
Launch the inference module of PeSTo by supplying commad line settings as described
in the help
"""

import os
import sys
import logging

from psb_shared import psb_config

ch = logging.StreamHandler()
LOGGER = logging.getLogger()
LOGGER.addHandler(ch)
from logging.handlers import RotatingFileHandler


log_format_string = '%(asctime)s %(levelname)-4s [%(filename)16s:%(lineno)d] %(message)s'
date_format_string = '%H:%M:%S'
log_formatter = logging.Formatter(log_format_string, date_format_string)
ch.setFormatter(log_formatter)

LOGGER.setLevel(logging.DEBUG)
ch.setLevel(logging.INFO)

import argparse


parser = psb_config.create_default_argument_parser(__doc__, os.path.dirname(os.path.dirname(__file__)))

parser.add_argument('--data_dir', 
          help="path where input PDB files are located", 
          type=str,
          required = True)

parser.add_argument('--model_file', 
          help="Path to a specific trained model .pt file",
          type=str,
          default = "/opt/PeSTo/model/save/i_v4_1_2021-09-07_11-21/model_ckpt.pt")

parser.add_argument('--pytorch_device', 
          help="cpu or gpu",
          type=str,
          default = "gpu",
          required = True)


parser.add_argument("-o", "--outdir", type=str,
                            help="Directory to use for output and results")
parser.add_argument("--uniquekey", type=str, required=False,
                            help="A gene/refseq/mutation/structure/chain/flavor unique identifer for this pathprox job")


args = parser.parse_args()


# cache_dir = "/tmp/cache"
log_filename = os.path.join(args.outdir, args.uniquekey + ".log")

LOGGER.info("Log file is %s" % log_filename)
need_roll = os.path.isfile(log_filename)
LOGGER.info("Command: %s" % ' '.join(sys.argv))

if args.debug:
    LOGGER.setLevel(logging.DEBUG)
    ch.setLevel(logging.DEBUG)

logger_fh = RotatingFileHandler(log_filename, backupCount=5)

logger_fh.setFormatter(log_formatter)
logger_fh.setLevel(logging.DEBUG if args.debug else logging.INFO)
LOGGER.addHandler(logger_fh)

if need_roll:
    logger_fh.doRollover()




LOGGER.info(f"data_dir = {args.data_dir}")
LOGGER.info(f"model_file = {args.model_file}")
model_dir = os.path.dirname(args.model_file)
LOGGER.info(f"model_dir = {model_dir}")
LOGGER.info(f"pytorch_dir = {model_dir}")





# In[ ]:


import sys
import h5py
import json
import numpy as np
import torch as pt
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from glob import glob

from src.dataset import StructuresDataset, collate_batch_features, select_by_sid, select_by_interface_types
from src.data_encoding import encode_structure, encode_features, extract_topology, categ_to_resnames, resname_to_categ
from src.structure import data_to_structure, encode_bfactor, concatenate_chains, split_by_chain
from src.structure_io import save_pdb, read_pdb
from src.scoring import bc_scoring, bc_score_names


# In[ ]:


# data parameters
# args.data_dir = "examples/issue_19_04_2023"
# args.data_dir = "alphatest"


# In[ ]:


# model parameters
# R3
#save_path = "model/save/i_v3_0_2021-05-27_14-27"  # 89
#save_path = "model/save/i_v3_1_2021-05-28_12-40"  # 90
# R4
#save_path = "model/save/i_v4_0_2021-09-07_11-20"  # 89
# save_path = "model/save/i_v4_1_2021-09-07_11-21"  # 91
# save_path = "./alphatest" 


# select saved model
# model_filepath = os.path.join(save_path, 'model_ckpt.pt')
#model_filepath = os.path.join(save_path, 'model.pt')


# In[ ]:


# add module to path
if model_dir not in sys.path:
    sys.path.insert(0, model_dir)
    
# load functions
from config import config_model, config_data
from data_handler import Dataset
from model import Model


# In[ ]:


# define device

# create model
model = Model(config_model)

model.load_state_dict(pt.load(args.model_file, map_location=pt.device(args.pytorch_device)))

# set model to inference
model = model.eval().to(args.pytorch_device)


# In[ ]:


# find pdb files and ignore already predicted oins
pdb_filepaths = glob(os.path.join(args.data_dir, "*.pdb"), recursive=True)
LOGGER.info("Processing the following pdb files %s" % pdb_filepaths)

pdb_filepaths = [fp for fp in pdb_filepaths if "_i" not in fp]

# create dataset loader with preprocessing
dataset = StructuresDataset(pdb_filepaths, with_preprocessing=True)

# debug print
LOGGER.debug("Dataset len = %d", len(dataset))


# In[ ]:


# run model on all subunits
with pt.no_grad():
    for subunits, filepath in tqdm(dataset):
        # concatenate all chains together
        structure = concatenate_chains(subunits)

        # encode structure and features
        X, M = encode_structure(structure)
        #q = pt.cat(encode_features(structure), dim=1)
        q = encode_features(structure)[0]

        # extract topology
        ids_topk, _, _, _, _ = extract_topology(X, 64)

        # pack data and setup sink (IMPORTANT)
        X, ids_topk, q, M = collate_batch_features([[X, ids_topk, q, M]])

        # run model
        z = model(X.to(args.pytorch_device), ids_topk.to(args.pytorch_device), q.to(args.pytorch_device), M.float().to(args.pytorch_device))

        # for all predictions
        for i in range(z.shape[1]):
            # prediction
            p = pt.sigmoid(z[:,i])

            # encode result
            numpy_result = p.cpu().numpy()
            with open('/tmp/output_%d' % i,'w') as f:
                for residue,interaction_probability in enumerate(numpy_result):
                    f.write("%f %f\n" % (residue,interaction_probability))

            # structure = encode_bfactor(structure, p.cpu().numpy())
            structure = encode_bfactor(structure, numpy_result)

            # save results
            output_filepath = filepath[:-4]+'_i{}.pdb'.format(i)
            LOGGER.info("Writing %s" % output_filepath)
            save_pdb(split_by_chain(structure), output_filepath)

