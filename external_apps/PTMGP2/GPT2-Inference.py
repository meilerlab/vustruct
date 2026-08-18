#!/usr/bin/env python
# coding: utf-8

# 2025 July.  Chris Moth modified the provided GPT2-Interence.ipynb notebook at
# https://github.com/pallucs/PTMGPT2 to create this batch mode .py
# script
# 
# This involved 
# jupyter nbconvert --to python
# Then integrating command line arguments
# 
# Important: 
# This script has no dependences on any PDBMap/VUStruct codes.  Thus the calling routines in VUStruct
# control the software
#
# See also:
# Shrestha, P., Kandel, J., Tayara, H. et al. 
# Post-translational modification prediction via prompt-based fine-tuning of 
# a GPT-2 model. Nat Commun 15, 6699 (2024). https://doi.org/10.1038/s41467-024-51071-9
#
# Requires a suitable nvidia GPU card
# Requires download of the various model path files (~5FB each)


#Import necessary libraries
import os
import torch
import evaluate
import numpy as np
import pandas as pd
import json
from torch.utils.data import DataLoader,Dataset
from transformers import DataCollatorForSeq2Seq
from transformers import AutoTokenizer, GPT2LMHeadModel,TrainingArguments, Trainer,GPT2Config
from sklearn.metrics import average_precision_score,matthews_corrcoef,f1_score, precision_score, recall_score, balanced_accuracy_score


# The following model directories are available
# and mjst be provided on program invocation
model_directories = [
    'Acetylation (K)',
    'Phosphorylation (Y)',
    'Methylation (R)',
    'Succinylation (K)',
    'Sumoylation (K)',
    'N-linked Glycosylation (N)',
    'Ubiquitination (K)',
    'O-linked Glycosylation (S,T)',
    'S-nitrosylation (C)',
    'Malonylation (K)',
    'Methylation (K)',
    'Phosphorylation (S,T)',
    'Glutathionylation (C)',
    'Glutarylation (K)',
    'Amidation (V)',
    'S-palmitoylation (C)',
    'Hydroxylation (P)',
    'Hydroxylation (K)',
    'Formylation (K)'
    ]

import logging
import argparse
import re

LOGGER=logging.getLogger();
logging.basicConfig(level=logging.INFO)
LOGGER.info('Starting GPT2-Inference.py')

parser = argparse.ArgumentParser(
                    prog='GPT2-Inference',
                    description='Call PTMGPT2 inference module to predict post translational modifications',
                    epilog='Integrated into VUStruct pipeline')

parser.add_argument('--model_root', required=True, type=str,
                    help='Explicit path to the parent directory containing all model directories')
parser.add_argument('--output_root', required=True, type=str,
                    help='Explicit path to the directory where all .json output files are written')

# parser.add_argument('--model_directory', required=True, type=str,
#                     help='Type of PTM to predict.  Must be one of:' + ','.join(model_directories))

parser.add_argument('--sequence', type=str,
                    help='String of amino acid letter codes on which to predict')
parser.add_argument('--fasta', type=str,
                    help='Fasta file of amino acid sequence(s) on which to predict')
parser.add_argument('--tokenizer_path', 
                    default = '/opt/PTMGPT2/Tokenizer/', # Default to the Tokenizer inside the container
                    help='Only change for special development cases')

args = parser.parse_args()

assert args.fasta or args.sequence,  "You must specify --sequence or --fasta"
assert not (args.fasta and args.sequence), "Do not specify both --sequence and --fasta"

# Do a sanity check to make sure our needed files are in place
# No longer needed since we're processing all directories every time
# assert args.model_directory in model_directories, f"'{args.model_directory}' is not one of: " + '\n  ' + '\n  '.join(model_directories)

def find_subsequences(sequence:str, chars:list, left=10, right=10):
    subsequences = []
    length = len(sequence)
    # Iterate through the sequence to find the character
    for i, c in enumerate(sequence):
        if c in chars:
            # Calculate the start and end indices for the subsequence
            start = max(0, i - left)  # Ensure start is not less than 0
            end = min(length, i + right + 1)  # Ensure end does not exceed the sequence length

            # Append the subsequence to the list
            subsequences.append({'Seq':sequence[start:end],
                                 'Pos':i+1,
                                 'text':f'<startoftext>SEQUENCE:{sequence[start:end]}\nLABEL:'
                                 })
    return subsequences


# In[3]:


def read_fasta(file_path):
    """
    Reads a FASTA file and returns a dictionary with sequence identifiers as keys
    and sequences as values.

    :param file_path: str, path to the FASTA file
    :return: dict, dictionary with sequence IDs as keys and sequences as values
    """
    sequences = {}
    sequence_id = None
    sequence_data = []

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence_id is not None:
                    sequences[sequence_id] = ''.join(sequence_data)
                sequence_id = line[1:]
                sequence_data = []
            else:
                sequence_data.append(line)

        # Add the last sequence
        if sequence_id is not None:
            sequences[sequence_id] = ''.join(sequence_data)

    return sequences


# ## Load model | Load tokenizer | Prediction API call

# In[4]:


def load_model(mdl_pth):
    """
    Loads a pre-trained GPT-2 model from the specified path.

    :param mdl_pth: str, path to the model directory.
    :return: GPT2LMHeadModel, the loaded GPT-2 model in evaluation mode on the CPU.
    """
    model_config = GPT2Config.from_pretrained(mdl_pth)
    model = GPT2LMHeadModel.from_pretrained(mdl_pth,config=model_config,ignore_mismatched_sizes=True)
    return model.cpu().eval()

def tokenize(sub_sequences,tokenizer):
    """
    Tokenizes the given subsequences using the specified tokenizer.

    :param sub_sequences: list of dicts, each containing a 'text' field with the subsequence to tokenize.
    :param tokenizer: AutoTokenizer, the tokenizer to use for tokenizing the subsequences.
    :return: dict, the tokenized subsequences with padding applied.
    """
    sub_sequences=[x['text'] for x in sub_sequences]
    encoded=tokenizer(sub_sequences,return_tensors='pt',padding='longest')
    return encoded

def inference(input_seq,tokenizer_pth,model_pth,chars:list):
    """
    Performs inference on the input sequence using a specified tokenizer and model, and extracts labels.

    :param input_seq: str, the input sequence to process.
    :param tokenizer_pth: str, path to the tokenizer directory.
    :param model_pth: str, path to the model directory.
    :param chars: list of str, characters to find subsequences for.
    :return: dict, a JSON-like dictionary containing the input sequence, model type, and labeled results.
    """
    tokenizer = AutoTokenizer.from_pretrained(tokenizer_pth,padding_side='left')
    model=load_model(model_pth)
    sub_sequences=find_subsequences(input_seq,chars=chars)
    inputs_encode=tokenize(sub_sequences=sub_sequences,tokenizer=tokenizer)
    predicted=model.generate(inputs_encode['input_ids'],attention_mask=inputs_encode['attention_mask'],do_sample=False,top_k=50,max_new_tokens=2,top_p=0.15,temperature=0.1,num_return_sequences=0,pad_token_id=50259)
    predicted_text=tokenizer.batch_decode(predicted,skip_special_tokens=True)
    predicted_labels=[x.split('LABEL:')[-1] for x in predicted_text]
    json_results={'Sequence':input_seq,
                'Type':model_pth,
                'Results':[]
                }
    for label,sub_seq in zip(predicted_labels,sub_sequences):
        json_results['Results'].append({sub_seq['Pos']:label})
    return json_results

fasta_sequences = []
if args.fasta:
    #Read Fasta File
    fasta_sequences=read_fasta(args.fasta) # '/opt/PTMGPT2/Data/sample.fasta')

for model_directory in model_directories:
    model_path = os.path.join(args.model_root,model_directory)

    required_model_files=['config.json', 'pytorch_model.bin']
    for required_model_file in required_model_files:
        assert os.path.isfile(os.path.join(model_path,required_model_file)), f"{required_model_file} file missing from {model_path}"


    # Parse the Res list from the model_directory (A,B) string at end
    p = re.compile('.* \((.*)\)$')
    m = p.match(model_directory)
    residue_letters = m.group(1).split(',')
    LOGGER.debug("Res = %s" % str(residue_letters))

    if args.fasta: # There could be multi sequences in one fasta - thought not int he VUStruct application
        for sequence_id, seq in fasta_sequences.items():
            # If we have a VUStruct uniprot ID in the fasta strings, then use that as the ID for json filename and ignore longer things
            if '|' in sequence_id:
                sequence_id = sequence_id.split('|')[0]

            LOGGER.info(f"inference(seq={seq[0:5]}..., tokenizer_path={args.tokenizer_path} , model_path={model_path},residue_letters={residue_letters})")
            result_dict=inference(seq, args.tokenizer_path ,model_path,residue_letters)
            output_fname = os.path.join(args.output_root,f"{sequence_id}_{model_directory}.json")
            LOGGER.info(f"Writing prediction results for {sequence_id} {model_path} written to {output_fname}")
            _Results = result_dict['Results']

            # The outputs come as list of {res: str} tuple-like dictionaries
            # So I use the next(iter... approach to get the first elements of the tuples
            positive_results = [next(iter(residue_result.keys())) for residue_result in _Results if 'POSITIVE'  ==  next(iter(residue_result.values()))]
            LOGGER.info("'POSITIVE' PTMs for %s: %s", model_directory,str(positive_results))
            with open(output_fname,'w') as json_f:
                json_f.write(json.dumps(result_dict))
    else: # The sequence was given on the command line
        LOGGER.info(f"inference(seq={ars.sequence[0:5]}..., tokenizer_path={args.tokenizer_path} , model_path={model_path},residue_letters={residue_letters})")
        result_dict = inference(args.sequence,args.tokenizer_path,model_path,residue_letters)
        output_fname = os.path.join(args.output_root,f"{model_directory}.json")
        LOGGER.info(f"Writing prediction results for {model_path} written to {output_fname}")
        positive_results = [next(iter(residue_result.keys())) for residue_result in _Results if 'POSITIVE'  ==  next(iter(residue_result.values()))]
        LOGGER.info("'POSITIVE' PTMs for %s: %s", model_directory,str(positive_results))
        with open(output_fname,'w') as json_f:
            json_f.write(json.dumps(result_dict))


LOGGER.info("End of modified GPT2-Interence.py")

