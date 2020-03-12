#!/usr/bin/env python
# 2018 Feb 2.  Chris Moth updated
# udn_pipeline to use new logging schemes
# config files, and pipeline2 compatability
# outputs
# Description    : Runs ddG calculations, or general sequence analsis, depending on input parameters
#                : spatial position relative to known pathogenic and neutral
#                : variation in protein structure.
###############################################################################
"""Runs various sequence and structure based analyses. 
If no structure is provided, only runs sequence-based analyses. Structure and sequence
based mutations must be provided if the structure numbering differs from the sequence numbering.
"""
## Package Dependenecies ##
# Standard
import sys,os,shutil,gzip,grp,stat
import logging

from logging.handlers import RotatingFileHandler
from logging import handlers
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


sh = logging.StreamHandler()
log_formatter = logging.Formatter('%(asctime)s %(levelname)-4s [%(filename)16s:%(lineno)4d] %(message)s',datefmt='%H:%M:%S')
ch = logging.StreamHandler()
ch.setFormatter(log_formatter)

# For this main routine, logger will be root
LOGGER = logging.getLogger()
LOGGER.addHandler(ch)
LOGGER.setLevel(logging.INFO)


import argparse,configparser
import pprint
import datetime
import glob
import re
import subprocess
import time
from collections import OrderedDict
from shutil import copyfile
from timeit import default_timer

from pdbmap import PDBMapSwiss

import inspect # For status updates



# Removed April 19 2018 sys.path.insert(0, '/dors/meilerlab/home/sliwosgr/gregit/')

# sys.path.insert(0,'/dors/meilerlab/home/sliwosgr/pdb_analysis_apps/gregit/')
udnpipepath = os.path.dirname(os.path.realpath(__file__))
# sys.path.insert(0, udnpipepath)
from udn_prepare import cleaner,parse_mutation,check_mutation
from udn_pipeline_classes import ddg_monomer,interface_analyzer,dssp,ligands,uniprot,predict_ss

#Note, paths should NOT end with /
# Was DSSP_PATH config_dict['dssp_exe'] = '/dors/capra_lab/projects/psb_collab/psb_pipeline/data/dssp/dssp_local.exe'
# now config_dict['pdb_dir'] + /structures/all/pdb - see code WAS PDB_PATH = '/dors/capra_lab/data/rcsb/structures/all/pdb'
# Now derifed from config_dict[pdb_dir'] BIOUNIT_PATH = '/dors/capra_lab/data/rcsb/biounit/coordinates/all'
# Now config_dict['modbase2016_dir'] # WAS  MODEL_PATH = '/dors/capra_lab/data/modbase/H_sapiens_2016/Homo_sapiens_2016/model'
# Chris Moth removed - because our chains are not here: config_dict['modbase2013_dir'] = '/dors/capra_lab/data/modbase/H_sapiens_2013/models/'
# config_dict['modbase2013_dir'] = '/dors/capra_lab/data/modbase/ModBase_H_sapiens_2013_GRCh37.70.pep.all/models/model/'
# Chris Moth found them ^^^^ here
RES_PATH = '/dors/capra_lab/projects/psb_collab/psb_pipeline/data/pdb_cmpd_res.idx'
# now config_dict['modbase2016_summary'] # WAS MODEL_METRICS = '/dors/capra_lab/data/modbase/H_sapiens_2016/Homo_sapiens_2016.summary.txt'
# now config_dict['modbase2013_summary']MODEL_METRICS_2013 = '/dors/capra_lab/data/modbase/H_sapiens_2013/H_sapiens_2013_GRCh37.70.pep.all.summary.txt'
SS_SCRIPTS = os.path.join(udnpipepath,"/udn_helper_scripts")
RES_QUALITY_CUTOFF = 2.5

from psb_shared import psb_config
cmdline_parser = psb_config.create_default_argument_parser(__doc__,os.path.dirname(os.path.dirname(__file__)))

cmdline_parser.add_argument('--structure', type=str, 
  help='structure to analyze (pdbid). If not provided no structural analyses performed. Must be paired with chain. If structure is 4mer pdbid, pdb file and biounits will be gathered (high quality). If not, first modbase library will be checked (low quality). If still not found, local directory will be checked (low quality).',required=False,default=None)
cmdline_parser.add_argument('-f', '--structure_file', type=str, 
  help='structure file - especially helpful in case of finding swissmodel structures',required=False,default=None)
cmdline_parser.add_argument(
  '--mutation', type=str, help='structure-numbered mutations to be analyzed, may be comma separated or file listing multiple mutations. \
   Must be in the form Arg30Pro or A30P. Will be inferred from sequence-based mutations if not provided or ignored if no structure provided', required=False,          default=None)
cmdline_parser.add_argument('-t', '--transcript', type=str, 
  help='transcript-numbered mutations to be analyzed, may be comma separated or indicate filename of linewise list. \
Must be in the form Arg30Pro or A30P. If not provided, will skip sequence-based analyses.', required=False, default=None)
cmdline_parser.add_argument('--gene', type=str, 
  help='Gene name to be analyzed. Required', required=True)
cmdline_parser.add_argument('--project', type=str, 
  help='Project designation. Default = "Misc".', required=False,default="Misc")
cmdline_parser.add_argument('--patient', type=str, 
  help='Patient identifier. Default = put results directly in gene subdirectory', required=False,default=None)
cmdline_parser.add_argument('--chain', type=str, 
  help='Chain in structure to be analyzed. Must be paired with pdbid and a single character', required=False,default=None)   
cmdline_parser.add_argument('-r', '--rosetta', type=str, 
  help='Path to Rosetta applications. default = /dors/meilerlab/apps/rosetta/rosetta-3.7/main/source/bin', required=False,default='/dors/meilerlab/apps/rosetta/rosetta-3.7/main/source/bin')   
cmdline_parser.add_argument("--outdir",type=str,
                      help="Directory to use for output and results")
cmdline_parser.add_argument("--uniquekey",type=str,required=True,
                      help="A gene/refseq/mutation/structure/chain/flavor unique identifer for this particular job")
cmdline_parser.add_argument('--ddg', type=int, 
  help='Number of ddg_monomer iterations to run. If 0, ddg_monomer is skipped. Default = 50', required=False, default=50)
cmdline_parser.add_argument('--quality', action='store_true',default=False,required=False,
  help='Overall automatic assignment of high quality to high res pdbids and low quality to modbase models low res pdbids. True forces high quality. False forces low quality.')
cmdline_parser.add_argument('--test', action='store_true',default=False,required=False,
  help='For testing. Redirect all file output to test directory.')

args,remaining_argv = cmdline_parser.parse_known_args()
if not args.outdir:
  if args.uniquekey:
     args.outdir = args.uniquekey
  else:
     args.outdir = 'UDNresults'
  timestamp = strftime("%Y-%m-%d")
  if not args.no_timestamp and timestamp not in args.outdir:
    args.outdir += "_%s"%timestamp
  LOGGER.info("The --outdir parameter is missing.  Set to %s"%args.outdir)

if not os.path.exists(args.outdir):
  try:
    os.makedirs(args.outdir)
  except:
    LOGGER.critical("Unable to create outdir: %s"%args.outdir)
    sys.exit(1)

set_capra_group_sticky(args.outdir)
     
# Write current version of this script to output directory
# pw_name = pwd.getpwuid( os.getuid() ).pw_name # example jsheehaj or mothcw

print(" Need to populate this logic!")
print((__file__.rstrip('c')))

# Initiate logging in the output directory
log_filename =  os.path.join(args.outdir,args.uniquekey + ".log") 

LOGGER.info("Log file is %s"%log_filename)
needRoll = os.path.isfile(log_filename)

fh = RotatingFileHandler(log_filename, maxBytes=(1048576*5), backupCount=7)
# formatter = logging.Formatter('%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s')
# fh.setFormatter(log_formatter)
fh.setLevel(logging.INFO)
LOGGER.addHandler(fh)

if needRoll:
  fh.doRollover()

LOGGER.info("Command: %s"%' '.join(sys.argv))

statusdir =  os.path.abspath(os.path.join(args.outdir,"status"))
LOGGER.info("Job status directory: %s"%statusdir)
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
    LOGGER.exception("Unable to delete file %s in status directory"%file_path)
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
  LOGGER.critical("Terminating in sys_exit_failure(): %s"%info)

  sys.exit(info)

def statusdir_info(info):
  __write_info_progress(info,"%s: %s\n"%(__file__,inspect.currentframe().f_back.f_lineno))


statusdir_info('Begun')

def makedirs_capra_lab(DEST_PATH, module_name):
    if not os.path.exists(DEST_PATH):
        os.makedirs(DEST_PATH)
    try:
        assert os.path.exists(DEST_PATH)
    except:
       sys_exit_failure("Fatal: Module %s failed to create destination path %s" % (module_name,DEST_PATH))
    set_capra_group_sticky(DEST_PATH)

required_config_items = ["dbhost","dbname","dbuser",
    # "dbpass",
    "pdb_dir",
    "sec2prim",
    "chimera_headless",
    "collaboration",
    "interpro_dir",
    "modbase2013_dir",
    "modbase2013_summary",
    "modbase2016_dir",
    "modbase2016_summary",
    "dssp_exe",
    "sprot",
    "swiss_dir",
    "swiss_summary"]

config,config_dict = psb_config.read_config_files(args,required_config_items)

missingKeys = [name for name in required_config_items if name not in config_dict]

config_dict_reduced = {x:config_dict[x] for x in required_config_items}

config_dict = config_dict_reduced

config_dict_shroud_password = {x:config_dict[x] for x in required_config_items}
# dbpass = config_dict.get('dbpass','?')
# config_dict_shroud_password['dbpass'] = '*' * len(dbpass)


# LOGGER.info("Command Line Arguments")
# pprint.pprint(vars(args))
LOGGER.info("Command Line Arguments:\n%s"%pprint.pformat(vars(args)))
LOGGER.info("Configuration File parameters:\n%s"%pprint.pformat(config_dict_shroud_password))

if (len(missingKeys)):
  LOGGER.critical('Can\'t proceed without configuration file options set for: %s'%str(missingKeys))

def run_ddg(DEST_PATH, rosetta_pdb, selected_chain,processed_mutations,prefix,quality):
    LOGGER.info("run_ddg() called with parameters:\n%s"%pprint.pformat(locals()))

    makedirs_capra_lab(DEST_PATH,'run_ddg')

    copyfile(rosetta_pdb,os.path.join(DEST_PATH,os.path.basename(rosetta_pdb)))

    try:
        os.chdir(DEST_PATH)
    except OSError:
        sys_exit_failure("Fatal: failed to enter %s" % DEST_PATH)
    LOGGER.info("curwd = %s"%os.getcwd())

    LOGGER.info("Running ddg_monomer")
    LOGGER.info("All generated files and results will be saved in %s\n" % DEST_PATH)

    #Run ddg monomer
    ddg = ddg_monomer(os.path.basename(rosetta_pdb),selected_chain,processed_mutations,args.rosetta)
    ddg_outcome,ddg_results = ddg.run(args.ddg,squality=quality)
    if not ddg_outcome:
        LOGGER.info("ddg_monomer -FAILED- with message: %s" % ddg_results)
    else: 
        finaloutfile = "%s_ddg.results" % prefix
        with open(finaloutfile,'w') as outfile:
            outfile.write("\n".join(ddg_results))
            outfile.write("\n")
        LOGGER.info("ddg_monomer success: Final results printed to %s" % os.path.join(DEST_PATH,finaloutfile))

def run_interface(DEST_PATH,rosetta_pdb,selected_chain,processed_mutations,prefix):
    LOGGER.info("run_interface() called with parameters:\n%s"%pprint.pformat(locals()))

    makedirs_capra_lab(DEST_PATH,'run_interface')

    try:
        os.chdir(DEST_PATH)
    except OSError:
        sys_exit_failure("Fatal: failed to enter %s" % DEST_PATH)
    LOGGER.info("curwd = %s"%os.getcwd())
    copyfile(rosetta_pdb,os.path.join(DEST_PATH,os.path.basename(rosetta_pdb)))   

    # import pdb; pdb.set_trace()
    LOGGER.info("All generated files and results can be found in %s" % DEST_PATH)
    interface = interface_analyzer(os.path.basename(rosetta_pdb), selected_chain, processed_mutations, args.rosetta)
    interface_outcome, interface_results = interface.run()
    if not interface_outcome:
        LOGGER.info("interface_analyzer -FAILED- with message %s" % interface_results)
    else:
        finaloutfile = "%s_interface.results" % prefix
        if os.path.isfile(finaloutfile):
            noheader = True
        else:
            noheader = False
        with open(finaloutfile,'a') as outfile:
            if noheader:
                outfile.write("\n".join(interface_results[1:]))
            else:
                outfile.write("\n".join(interface_results))            
            outfile.write("\n")
        LOGGER.info("interface_analyzer success: Final results printed to %s" % os.path.join(DEST_PATH,finaloutfile))

def run_dssp(DEST_PATH,dssp_pdb,selected_chain,processed_mutations,prefix):
    LOGGER.info("run_dssp() called with parameters:\n%s"%pprint.pformat(locals()))

    makedirs_capra_lab(DEST_PATH,'run_dssp')

    try:
        os.chdir(DEST_PATH)
    except OSError:
        sys_exit_failure("Fatal: failed to enter %s" % DEST_PATH)
    copyfile(dssp_pdb,os.path.join(DEST_PATH,os.path.basename(dssp_pdb)))   

    LOGGER.info("Running DSSP Analyzer on %s" % os.path.basename(dssp_pdb))
    LOGGER.info("All generated files and results can be found in %s\n" % DEST_PATH)
    interface = dssp(os.path.basename(dssp_pdb), selected_chain, processed_mutations, config_dict['dssp_exe'])
    interface_outcome, interface_results = interface.run()
    if not interface_outcome:
        LOGGER.info("DSSP -FAILED- with message %s" % interface_results)
    else:
        finaloutfile = "%s_dssp.results" % prefix
        if os.path.isfile(finaloutfile):
            noheader = True
        else:
            noheader = False
        with open(finaloutfile,'a') as outfile:
            if noheader:
                outfile.write("\n".join(interface_results[1:]))
            else:
                outfile.write("\n".join(interface_results))            
            outfile.write("\n")
        LOGGER.info("DSSP success: Final results printed to %s" % os.path.join(DEST_PATH,finaloutfile))

def run_ligands(DEST_PATH,raw_pdb,selected_chain,processed_mutations,prefix):
    LOGGER.info("run_ligands() called with parameters:\n%s"%pprint.pformat(locals()))

    makedirs_capra_lab(DEST_PATH,'run_ligands')

    try:
        os.chdir(DEST_PATH)
    except OSError:
        sys_exit_failure("Fatal: failed to enter %s" % DEST_PATH)
    copyfile(raw_pdb,os.path.join(DEST_PATH,os.path.basename(raw_pdb)))
    ligdist = ligands(os.path.basename(raw_pdb),selected_chain,processed_mutations)
    ligdist_outcome,ligdist_results = ligdist.run()
    if not ligdist_outcome:
        LOGGER.info("Ligand distances -FAILED- with message %" % ligdist_results)
    else:
        finaloutfile = "%s_ligands.results" % prefix
        if os.path.isfile(finaloutfile):
            noheader = True
        else:
            noheader = False
        with open(finaloutfile,'a') as outfile:
            if noheader:
                outfile.write("\n".join(ligdist_results[1:]))
            else:
                outfile.write("\n".join(ligdist_results))
            outfile.write("\n")
        LOGGER.info("Ligand distances success: final results printed to %s" % os.path.join(DEST_PATH,finaloutfile))

def run_unp(DEST_PATH,gene,processed_mutations,prefix):
    LOGGER.info("run_unp() called with parameters:\n%s"%pprint.pformat(locals()))
    makedirs_capra_lab(DEST_PATH,'run_unp')

    save_cwd = os.getcwd()

    try:
        os.chdir(DEST_PATH)
    except OSError:
        sys_exit_failure("Fatal: failed to enter %s" % DEST_PATH)

    # May 2018 - removed this exit.  Caller can decide whether or not they want a rerun
    # if os.path.isfile("%s_unp.results" % prefix) and not args.test:
    #    LOGGER.info("Uniprot information already present for %s" % prefix)
    #    return
    unp = uniprot(gene,processed_mutations)
    unp_outcome,unp_results = unp.run()
    if not unp_outcome:
        LOGGER.info("Uniprot annotation -FAILED- with message: %s" % unp_results)
        finaloutfile = "%s_unp.results" % prefix
        with open(finaloutfile,'a') as outfile:
            outfile.write(unp_results)
            outfile.write("\n")
    else:
        finaloutfile = "%s_unp.results" % prefix
        with open(finaloutfile,'w') as outfile:
            outfile.write("\n".join(unp_results))
            outfile.write("\n")
        LOGGER.info("Uniprot annotation success: final results printed to %s" % os.path.join(DEST_PATH, finaloutfile))

    LOGGER.info("Running secondary structure prediction")
    fastafile = "%s.fasta" % gene
    if not os.path.isfile(fastafile):
        LOGGER.info("Fasta file %s not found, unable to run ss predictions" % fastafile)
        os.chdir(save_cwd)
        return
    ss = predict_ss(fastafile,processed_mutations,scripts_path = SS_SCRIPTS)
    ss_outcome,ss_results = ss.run()
    ssoutfile = "%s_ss.results" % prefix
    if not ss_outcome:
        LOGGER.info("SS prediction -FAILED- with message: %s" % ss_results)
        with open(ssoutfile,'a') as outfile:
            outfile.write(ss_results)
            outfile.write("\n")
    else:
        with open(ssoutfile,'w') as outfile:
            outfile.write("\n".join("\t".join(x) for x in ss_results))
            outfile.write("\n")
        LOGGER.info("SS prediction success: final results printed to %s/%s" % (DEST_PATH,ssoutfile))
        visoutfile = "%s_ss.png" % prefix
        if ss.make_pic(visoutfile):
            LOGGER.info("SS visualization created %s/%s" % (DEST_PATH,visoutfile))
        else:
            LOGGER.info("Failed to generate visual file")
    os.chdir(save_cwd)
           
def get_mutations(mutation):
    with open(mutation) as infile:
        raw = [x.strip() for x in infile.readlines() if x.strip()!='']
    return raw       

def get_resolution_from_RES_PATH(pdbid):
    LOGGER.info("get_resolution() locating resolution of pdbid: %s in file %s"%(pdbid,RES_PATH))
    with open(RES_PATH) as infile:
        for line in infile:
            pdb = line.split(';')[0].strip()
            if pdb.upper() == pdbid.upper(): # Found pdbid in the file
                try:
                    retval =  float(line.split(';')[1].strip())
                    LOGGER.info("Returning found resolution=%f"%retval)
                    return retval
                except:
                    LOGGER.info("Resolution not found in line %s after :  Returning -1.0"%line)
                    return -1.0
    LOGGER.info("%s not found in pdb resolution file. Assuming low quality. Returning -1.0" % pdbid)
    return -1.



def evaluate_Modbase(model,old=False):
    """Modbase models are evaluated for whether or not they are of high enough quality
    To run in ddg monomer.
    Quality metrics must all be true, including
    modpipe: modpipe quality score (cumulative overall score from modbase pipeline)
    ztsvmod: RMSD of model to template.
    #seqid: sequence id of model seq and template seq """

    LOGGER.info("evaluate_Modbase() called with parametersi:\n%s"%pprint.pformat(locals()))

    modpipe = True #Removed filter for low modpipe scores (seqid and rmsd to template should be enough)
    tsvmod = False
    seqid = False
    #If it's a 2013 model, get the scores from the summary file to avoid xz nonsense
    if old:
        found = False
        with open(config_dict['modbase2013_summary']) as infile:
            for line in infile:
                cl = line.split()
                cm = cl[1]
                if cm!=os.path.splitext(os.path.basename(model))[0]: continue
                
                
    with gzip.open(model,'rt') as infile:
        for line in infile:
            if line.startswith('REMARK 220 SEQUENCE IDENTITY'):
                try:
                    sid = float(line.strip().split()[4])
                except:
                    sid = 0
                if sid>=40:
                    seqid = True
                LOGGER.info("Sequence Identity: %f"%sid)

            elif line.startswith('REMARK 220 TSVMOD RMSD'):
                try:
                    rmsd = float(line.strip().split()[4])
                except:
                    rmsd = 1000
                if rmsd<=4:
                    tsvmod = True
                LOGGER.info("TSVMOD RMSD: %f"%rmsd)

            elif line.startswith('REMARK 220 MODPIPE QUALITY SCORE'):
                try:
                    qual = float(line.strip().split()[5])
                except:
                    qual = 0
                if qual>=1.1: 
                    modpipe = True
                LOGGER.info("Modpipe quality score %f"%qual)

    LOGGER.info("modpipe: %s, tsvmod: %s, seqid: %s"%(str(modpipe),str(tsvmod),str(seqid)))
    if modpipe and tsvmod and seqid:
        return True
    else:
        return False

def evaluate_Swiss(modelid):
  """Swiss models are assumed to be of reasonable fundamental quality.  We just do a check for sid >= 40 and move on"""
  metrics = PDBMapSwiss.load_REMARK3_metrics(modelid)
  if 'sid' in metrics:
    sid = float(metrics['sid']) 
    threshold = 40.0
    LOGGER.info("Swiss Model sid=%f Minimum Acceptable: %f"%(sid,threshold))
    return float(metrics['sid']) >= threshold
  LOGGER.info("Swiss Model %s lacks REMARK for sid (seq identity)"%modelid)
  return False
    
##### MAIN
# import pdb; pdb.set_trace()
#Process arguments and generate constants based on arguments
# Chris Moth replaced wtih config_dict['outdir'] wAS 
# ROOT_PATH = '/dors/capra_lab/projects/psb_collab/%s' % args.project if not test else '/dors/meilerlab/home/sliwosgr/psb_collab/%s' % args.roject
if bool(args.structure) != bool(args.chain):
    sys_exit_failure("Fatal: Structure and Chain must be both present or both absent")   
if args.chain is not None:
    try:
        assert len(args.chain) == 1
    except AssertionError:
        sys_exit_failure("Fatal: Bad chain selection, must be 1 character")
if args.mutation is None and args.transcript is None:
    sys_exit_failure("Fatal: must provide at least one args.mutation source")
if not args.patient:
    args.patient = ""
else:
    args.patient = "%s/" % args.patient   
if args.mutation is None:
    args.mutation = args.transcript

   
##Parse args.mutation input
if os.path.isfile(args.mutation):
    raw = get_args.mutations(args.mutation)
    LOGGER.info("%d mutations read from file %s",len(raw),args.mutation)
else:
    raw = args.mutation.split(",")
    LOGGER.info("%d mutation(s) read from command line",len(raw))

processed_mutations = []
processed_trans = []

#Check structure args.mutations for quality if they will be needed
if args.structure is not None:
    for item in raw:
        current_mut = parse_mutation(item)
        if current_mut is None:
            sys_exit_failure("Fatal: Unable to parse args.mutation input %s" % item)
        elif current_mut[0]==current_mut[2]:
            LOGGER.info("Warning: %s has identical start and end aa's, skipping" % current_mut)
            continue
        else:
            processed_mutations.append(current_mut)
    try:
        assert len(processed_mutations)>0
    except AssertionError:
        sys_exit_failure("Fatal: No relevant mutations to run")

LOGGER.info("%d processed mutations"%len(processed_mutations))

#Check sequence args.mutations for quality if they will be needed
if args.transcript is not None:
    if os.path.isfile(args.transcript):
        rawtrans = get_args.mutations(args.transcript)
    else:
        rawtrans = args.transcript.split(",")   
    for item in rawtrans:
        current_mut = parse_mutation(item)
        if current_mut is None:
            sys_exit_failure("Fatal: Unable to parse args.mutation input %s" % item)
        elif current_mut[0]==current_mut[2]:
            LOGGER.info("Warning: %s has identical start and end aa's, skipping" % current_mut)
            continue
        else:
            processed_trans.append(current_mut)
    logging.getLogger("%d transcripts will be analyzed"%len(processed_trans))

if len(processed_trans)==0:
    LOGGER.info("Skipping transcript analyses")
if args.structure is None:
    LOGGER.info("Skipping structure analyses")

#Get the sequence level stuff out of the way in case no args.structure is provided
if len(processed_trans)>0:
    seq_prefix = "%s_%s" % (args.gene,args.transcript)    
    sequence_path = os.path.abspath(os.path.join(args.outdir,"sequence"))
    # ABOVE Path used to also have ,args.patient,args.gene,seq_prefix)
    run_unp(sequence_path,args.gene,processed_trans,seq_prefix)
if args.structure==None:
  msg = "No structure supplied on command line.  Exiting with success ."
  statusdir_info(msg)
  # Mark this analysis as complete
  LOGGER.info(msg)
  open("%s/complete"%statusdir,'w').close()
  sys.exit(0)

isModbase = False # Model files are subjected to quality check before ddg
# Swiss model is special case because the filename has a GUID in it
# and swiss in filename...  Eliminate that one first!
isSwiss = False
pdbfile = None
#Check if args.structure is local file and set quality = command line and biounits to local file
#DDG and other routines require full rooted path name to the structure file
# import pdb; pdb.set_trace()
is_local_pdb = False
if os.path.isfile(args.structure):
  if args.quality is None:
    args.quality = False
  if os.path.isfile(os.path.abspath(args.structure)):
    pdbfile = os.path.abspath(args.structure)
  else:
    if args.structure.startswith("/"):
      pdbfile = args.structure
    else:
      pdbfile = os.path.normpath("%s/%s" % (cwd,args.structure))
  biounits = [pdbfile]
  isModbase = False
  LOGGER.info("%s found as pdbfile=%s"%(args.structure,pdbfile))
  is_local_pdb = True
else: # The normal case, which is that the structure must be tracked down in our rcsb, modbase13, modbase16, or swiss respositories
  PDBMapSwiss.load_swiss_INDEX_JSON(config_dict['swiss_dir'],config_dict['swiss_summary']);

  try:
     pdbfile = PDBMapSwiss.get_coord_file(args.structure)
  except:
    pass

  isSwiss = (pdbfile and len(pdbfile) > 1)
  if isSwiss:
    isModbase = False
    biounits = ''
    modelfile = ''
    modelfile2013 = ''
    LOGGER.info("Swissmodel %s found as pdbfile=%s"%(args.structure,pdbfile))
  else:
    pdbfile = os.path.join(config_dict['pdb_dir'],"structures","divided","pdb",args.structure.lower()[1:3],"pdb%s.ent.gz"%args.structure.lower())
    #Check if args.structure is a 'real` rcsb deposited entry and get biounits and set quality accordingly
    if os.path.isfile(pdbfile):
      if not args.quality:
        res = get_resolution_from_RES_PATH(args.structure)
        if res<=RES_QUALITY_CUTOFF and res >0:
            LOGGER.info("pdb is high resolution, setting quality high")
            args.quality = True
        else:
            LOGGER.info("pdb is low resolution, setting quality low")
            args.quality = False
  #    if quality is None:
  #        quality = True
      biounits = glob.glob("%s/%s.pdb?.gz" % (os.path.join(config_dict['pdb_dir'],"biounit/PDB/divided",args.structure.lower()[1:3]), args.structure.lower()))        
      isModbase = False
    else: # The file must be some kind of modbase file - or it simply is a typo  Proceed as if modbase file
      biounits = ''
      modelfile = os.path.join(config_dict['modbase2016_dir'],"%s.pdb.gz"%args.structure)
      # modelfile2013 = "%s/%s.pdb.xz"
      # ChriMoth changing to gz
      modelfile2013 = os.path.join(config_dict['modbase2013_dir'],"%s.pdb.gz"%args.structure)
      # import pdb; pdb.set_trace()
  
      #Lastly check if it's modbase file and decide whether to run ddg monomer based on model stats.
      if os.path.isfile(modelfile):
        isModbase = True
        pdbfile = modelfile
      elif os.path.isfile(modelfile2013):
        isModbase = True
        pdbfile = modelfile2013
      else:
        sys_exit_failure("Fatal: Could not locate structure %s anywhere in filesystem.  Tried Swiss and :\n%s\n%s\n%s\n%s\n"%(args.structure,args.structure,pdbfile,modelfile,modelfile2013))
  
if (isModbase):
  biounits = [pdbfile]
  runddg = evaluate_Modbase(pdbfile,True)
  if not runddg:
    LOGGER.info("Low quality modbase file %s, skipping ddg_monomer"%pdbfile)
    args.ddg = 0

elif (isSwiss):
  biounits = [pdbfile]
  runddg = evaluate_Swiss(args.structure)
  if not runddg:
    LOGGER.info("Low quality swiss model %s, skipping ddg_monomer"%pdbfile)
    args.ddg = 0


suffix = datetime.datetime.now().strftime("%d%m%y%H%M%S")       

# if pdbfile[-3:] == '.gz':
#     pdbtext = os.path.basename(os.path.splitext(pdbfile[:-3])[0])
# else:
#     pdbtext = os.path.basename(os.path.splitext(pdbfile)[0])
# print pdbfile
biounittext = []
for item in biounits:
    if item[-3:] == '.gz':
        biounittext.append(os.path.basename(os.path.splitext(item)[0]))
    else:
        biounittext.append(os.path.basename(item))

#Destination path conventions
if is_local_pdb and args.structure.endswith('.pdb'):
    struct_prefix = "%s_%s_%s" % (os.path.basename(args.structure)[0:-4],args.chain,args.mutation)
else:  
    struct_prefix = "%s_%s_%s" % (os.path.basename(args.structure),args.chain,args.mutation)
cleaned_pdb_path = os.path.abspath(os.path.join(args.outdir,"clean_pdbs"))
ddg_path = os.path.abspath(os.path.join(args.outdir,"ddg"))
interface_path = os.path.abspath(os.path.join(args.outdir,"interface"))
ligand_path = os.path.abspath(os.path.join(args.outdir,"ligands"))

#Clean pdbs for use in Rosetta/DSSP
makedirs_capra_lab(cleaned_pdb_path,'Clean pdb path for main')

try:
    assert os.path.exists(cleaned_pdb_path)
except AssertionError:
    sys_exit_failure("Fatal: failed to create destination path %s" % cleaned_pdb_path)

pdbfile_abspath = os.path.abspath(pdbfile)
if pdbfile_abspath != pdbfile:
  pdbfile = pdbfile_abspath
  LOGGER.info("Adding full absolute path.  pdbfile now %s"%pdbfile);
pdbfile_copy_for_cleaning = os.path.join(cleaned_pdb_path, os.path.basename(pdbfile))
copyfile(pdbfile,os.path.join(cleaned_pdb_path, os.path.basename(pdbfile)))
LOGGER.info("Successful copy (cp) of %s to %s",pdbfile,pdbfile_copy_for_cleaning)

try:
    os.chdir(cleaned_pdb_path)
except OSError:
    sys_exit_failure("Fatal: failed to enter %s" % cleaned_pdb_path)

LOGGER.info("In preparation for cleaning, curwd is now: %s"%os.getcwd())

LOGGER.info("Cleaning %s for use in Rosetta" % os.path.basename(pdbfile))
preparation = cleaner(os.path.basename(pdbfile),verbose=False,outpath=".")
rosetta_pdb = preparation.write(bbcheck=True)

bio_rosetta = []
bio_dssp = []
for x in range(len(biounits)):
    if biounits[x]!=pdbfile:
        copyfile(biounits[x],"%s/%s" % (cleaned_pdb_path, os.path.basename(biounits[x])))
    LOGGER.info("Cleaning %s for use in Rosetta" % os.path.basename(biounits[x]))
    xprep = cleaner(os.path.basename(biounits[x]),verbose=False,outpath=cleaned_pdb_path)
    bio_rosetta.append(xprep.write(True))
    LOGGER.info("Cleaning %s for use in DSSP" % os.path.basename(biounits[x]))
    bio_dssp.append(xprep.write())
    
with open(rosetta_pdb) as infile:
    rawpdb = infile.readlines()

#Verify that args.mutations match the given residue #'s and that each args.mutation has a unique residue #
try:
    assert len(rawpdb)>0
except AssertionError:
    sys_exit_failure("Fatal: prepared %s is empty" % rosetta_pdb)
processed_mutations = check_mutation(processed_mutations,rawpdb,args.chain)

if processed_mutations is None:
    sys_exit_failure("Fatal: chain %s not found in %s" % (args.chain,args.structure))
try:
    assert len(processed_mutations)>0
except AssertionError:
    sys_exit_failure("Fatal: No relevant mutations to run")

LOGGER.info("%d mutations have been processed and verified" % len(processed_mutations))


#Run ddg_monomer if not silenced
if args.ddg>0:
    run_ddg(ddg_path,rosetta_pdb,args.chain,processed_mutations,struct_prefix,args.quality)
    
else:
    LOGGER.info("Skipping ddg_monomer")
    
#Run ligands
run_ligands(ligand_path, pdbfile, args.chain, processed_mutations, struct_prefix)

#Run interface score and DSSP interface for each biounut
LOGGER.info("Running %d biounit interfaces" % len(biounits))

for unit in range(len(bio_rosetta)):
    run_interface(interface_path,bio_rosetta[unit],args.chain,processed_mutations,struct_prefix)
    
    run_dssp(interface_path,bio_dssp[unit],args.chain,processed_mutations,struct_prefix)
    

LOGGER.info("All biounit interfaces complete.")
statusdir_info("All biounit interfaces complete.")
# Mark this analysis as complete
open("%s/complete"%statusdir,'w').close()

LOGGER.info("All biounit interfaces complete.")
