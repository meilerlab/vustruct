#!/usr/bin/env python2.7

ClinvarFilename = '/tmp/allClinvar.csv'  # Welcome to change
COSMICFilename = '/tmp/allCOSMIC.csv'    # Welcome to change

import os,sys
import pandas as pd
import re

df_COSMIC = pd.DataFrame()
df_Clinvar = pd.DataFrame()

rootdir = "."
# Search all the stuff in the pipelineV1 PathProx output area
for root, subFolders, files in os.walk(rootdir):
  for f in files:
    if "_exac_D_summary.csv" in f:  # Have we found a pathprox summary file
      filename = os.path.join(root,f)
      print filename
      df = pd.DataFrame.from_csv(filename,sep='\t')
      if len(df) != 1:
        print "Something profoundly wrong in %s"%filename
        sys.exit(1)
      # Fish out the struct id from the filename
      m = re.search("^\./.*/(.*)/PathProx.*$",filename)
      if m:
        df['Struct ID'] = m.group(1)
      else:
        print "%s does not seem to contain a struct id"
        sys.exit(1)

      # Fish out the refseq from the filename
      m = re.search("^\./.*_(NM_{0,1}.*|RefSeqNotFound_UsingGeneOnly.*)_[A-Z][0-9]+[A-Z]/.*$",filename)
      if m:
        df['refseq'] = m.group(1)
      else:
        print "%s does not seem to contain a refseq"
        sys.exit(1)

      # Fish out the mutation from the filename
      m = re.search("^\./.*_(NM_{0,1}.*|RefSeqNotFound_UsingGeneOnly)_([A-Z][0-9]+[A-Z])/.*$",filename)
      if m:
        mut = m.group(2)
        df['mutation'] = mut
        # We need to rename some columns in the dataframe - to get rid of the mutation-specific that Mike added to the column name - uugh
        df = df.rename(columns={"%s_pathprox"%mut: "pathprox", "%s_pathcon"%mut: "pathcon", "%s_neutcon"%mut: "neutcon"})
        # However... the command above only works for models or structures where the mutation is in the same position as the transcript
        # for a lot of pdb structures, the mutation could be in a different "pdb" reference spot - so figure that out as well
        for c in df.columns:
          m = re.search("^([A-Z][0-9]+[A-Z])_(pathprox|pathcon|neutcon)",c)
          if m:
            mut = m.group(1)
            # print "************** got a mutation %s in %s"%(mut,filename)
            df = df.rename(columns={"%s_pathprox"%mut: "pathprox", "%s_pathcon"%mut: "pathcon", "%s_neutcon"%mut: "neutcon"})
      else:
        print "%s does not seem to contain a mutation"
        sys.exit(1)

      if 'linvar_exac' in f:
        df_Clinvar = df_Clinvar.append(df)
      elif 'osmic_exac' in f:
        df_COSMIC = df_COSMIC.append(df)
      else:
        print "%s does not seem to be either cosmic of clinvar!"
        sys.exit(1)
    

df_Clinvar.sort_values('pathprox',ascending = False).to_csv(ClinvarFilename,sep='\t')
df_COSMIC.sort_values('pathprox',ascending = False).to_csv(COSMICFilename,sep='\t')

print 'That which ye seek is in %s and %s'%(ClinvarFilename,COSMICFilename)
      
  # sys.exit(0)
