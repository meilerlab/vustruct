#!/usr/bin/env python
import gzip
import re
import sys
import os
from sys import argv, stderr, stdout
from os import popen, system
from os.path import exists, basename

import logging

# In case we are called by a client that has not configured logging,
# NullHandler ensures that does not fail
class NullHandler(logging.Handler):
  def emit(self, record):
    pass

LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(NullHandler())

MODRES={
'0CS':'ALA',                    ##  0CS ALA  3-[(S)-HYDROPEROXYSULFINYL]-L-ALANINE
'1AB':'PRO',                    ##  1AB PRO  1,4-DIDEOXY-1,4-IMINO-D-ARABINITOL
'1LU':'LEU',                    ##  1LU LEU  4-METHYL-PENTANOIC ACID-2-OXYL GROUP
'1PA':'PHE',                    ##  1PA PHE  PHENYLMETHYLACETIC ACID ALANINE
'1TQ':'TRP',                    ##  1TQ TRP  6-(FORMYLAMINO)-7-HYDROXY-L-TRYPTOPHAN
'1TY':'TYR',                    ##  1TY TYR
'23F':'PHE',                    ##  23F PHE  (2Z)-2-AMINO-3-PHENYLACRYLIC ACID
'23S':'TRP',                    ##  23S TRP  MODIFIED TRYPTOPHAN
'2BU':'ALA',                    ##  2BU ADE
'2ML':'LEU',                    ##  2ML LEU  2-METHYLLEUCINE
'2MR':'ARG',                    ##  2MR ARG  N3, N4-DIMETHYLARGININE
'2MT':'PRO',                    ##  2MT PRO
'2OP':'ALA',                    ##  2OP (2S  2-HYDROXYPROPANAL
'2TY':'TYR',                    ##  2TY TYR
'32S':'TRP',                    ##  32S TRP  MODIFIED TRYPTOPHAN
'32T':'TRP',                    ##  32T TRP  MODIFIED TRYPTOPHAN
'3AH':'HIS',                    ##  3AH HIS
'3MD':'ASP',                    ##  3MD ASP  2S,3S-3-METHYLASPARTIC ACID
'3TY':'TYR',                    ##  3TY TYR  MODIFIED TYROSINE
'4DP':'TRP',                    ##  4DP TRP
'4F3':'ALA',                    ##  4F3 ALA  CYCLIZED
'4FB':'PRO',                    ##  4FB PRO  (4S)-4-FLUORO-L-PROLINE
'4FW':'TRP',                    ##  4FW TRP  4-FLUOROTRYPTOPHANE
'4HT':'TRP',                    ##  4HT TRP  4-HYDROXYTRYPTOPHAN
'4IN':'TRP',                    ##  4IN TRP  4-AMINO-L-TRYPTOPHAN
'4PH':'PHE',                    ##  4PH PHE  4-METHYL-L-PHENYLALANINE
'5CS':'CYS',                    ##  5CS CYS
'6CL':'LYS',                    ##  6CL LYS  6-CARBOXYLYSINE
'6CW':'TRP',                    ##  6CW TRP  6-CHLORO-L-TRYPTOPHAN
'A0A':'ASP',                    ##  A0A ASP  ASPARTYL-FORMYL MIXED ANHYDRIDE
'AA4':'ALA',                    ##  AA4 ALA  2-AMINO-5-HYDROXYPENTANOIC ACID
'AAR':'ARG',                    ##  AAR ARG  ARGININEAMIDE
'AB7':'GLU',                    ##  AB7 GLU  ALPHA-AMINOBUTYRIC ACID
'ABA':'ALA',                    ##  ABA ALA  ALPHA-AMINOBUTYRIC ACID
'ACB':'ASP',                    ##  ACB ASP  3-METHYL-ASPARTIC ACID
'ACL':'ARG',                    ##  ACL ARG  DEOXY-CHLOROMETHYL-ARGININE
'ACY':'GLY',                    ##  ACY GLY  POST-TRANSLATIONAL MODIFICATION
'AEI':'THR',                    ##  AEI THR  ACYLATED THR
'AFA':'ASN',                    ##  AFA ASN  N-[7-METHYL-OCT-2,4-DIENOYL]ASPARAGINE
'AGM':'ARG',                    ##  AGM ARG  4-METHYL-ARGININE
'AGT':'CYS',                    ##  AGT CYS  AGMATINE-CYSTEINE ADDUCT
'AHB':'ASN',                    ##  AHB ASN  BETA-HYDROXYASPARAGINE
'AHO':'ALA',                    ##  AHO ALA  N-ACETYL-N-HYDROXY-L-ORNITHINE
'AHP':'ALA',                    ##  AHP ALA  2-AMINO-HEPTANOIC ACID
'AIB':'ALA',                    ##  AIB ALA  ALPHA-AMINOISOBUTYRIC ACID
'AKL':'ASP',                    ##  AKL ASP  3-AMINO-5-CHLORO-4-OXOPENTANOIC ACID
'ALA':'ALA',                    ##  ALA ALA
'ALC':'ALA',                    ##  ALC ALA  2-AMINO-3-CYCLOHEXYL-PROPIONIC ACID
'ALG':'ARG',                    ##  ALG ARG  GUANIDINOBUTYRYL GROUP
'ALM':'ALA',                    ##  ALM ALA  1-METHYL-ALANINAL
'ALN':'ALA',                    ##  ALN ALA  NAPHTHALEN-2-YL-3-ALANINE
'ALO':'THR',                    ##  ALO THR  ALLO-THREONINE
'ALS':'ALA',                    ##  ALS ALA  2-AMINO-3-OXO-4-SULFO-BUTYRIC ACID
'ALT':'ALA',                    ##  ALT ALA  THIOALANINE
'ALY':'LYS',                    ##  ALY LYS  N(6)-ACETYLLYSINE
'AME':'MET',                    ##  AME MET  ACETYLATED METHIONINE
'AP7':'ALA',                    ##  AP7 ADE
'APH':'ALA',                    ##  APH ALA  P-AMIDINOPHENYL-3-ALANINE
'API':'LYS',                    ##  API LYS  2,6-DIAMINOPIMELIC ACID
'APK':'LYS',                    ##  APK LYS
'AR2':'ARG',                    ##  AR2 ARG  ARGINYL-BENZOTHIAZOLE-6-CARBOXYLIC ACID
'AR4':'GLU',                    ##  AR4 GLU
'ARG':'ARG',                    ##  ARG ARG
'ARM':'ARG',                    ##  ARM ARG  DEOXY-METHYL-ARGININE
'ARO':'ARG',                    ##  ARO ARG  C-GAMMA-HYDROXY ARGININE
'ASA':'ASP',                    ##  ASA ASP  ASPARTIC ALDEHYDE
'ASB':'ASP',                    ##  ASB ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
'ASI':'ASP',                    ##  ASI ASP  L-ISO-ASPARTATE
'ASK':'ASP',                    ##  ASK ASP  DEHYDROXYMETHYLASPARTIC ACID
'ASL':'ASP',                    ##  ASL ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
'ASN':'ASN',                    ##  ASN ASN
'ASP':'ASP',                    ##  ASP ASP
'AYA':'ALA',                    ##  AYA ALA  N-ACETYLALANINE
'AYG':'ALA',                    ##  AYG ALA
'AZK':'LYS',                    ##  AZK LYS  (2S)-2-AMINO-6-TRIAZANYLHEXAN-1-OL
'B2A':'ALA',                    ##  B2A ALA  ALANINE BORONIC ACID
'B2F':'PHE',                    ##  B2F PHE  PHENYLALANINE BORONIC ACID
'B2I':'ILE',                    ##  B2I ILE  ISOLEUCINE BORONIC ACID
'B2V':'VAL',                    ##  B2V VAL  VALINE BORONIC ACID
'B3A':'ALA',                    ##  B3A ALA  (3S)-3-AMINOBUTANOIC ACID
'B3D':'ASP',                    ##  B3D ASP  3-AMINOPENTANEDIOIC ACID
'B3E':'GLU',                    ##  B3E GLU  (3S)-3-AMINOHEXANEDIOIC ACID
'B3K':'LYS',                    ##  B3K LYS  (3S)-3,7-DIAMINOHEPTANOIC ACID
'B3S':'SER',                    ##  B3S SER  (3R)-3-AMINO-4-HYDROXYBUTANOIC ACID
'B3X':'ASN',                    ##  B3X ASN  (3S)-3,5-DIAMINO-5-OXOPENTANOIC ACID
'B3Y':'TYR',                    ##  B3Y TYR
'BAL':'ALA',                    ##  BAL ALA  BETA-ALANINE
'BBC':'CYS',                    ##  BBC CYS
'BCS':'CYS',                    ##  BCS CYS  BENZYLCYSTEINE
'BCX':'CYS',                    ##  BCX CYS  BETA-3-CYSTEINE
'BFD':'ASP',                    ##  BFD ASP  ASPARTATE BERYLLIUM FLUORIDE
'BG1':'SER',                    ##  BG1 SER
'BHD':'ASP',                    ##  BHD ASP  BETA-HYDROXYASPARTIC ACID
'BIF':'PHE',                    ##  BIF PHE
'BLE':'LEU',                    ##  BLE LEU  LEUCINE BORONIC ACID
'BLY':'LYS',                    ##  BLY LYS  LYSINE BORONIC ACID
'BMT':'THR',                    ##  BMT THR
'BNN':'ALA',                    ##  BNN ALA  ACETYL-P-AMIDINOPHENYLALANINE
'BOR':'ARG',                    ##  BOR ARG
'BPE':'CYS',                    ##  BPE CYS
'BTR':'TRP',                    ##  BTR TRP  6-BROMO-TRYPTOPHAN
'BUC':'CYS',                    ##  BUC CYS  S,S-BUTYLTHIOCYSTEINE
'BUG':'LEU',                    ##  BUG LEU  TERT-LEUCYL AMINE
'C12':'ALA',                    ##  C12 ALA
'C1X':'LYS',                    ##  C1X LYS  MODIFIED LYSINE
'C3Y':'CYS',                    ##  C3Y CYS  MODIFIED CYSTEINE
'C5C':'CYS',                    ##  C5C CYS  S-CYCLOPENTYL THIOCYSTEINE
'C6C':'CYS',                    ##  C6C CYS  S-CYCLOHEXYL THIOCYSTEINE
'C99':'ALA',                    ##  C99 ALA
'CAB':'ALA',                    ##  CAB ALA  4-CARBOXY-4-AMINOBUTANAL
'CAF':'CYS',                    ##  CAF CYS  S-DIMETHYLARSINOYL-CYSTEINE
'CAS':'CYS',                    ##  CAS CYS  S-(DIMETHYLARSENIC)CYSTEINE
'CCS':'CYS',                    ##  CCS CYS  CARBOXYMETHYLATED CYSTEINE
'CGU':'GLU',                    ##  CGU GLU  CARBOXYLATION OF THE CG ATOM
'CH6':'ALA',                    ##  CH6 ALA
'CH7':'ALA',                    ##  CH7 ALA
'CHG':'GLY',                    ##  CHG GLY  CYCLOHEXYL GLYCINE
'CHP':'GLY',                    ##  CHP GLY  3-CHLORO-4-HYDROXYPHENYLGLYCINE
'CHS':'PHE',                    ##  CHS PHE  4-AMINO-5-CYCLOHEXYL-3-HYDROXY-PENTANOIC AC
'CIR':'ARG',                    ##  CIR ARG  CITRULLINE
'CLB':'ALA',                    ##  CLB ALA
'CLD':'ALA',                    ##  CLD ALA
'CLE':'LEU',                    ##  CLE LEU  LEUCINE AMIDE
'CLG':'LYS',                    ##  CLG LYS
'CLH':'LYS',                    ##  CLH LYS
'CLV':'ALA',                    ##  CLV ALA
'CME':'CYS',                    ##  CME CYS  MODIFIED CYSTEINE
'CML':'CYS',                    ##  CML CYS
'CMT':'CYS',                    ##  CMT CYS  O-METHYLCYSTEINE
'CQR':'ALA',                    ##  CQR ALA
'CR2':'ALA',                    ##  CR2 ALA  POST-TRANSLATIONAL MODIFICATION
'CR5':'ALA',                    ##  CR5 ALA
'CR7':'ALA',                    ##  CR7 ALA
'CR8':'ALA',                    ##  CR8 ALA
'CRK':'ALA',                    ##  CRK ALA
'CRO':'THR',                    ##  CRO THR  CYCLIZED
'CRQ':'TYR',                    ##  CRQ TYR
'CRW':'ALA',                    ##  CRW ALA
'CRX':'ALA',                    ##  CRX ALA
'CS1':'CYS',                    ##  CS1 CYS  S-(2-ANILINYL-SULFANYL)-CYSTEINE
'CS3':'CYS',                    ##  CS3 CYS
'CS4':'CYS',                    ##  CS4 CYS
'CSA':'CYS',                    ##  CSA CYS  S-ACETONYLCYSTEIN
'CSB':'CYS',                    ##  CSB CYS  CYS BOUND TO LEAD ION
'CSD':'CYS',                    ##  CSD CYS  3-SULFINOALANINE
'CSE':'CYS',                    ##  CSE CYS  SELENOCYSTEINE
'CSI':'ALA',                    ##  CSI ALA
'CSO':'CYS',                    ##  CSO CYS  INE S-HYDROXYCYSTEINE
'CSR':'CYS',                    ##  CSR CYS  S-ARSONOCYSTEINE
'CSS':'CYS',                    ##  CSS CYS  1,3-THIAZOLE-4-CARBOXYLIC ACID
'CSU':'CYS',                    ##  CSU CYS  CYSTEINE-S-SULFONIC ACID
'CSW':'CYS',                    ##  CSW CYS  CYSTEINE-S-DIOXIDE
'CSX':'CYS',                    ##  CSX CYS  OXOCYSTEINE
'CSY':'ALA',                    ##  CSY ALA  MODIFIED TYROSINE COMPLEX
'CSZ':'CYS',                    ##  CSZ CYS  S-SELANYL CYSTEINE
'CTH':'THR',                    ##  CTH THR  4-CHLOROTHREONINE
'CWR':'ALA',                    ##  CWR ALA
'CXM':'MET',                    ##  CXM MET  N-CARBOXYMETHIONINE
'CY0':'CYS',                    ##  CY0 CYS  MODIFIED CYSTEINE
'CY1':'CYS',                    ##  CY1 CYS  ACETAMIDOMETHYLCYSTEINE
'CY3':'CYS',                    ##  CY3 CYS  2-AMINO-3-MERCAPTO-PROPIONAMIDE
'CY4':'CYS',                    ##  CY4 CYS  S-BUTYRYL-CYSTEIN
'CY7':'CYS',                    ##  CY7 CYS  MODIFIED CYSTEINE
'CYD':'CYS',                    ##  CYD CYS
'CYF':'CYS',                    ##  CYF CYS  FLUORESCEIN LABELLED CYS380 (P14)
'CYG':'CYS',                    ##  CYG CYS
'CYJ':'LYS',                    ##  CYJ LYS  MODIFIED LYSINE
'CYQ':'CYS',                    ##  CYQ CYS
'CYR':'CYS',                    ##  CYR CYS
'CYS':'CYS',                    ##  CYS CYS
'CZ2':'CYS',                    ##  CZ2 CYS  S-(DIHYDROXYARSINO)CYSTEINE
'CZZ':'CYS',                    ##  CZZ CYS  THIARSAHYDROXY-CYSTEINE
'DA2':'ARG',                    ##  DA2 ARG  MODIFIED ARGININE
'DAB':'ALA',                    ##  DAB ALA  2,4-DIAMINOBUTYRIC ACID
'DAH':'PHE',                    ##  DAH PHE  3,4-DIHYDROXYDAHNYLALANINE
'DAL':'ALA',                    ##  DAL ALA  D-ALANINE
'DAM':'ALA',                    ##  DAM ALA  N-METHYL-ALPHA-BETA-DEHYDROALANINE
'DAR':'ARG',                    ##  DAR ARG  D-ARGININE
'DAS':'ASP',                    ##  DAS ASP  D-ASPARTIC ACID
'DBU':'ALA',                    ##  DBU ALA  (2E)-2-AMINOBUT-2-ENOIC ACID
'DBY':'TYR',                    ##  DBY TYR  3,5 DIBROMOTYROSINE
'DBZ':'ALA',                    ##  DBZ ALA  3-(BENZOYLAMINO)-L-ALANINE
'DCL':'LEU',                    ##  DCL LEU  2-AMINO-4-METHYL-PENTANYL GROUP
'DCY':'CYS',                    ##  DCY CYS  D-CYSTEINE
'DDE':'HIS',                    ##  DDE HIS
'DGL':'GLU',                    ##  DGL GLU  D-GLU
'DGN':'GLN',                    ##  DGN GLN  D-GLUTAMINE
'DHA':'ALA',                    ##  DHA ALA  2-AMINO-ACRYLIC ACID
'DHI':'HIS',                    ##  DHI HIS  D-HISTIDINE
'DHL':'SER',                    ##  DHL SER  POST-TRANSLATIONAL MODIFICATION
'DIL':'ILE',                    ##  DIL ILE  D-ISOLEUCINE
'DIV':'VAL',                    ##  DIV VAL  D-ISOVALINE
'DLE':'LEU',                    ##  DLE LEU  D-LEUCINE
'DLS':'LYS',                    ##  DLS LYS  DI-ACETYL-LYSINE
'DLY':'LYS',                    ##  DLY LYS  D-LYSINE
'DMH':'ASN',                    ##  DMH ASN  N4,N4-DIMETHYL-ASPARAGINE
'DMK':'ASP',                    ##  DMK ASP  DIMETHYL ASPARTIC ACID
'DNE':'LEU',                    ##  DNE LEU  D-NORLEUCINE
'DNG':'LEU',                    ##  DNG LEU  N-FORMYL-D-NORLEUCINE
'DNL':'LYS',                    ##  DNL LYS  6-AMINO-HEXANAL
'DNM':'LEU',                    ##  DNM LEU  D-N-METHYL NORLEUCINE
'DPH':'PHE',                    ##  DPH PHE  DEAMINO-METHYL-PHENYLALANINE
'DPL':'PRO',                    ##  DPL PRO  4-OXOPROLINE
'DPN':'PHE',                    ##  DPN PHE  D-CONFIGURATION
'DPP':'ALA',                    ##  DPP ALA  DIAMMINOPROPANOIC ACID
'DPQ':'TYR',                    ##  DPQ TYR  TYROSINE DERIVATIVE
'DPR':'PRO',                    ##  DPR PRO  D-PROLINE
'DSE':'SER',                    ##  DSE SER  D-SERINE N-METHYLATED
'DSG':'ASN',                    ##  DSG ASN  D-ASPARAGINE
'DSN':'SER',                    ##  DSN SER  D-SERINE
'DTH':'THR',                    ##  DTH THR  D-THREONINE
'DTR':'TRP',                    ##  DTR TRP  D-TRYPTOPHAN
'DTY':'TYR',                    ##  DTY TYR  D-TYROSINE
'DVA':'VAL',                    ##  DVA VAL  D-VALINE
'DYG':'ALA',                    ##  DYG ALA
'DYS':'CYS',                    ##  DYS CYS
'EFC':'CYS',                    ##  EFC CYS  S,S-(2-FLUOROETHYL)THIOCYSTEINE
'ESB':'TYR',                    ##  ESB TYR
'ESC':'MET',                    ##  ESC MET  2-AMINO-4-ETHYL SULFANYL BUTYRIC ACID
'FCL':'PHE',                    ##  FCL PHE  3-CHLORO-L-PHENYLALANINE
'FGL':'ALA',                    ##  FGL ALA  2-AMINOPROPANEDIOIC ACID
'FGP':'SER',                    ##  FGP SER
'FHL':'LYS',                    ##  FHL LYS  MODIFIED LYSINE
'FLE':'LEU',                    ##  FLE LEU  FUROYL-LEUCINE
'FLT':'TYR',                    ##  FLT TYR  FLUOROMALONYL TYROSINE
'FME':'MET',                    ##  FME MET  FORMYL-METHIONINE
'FOE':'CYS',                    ##  FOE CYS
'FOG':'PHE',                    ##  FOG PHE  PHENYLALANINOYL-[1-HYDROXY]-2-PROPYLENE
'FOR':'MET',                    ##  FOR MET
'FRF':'PHE',                    ##  FRF PHE  PHE FOLLOWED BY REDUCED PHE
'FTR':'TRP',                    ##  FTR TRP  FLUOROTRYPTOPHANE
'FTY':'TYR',                    ##  FTY TYR  DEOXY-DIFLUOROMETHELENE-PHOSPHOTYROSINE
'GHG':'GLN',                    ##  GHG GLN  GAMMA-HYDROXY-GLUTAMINE
'GHP':'GLY',                    ##  GHP GLY  4-HYDROXYPHENYLGLYCINE
'GL3':'GLY',                    ##  GL3 GLY  POST-TRANSLATIONAL MODIFICATION
'GLH':'GLN',                    ##  GLH GLN
'GLN':'GLN',                    ##  GLN GLN
'GLU':'GLU',                    ##  GLU GLU
'GLY':'GLY',                    ##  GLY GLY
'GLZ':'GLY',                    ##  GLZ GLY  AMINO-ACETALDEHYDE
'GMA':'GLU',                    ##  GMA GLU  1-AMIDO-GLUTAMIC ACID
'GMU':'ALA',                    ##  GMU 5MU
'GPL':'LYS',                    ##  GPL LYS  LYSINE GUANOSINE-5'-MONOPHOSPHATE
'GT9':'CYS',                    ##  GT9 CYS  SG ALKYLATED
'GVL':'SER',                    ##  GVL SER  SERINE MODIFED WITH PHOSPHOPANTETHEINE
'GYC':'CYS',                    ##  GYC CYS
'GYS':'GLY',                    ##  GYS GLY
'H5M':'PRO',                    ##  H5M PRO  TRANS-3-HYDROXY-5-METHYLPROLINE
'HHK':'ALA',                    ##  HHK ALA  (2S)-2,8-DIAMINOOCTANOIC ACID
'HIA':'HIS',                    ##  HIA HIS  L-HISTIDINE AMIDE
'HIC':'HIS',                    ##  HIC HIS  4-METHYL-HISTIDINE
'HIP':'HIS',                    ##  HIP HIS  ND1-PHOSPHONOHISTIDINE
'HIQ':'HIS',                    ##  HIQ HIS  MODIFIED HISTIDINE
'HIS':'HIS',                    ##  HIS HIS
'HLU':'LEU',                    ##  HLU LEU  BETA-HYDROXYLEUCINE
'HMF':'ALA',                    ##  HMF ALA  2-AMINO-4-PHENYL-BUTYRIC ACID
'HMR':'ARG',                    ##  HMR ARG  BETA-HOMOARGININE
'HPE':'PHE',                    ##  HPE PHE  HOMOPHENYLALANINE
'HPH':'PHE',                    ##  HPH PHE  PHENYLALANINOL GROUP
'HPQ':'PHE',                    ##  HPQ PHE  HOMOPHENYLALANINYLMETHANE
'HRG':'ARG',                    ##  HRG ARG  L-HOMOARGININE
'HSE':'SER',                    ##  HSE SER  L-HOMOSERINE
'HSL':'SER',                    ##  HSL SER  HOMOSERINE LACTONE
'HSO':'HIS',                    ##  HSO HIS  HISTIDINOL
'HTI':'CYS',                    ##  HTI CYS
'HTR':'TRP',                    ##  HTR TRP  BETA-HYDROXYTRYPTOPHANE
'HY3':'PRO',                    ##  HY3 PRO  3-HYDROXYPROLINE
'HYP':'PRO',                    ##  HYP PRO  4-HYDROXYPROLINE
'IAM':'ALA',                    ##  IAM ALA  4-[(ISOPROPYLAMINO)METHYL]PHENYLALANINE
'IAS':'ASP',                    ##  IAS ASP  ASPARTYL GROUP
'IGL':'ALA',                    ##  IGL ALA  ALPHA-AMINO-2-INDANACETIC ACID
'IIL':'ILE',                    ##  IIL ILE  ISO-ISOLEUCINE
'ILE':'ILE',                    ##  ILE ILE
'ILG':'GLU',                    ##  ILG GLU  GLU LINKED TO NEXT RESIDUE VIA CG
'ILX':'ILE',                    ##  ILX ILE  4,5-DIHYDROXYISOLEUCINE
'IML':'ILE',                    ##  IML ILE  N-METHYLATED
'IPG':'GLY',                    ##  IPG GLY  N-ISOPROPYL GLYCINE
'IT1':'LYS',                    ##  IT1 LYS
'IYR':'TYR',                    ##  IYR TYR  3-IODO-TYROSINE
'KCX':'LYS',                    ##  KCX LYS  CARBAMOYLATED LYSINE
'KGC':'LYS',                    ##  KGC LYS
'KOR':'CYS',                    ##  KOR CYS  MODIFIED CYSTEINE
'KST':'LYS',                    ##  KST LYS  N~6~-(5-CARBOXY-3-THIENYL)-L-LYSINE
'KYN':'ALA',                    ##  KYN ALA  KYNURENINE
'LA2':'LYS',                    ##  LA2 LYS
'LAL':'ALA',                    ##  LAL ALA  N,N-DIMETHYL-L-ALANINE
'LCK':'LYS',                    ##  LCK LYS
'LCX':'LYS',                    ##  LCX LYS  CARBAMYLATED LYSINE
'LDH':'LYS',                    ##  LDH LYS  N~6~-ETHYL-L-LYSINE
'LED':'LEU',                    ##  LED LEU  POST-TRANSLATIONAL MODIFICATION
'LEF':'LEU',                    ##  LEF LEU  2-5-FLUOROLEUCINE
'LET':'LYS',                    ##  LET LYS  ODIFIED LYSINE
'LEU':'LEU',                    ##  LEU LEU
'LLP':'LYS',                    ##  LLP LYS
'LLY':'LYS',                    ##  LLY LYS  NZ-(DICARBOXYMETHYL)LYSINE
'LME':'GLU',                    ##  LME GLU  (3R)-3-METHYL-L-GLUTAMIC ACID
'LNT':'LEU',                    ##  LNT LEU
'LPD':'PRO',                    ##  LPD PRO  L-PROLINAMIDE
'LSO':'LYS',                    ##  LSO LYS  MODIFIED LYSINE
'LYM':'LYS',                    ##  LYM LYS  DEOXY-METHYL-LYSINE
'LYN':'LYS',                    ##  LYN LYS  2,6-DIAMINO-HEXANOIC ACID AMIDE
'LYP':'LYS',                    ##  LYP LYS  N~6~-METHYL-N~6~-PROPYL-L-LYSINE
'LYR':'LYS',                    ##  LYR LYS  MODIFIED LYSINE
'LYS':'LYS',                    ##  LYS LYS
'LYX':'LYS',                    ##  LYX LYS  N''-(2-COENZYME A)-PROPANOYL-LYSINE
'LYZ':'LYS',                    ##  LYZ LYS  5-HYDROXYLYSINE
'M0H':'CYS',                    ##  M0H CYS  S-(HYDROXYMETHYL)-L-CYSTEINE
'M2L':'LYS',                    ##  M2L LYS
'M3L':'LYS',                    ##  M3L LYS  N-TRIMETHYLLYSINE
'MAA':'ALA',                    ##  MAA ALA  N-METHYLALANINE
'MAI':'ARG',                    ##  MAI ARG  DEOXO-METHYLARGININE
'MBQ':'TYR',                    ##  MBQ TYR
'MC1':'SER',                    ##  MC1 SER  METHICILLIN ACYL-SERINE
'MCL':'LYS',                    ##  MCL LYS  NZ-(1-CARBOXYETHYL)-LYSINE
'MCS':'CYS',                    ##  MCS CYS  MALONYLCYSTEINE
'MDO':'ALA',                    ##  MDO ALA
'MEA':'PHE',                    ##  MEA PHE  N-METHYLPHENYLALANINE
'MEG':'GLU',                    ##  MEG GLU  (2S,3R)-3-METHYL-GLUTAMIC ACID
'MEN':'ASN',                    ##  MEN ASN  GAMMA METHYL ASPARAGINE
'MET':'MET',                    ##  MET MET
'MEU':'GLY',                    ##  MEU GLY  O-METHYL-GLYCINE
'MFC':'ALA',                    ##  MFC ALA  CYCLIZED
'MGG':'ARG',                    ##  MGG ARG  MODIFIED D-ARGININE
'MGN':'GLN',                    ##  MGN GLN  2-METHYL-GLUTAMINE
'MHL':'LEU',                    ##  MHL LEU  N-METHYLATED, HYDROXY
'MHO':'MET',                    ##  MHO MET  POST-TRANSLATIONAL MODIFICATION
'MHS':'HIS',                    ##  MHS HIS  1-N-METHYLHISTIDINE
'MIS':'SER',                    ##  MIS SER  MODIFIED SERINE
'MLE':'LEU',                    ##  MLE LEU  N-METHYLATED
'MLL':'LEU',                    ##  MLL LEU  METHYL L-LEUCINATE
'MLY':'LYS',                    ##  MLY LYS  METHYLATED LYSINE
'MLZ':'LYS',                    ##  MLZ LYS  N-METHYL-LYSINE
'MME':'MET',                    ##  MME MET  N-METHYL METHIONINE
'MNL':'LEU',                    ##  MNL LEU  4,N-DIMETHYLNORLEUCINE
'MNV':'VAL',                    ##  MNV VAL  N-METHYL-C-AMINO VALINE
'MPQ':'GLY',                    ##  MPQ GLY  N-METHYL-ALPHA-PHENYL-GLYCINE
'MSA':'GLY',                    ##  MSA GLY  (2-S-METHYL) SARCOSINE
'MSE':'MET',                    ##  MSE MET  ELENOMETHIONINE
'MSO':'MET',                    ##  MSO MET  METHIONINE SULFOXIDE
'MTY':'PHE',                    ##  MTY PHE  3-HYDROXYPHENYLALANINE
'MVA':'VAL',                    ##  MVA VAL  N-METHYLATED
'N10':'SER',                    ##  N10 SER  O-[(HEXYLAMINO)CARBONYL]-L-SERINE
'NAL':'ALA',                    ##  NAL ALA  BETA-(2-NAPHTHYL)-ALANINE
'NAM':'ALA',                    ##  NAM ALA  NAM NAPTHYLAMINOALANINE
'NBQ':'TYR',                    ##  NBQ TYR
'NC1':'SER',                    ##  NC1 SER  NITROCEFIN ACYL-SERINE
'NCB':'ALA',                    ##  NCB ALA  CHEMICAL MODIFICATION
'NEP':'HIS',                    ##  NEP HIS  N1-PHOSPHONOHISTIDINE
'NFA':'PHE',                    ##  NFA PHE  MODIFIED PHENYLALANINE
'NIY':'TYR',                    ##  NIY TYR  META-NITRO-TYROSINE
'NLE':'LEU',                    ##  NLE LEU  NORLEUCINE
'NLN':'LEU',                    ##  NLN LEU  NORLEUCINE AMIDE
'NLO':'LEU',                    ##  NLO LEU  O-METHYL-L-NORLEUCINE
'NMC':'GLY',                    ##  NMC GLY  N-CYCLOPROPYLMETHYL GLYCINE
'NMM':'ARG',                    ##  NMM ARG  MODIFIED ARGININE
'NPH':'CYS',                    ##  NPH CYS
'NRQ':'ALA',                    ##  NRQ ALA
'NVA':'VAL',                    ##  NVA VAL  NORVALINE
'NYC':'ALA',                    ##  NYC ALA
'NYS':'CYS',                    ##  NYS CYS
'NZH':'HIS',                    ##  NZH HIS
'OAS':'SER',                    ##  OAS SER  O-ACETYLSERINE
'OBS':'LYS',                    ##  OBS LYS  MODIFIED LYSINE
'OCS':'CYS',                    ##  OCS CYS  CYSTEINE SULFONIC ACID
'OCY':'CYS',                    ##  OCY CYS  HYDROXYETHYLCYSTEINE
'OHI':'HIS',                    ##  OHI HIS  3-(2-OXO-2H-IMIDAZOL-4-YL)-L-ALANINE
'OHS':'ASP',                    ##  OHS ASP  O-(CARBOXYSULFANYL)-4-OXO-L-HOMOSERINE
'OLT':'THR',                    ##  OLT THR  O-METHYL-L-THREONINE
'OMT':'MET',                    ##  OMT MET  METHIONINE SULFONE
'OPR':'ARG',                    ##  OPR ARG  C-(3-OXOPROPYL)ARGININE
'ORN':'ALA',                    ##  ORN ALA  ORNITHINE
'ORQ':'ARG',                    ##  ORQ ARG  N~5~-ACETYL-L-ORNITHINE
'OSE':'SER',                    ##  OSE SER  O-SULFO-L-SERINE
'OTY':'TYR',                    ##  OTY TYR
'OXX':'ASP',                    ##  OXX ASP  OXALYL-ASPARTYL ANHYDRIDE
'P1L':'CYS',                    ##  P1L CYS  S-PALMITOYL CYSTEINE
'P2Y':'PRO',                    ##  P2Y PRO  (2S)-PYRROLIDIN-2-YLMETHYLAMINE
'PAQ':'TYR',                    ##  PAQ TYR  SEE REMARK 999
'PAT':'TRP',                    ##  PAT TRP  ALPHA-PHOSPHONO-TRYPTOPHAN
'PBB':'CYS',                    ##  PBB CYS  S-(4-BROMOBENZYL)CYSTEINE
'PBF':'PHE',                    ##  PBF PHE  PARA-(BENZOYL)-PHENYLALANINE
'PCA':'PRO',                    ##  PCA PRO  5-OXOPROLINE
'PCS':'PHE',                    ##  PCS PHE  PHENYLALANYLMETHYLCHLORIDE
'PEC':'CYS',                    ##  PEC CYS  S,S-PENTYLTHIOCYSTEINE
'PF5':'PHE',                    ##  PF5 PHE  2,3,4,5,6-PENTAFLUORO-L-PHENYLALANINE
'PFF':'PHE',                    ##  PFF PHE  4-FLUORO-L-PHENYLALANINE
'PG1':'SER',                    ##  PG1 SER  BENZYLPENICILLOYL-ACYLATED SERINE
'PG9':'GLY',                    ##  PG9 GLY  D-PHENYLGLYCINE
'PHA':'PHE',                    ##  PHA PHE  PHENYLALANINAL
'PHD':'ASP',                    ##  PHD ASP  2-AMINO-4-OXO-4-PHOSPHONOOXY-BUTYRIC ACID
'PHE':'PHE',                    ##  PHE PHE
'PHI':'PHE',                    ##  PHI PHE  IODO-PHENYLALANINE
'PHL':'PHE',                    ##  PHL PHE  L-PHENYLALANINOL
'PHM':'PHE',                    ##  PHM PHE  PHENYLALANYLMETHANE
'PIA':'ALA',                    ##  PIA ALA  FUSION OF ALA 65, TYR 66, GLY 67
'PLE':'LEU',                    ##  PLE LEU  LEUCINE PHOSPHINIC ACID
'PM3':'PHE',                    ##  PM3 PHE
'POM':'PRO',                    ##  POM PRO  CIS-5-METHYL-4-OXOPROLINE
'PPH':'LEU',                    ##  PPH LEU  PHENYLALANINE PHOSPHINIC ACID
'PPN':'PHE',                    ##  PPN PHE  THE LIGAND IS A PARA-NITRO-PHENYLALANINE
'PR3':'CYS',                    ##  PR3 CYS  INE DTT-CYSTEINE
'PRO':'PRO',                    ##  PRO PRO
'PRQ':'PHE',                    ##  PRQ PHE  PHENYLALANINE
'PRR':'ALA',                    ##  PRR ALA  3-(METHYL-PYRIDINIUM)ALANINE
'PRS':'PRO',                    ##  PRS PRO  THIOPROLINE
'PSA':'PHE',                    ##  PSA PHE
'PSH':'HIS',                    ##  PSH HIS  1-THIOPHOSPHONO-L-HISTIDINE
'PTH':'TYR',                    ##  PTH TYR  METHYLENE-HYDROXY-PHOSPHOTYROSINE
'PTM':'TYR',                    ##  PTM TYR  ALPHA-METHYL-O-PHOSPHOTYROSINE
'PTR':'TYR',                    ##  PTR TYR  O-PHOSPHOTYROSINE
'PYA':'ALA',                    ##  PYA ALA  3-(1,10-PHENANTHROL-2-YL)-L-ALANINE
'PYC':'ALA',                    ##  PYC ALA  PYRROLE-2-CARBOXYLATE
'PYR':'SER',                    ##  PYR SER  CHEMICALLY MODIFIED
'PYT':'ALA',                    ##  PYT ALA  MODIFIED ALANINE
'PYX':'CYS',                    ##  PYX CYS  S-[S-THIOPYRIDOXAMINYL]CYSTEINE
'R1A':'CYS',                    ##  R1A CYS
'R1B':'CYS',                    ##  R1B CYS
'R1F':'CYS',                    ##  R1F CYS
'R7A':'CYS',                    ##  R7A CYS
'RC7':'ALA',                    ##  RC7 ALA
'RCY':'CYS',                    ##  RCY CYS
'S1H':'SER',                    ##  S1H SER  1-HEXADECANOSULFONYL-O-L-SERINE
'SAC':'SER',                    ##  SAC SER  N-ACETYL-SERINE
'SAH':'CYS',                    ##  SAH CYS  S-ADENOSYL-L-HOMOCYSTEINE
'SAR':'GLY',                    ##  SAR GLY  SARCOSINE
'SBD':'SER',                    ##  SBD SER
'SBG':'SER',                    ##  SBG SER  MODIFIED SERINE
'SBL':'SER',                    ##  SBL SER
'SC2':'CYS',                    ##  SC2 CYS  N-ACETYL-L-CYSTEINE
'SCH':'CYS',                    ##  SCH CYS  S-METHYL THIOCYSTEINE GROUP
'SCS':'CYS',                    ##  SCS CYS  MODIFIED CYSTEINE
'SCY':'CYS',                    ##  SCY CYS  CETYLATED CYSTEINE
'SDP':'SER',                    ##  SDP SER
'SEB':'SER',                    ##  SEB SER  O-BENZYLSULFONYL-SERINE
'SEC':'ALA',                    ##  SEC ALA  2-AMINO-3-SELENINO-PROPIONIC ACID
'SEL':'SER',                    ##  SEL SER  2-AMINO-1,3-PROPANEDIOL
'SEP':'SER',                    ##  SEP SER  E PHOSPHOSERINE
'SER':'SER',                    ##  SER SER
'SET':'SER',                    ##  SET SER  AMINOSERINE
'SGB':'SER',                    ##  SGB SER  MODIFIED SERINE
'SGR':'SER',                    ##  SGR SER  MODIFIED SERINE
'SHC':'CYS',                    ##  SHC CYS  S-HEXYLCYSTEINE
'SHP':'GLY',                    ##  SHP GLY  (4-HYDROXYMALTOSEPHENYL)GLYCINE
'SIC':'ALA',                    ##  SIC ALA
'SLZ':'LYS',                    ##  SLZ LYS  L-THIALYSINE
'SMC':'CYS',                    ##  SMC CYS  POST-TRANSLATIONAL MODIFICATION
'SME':'MET',                    ##  SME MET  METHIONINE SULFOXIDE
'SMF':'PHE',                    ##  SMF PHE  4-SULFOMETHYL-L-PHENYLALANINE
'SNC':'CYS',                    ##  SNC CYS  S-NITROSO CYSTEINE
'SNN':'ASP',                    ##  SNN ASP  POST-TRANSLATIONAL MODIFICATION
'SOC':'CYS',                    ##  SOC CYS  DIOXYSELENOCYSTEINE
'SOY':'SER',                    ##  SOY SER  OXACILLOYL-ACYLATED SERINE
'SUI':'ALA',                    ##  SUI ALA
'SUN':'SER',                    ##  SUN SER  TABUN CONJUGATED SERINE
'SVA':'SER',                    ##  SVA SER  SERINE VANADATE
'SVV':'SER',                    ##  SVV SER  MODIFIED SERINE
'SVX':'SER',                    ##  SVX SER  MODIFIED SERINE
'SVY':'SER',                    ##  SVY SER  MODIFIED SERINE
'SVZ':'SER',                    ##  SVZ SER  MODIFIED SERINE
'SXE':'SER',                    ##  SXE SER  MODIFIED SERINE
'TBG':'GLY',                    ##  TBG GLY  T-BUTYL GLYCINE
'TBM':'THR',                    ##  TBM THR
'TCQ':'TYR',                    ##  TCQ TYR  MODIFIED TYROSINE
'TEE':'CYS',                    ##  TEE CYS  POST-TRANSLATIONAL MODIFICATION
'TH5':'THR',                    ##  TH5 THR  O-ACETYL-L-THREONINE
'THC':'THR',                    ##  THC THR  N-METHYLCARBONYLTHREONINE
'THR':'THR',                    ##  THR THR
'TIH':'ALA',                    ##  TIH ALA  BETA(2-THIENYL)ALANINE
'TMD':'THR',                    ##  TMD THR  N-METHYLATED, EPSILON C ALKYLATED
'TNB':'CYS',                    ##  TNB CYS  S-(2,3,6-TRINITROPHENYL)CYSTEINE
'TOX':'TRP',                    ##  TOX TRP
'TPL':'TRP',                    ##  TPL TRP  TRYTOPHANOL
'TPO':'THR',                    ##  TPO THR  HOSPHOTHREONINE
'TPQ':'ALA',                    ##  TPQ ALA  2,4,5-TRIHYDROXYPHENYLALANINE
'TQQ':'TRP',                    ##  TQQ TRP
'TRF':'TRP',                    ##  TRF TRP  N1-FORMYL-TRYPTOPHAN
'TRN':'TRP',                    ##  TRN TRP  AZA-TRYPTOPHAN
'TRO':'TRP',                    ##  TRO TRP  2-HYDROXY-TRYPTOPHAN
'TRP':'TRP',                    ##  TRP TRP
'TRQ':'TRP',                    ##  TRQ TRP
'TRW':'TRP',                    ##  TRW TRP
'TRX':'TRP',                    ##  TRX TRP  6-HYDROXYTRYPTOPHAN
'TTQ':'TRP',                    ##  TTQ TRP  6-AMINO-7-HYDROXY-L-TRYPTOPHAN
'TTS':'TYR',                    ##  TTS TYR
'TY2':'TYR',                    ##  TY2 TYR  3-AMINO-L-TYROSINE
'TY3':'TYR',                    ##  TY3 TYR  3-HYDROXY-L-TYROSINE
'TYB':'TYR',                    ##  TYB TYR  TYROSINAL
'TYC':'TYR',                    ##  TYC TYR  L-TYROSINAMIDE
'TYI':'TYR',                    ##  TYI TYR  3,5-DIIODOTYROSINE
'TYN':'TYR',                    ##  TYN TYR  ADDUCT AT HYDROXY GROUP
'TYO':'TYR',                    ##  TYO TYR
'TYQ':'TYR',                    ##  TYQ TYR  AMINOQUINOL FORM OF TOPA QUINONONE
'TYR':'TYR',                    ##  TYR TYR
'TYS':'TYR',                    ##  TYS TYR  INE SULPHONATED TYROSINE
'TYT':'TYR',                    ##  TYT TYR
'TYX':'CYS',                    ##  TYX CYS  S-(2-ANILINO-2-OXOETHYL)-L-CYSTEINE
'TYY':'TYR',                    ##  TYY TYR  IMINOQUINONE FORM OF TOPA QUINONONE
'TYZ':'ARG',                    ##  TYZ ARG  PARA ACETAMIDO BENZOIC ACID
'UMA':'ALA',                    ##  UMA ALA
'VAD':'VAL',                    ##  VAD VAL  DEAMINOHYDROXYVALINE
'VAF':'VAL',                    ##  VAF VAL  METHYLVALINE
'VAL':'VAL',                    ##  VAL VAL
'VDL':'VAL',                    ##  VDL VAL  (2R,3R)-2,3-DIAMINOBUTANOIC ACID
'VLL':'VAL',                    ##  VLL VAL  (2S)-2,3-DIAMINOBUTANOIC ACID
'VME':'VAL',                    ##  VME VAL  O- METHYLVALINE
'X9Q':'ALA',                    ##  X9Q ALA
'XX1':'LYS',                    ##  XX1 LYS  N~6~-7H-PURIN-6-YL-L-LYSINE
'XXY':'ALA',                    ##  XXY ALA
'XYG':'ALA',                    ##  XYG ALA
'YCM':'CYS',                    ##  YCM CYS  S-(2-AMINO-2-OXOETHYL)-L-CYSTEINE
'YOF':'TYR' }                   ##  YOF TYR  3-FLUOROTYROSINE

AA = {"R":"ARG","ARG":"R","H":"HIS","HIS":"H","K":"LYS","LYS":"K",
      "D":"ASP","ASP":"D","E":"GLU","GLU":"E","S":"SER","SER":"S",
      "T":"THR","THR":"T","N":"ASN","ASN":"N","Q":"GLN","GLN":"Q",
      "C":"CYS","CYS":"C","U":"SEC","SEC":"U","G":"GLY","GLY":"G",
      "P":"PRO","PRO":"P","A":"ALA","ALA":"A","V":"VAL","VAL":"V",
      "I":"ILE","ILE":"I","L":"LEU","LEU":"L","M":"MET","MET":"M",
      "F":"PHE","PHE":"F","Y":"TYR","TYR":"Y","W":"TRP","TRP":"W"
      }
      
LONGER_NAMES={'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
              'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
              'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
              'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
              'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V','TYS':'Y'}

#Bad ligands are likely solvent non-selective ligands not important
BAD_LIGANDS = ['HOH','NAG','ACT','PEG','EDO','ACE']
#Bad Rosetta ligands are ones that are problematic for Rosetta for non-obvious reasons
#For T19,T42,T16,DP7,DI2,DI3,DI4,DI5,S2C,368,HDB,412,427,ABH,34B,39B,39E,A48,C08,BON,5AB,4U7,3T4 Rosetta doesn't support atom type B. This may be a flag and if found, they should be fine
#UNL is unknown ligand and will never work with Rosetta
#UNX is unknown atom or ion
BAD_ROSETTA_LIGANDS = ['0TU','SCN','PHG','AZI','CMO','URE','T19','T42','T16','AUC','DP7','DI2','DI3','DI4','DI5','HG2','S2C','GOL','368','HDB','412','427','ABH','34B','39B','39E','CCN','CO2','MNC','A48','C08','BE7','UNL','BON','5AB','4U7','3T4','UNX']
DNA = ['DT','DG','DC','DA']

def check_ligand(ligand,resname):
    '''Handle strange ligand edge-cases'''
    passed = False
    #Have seen cases where sulfure atom is missing (1F5K)
    if resname=='SO4':
        for atom in ligand:
            atomname = atom[12:16]
            if atomname.strip()=='S':
                passed = True
    #Have seen cases where HGB is missing heavy atoms (1F2W)
    elif resname=='HGB':
        if len(ligand)<11:
            passed = False
        else:
            passed = True
    #Example of ACD missing heavy atoms (1GNJ)
    elif resname=='ACD':
        if len(ligand)<22:
            passed = False
        else:
            passed = True
    #Example of CPS missing heavy atoms (1XU7)
    elif resname=='CPS':
        if len(ligand)<42:
            passed = False
        else:
            passed = True
    #Example of COA missing heavy atoms (2JFK)
    elif resname=='COA':
        if len(ligand)<48:
            passed = False
        else:
            passed = True
    elif resname=='IDS':
        if len(ligand)<17:
            passed = False
        else:
            passed = True
    #Example (3V7X)
    elif resname=='MBO':
        if len(ligand)<10:
            passed = False
        else:
            passed = True
    #Example (4DO6)
    elif resname=='2PE':
        if len(ligand)<28:
            passed = False
        else:
            passed = True
    #These cause issues with Rosetta and it's not immediately clear why:
    elif resname in BAD_ROSETTA_LIGANDS:
        passed = False
    else:
        passed = True
    if not passed:
        print("Warning: skipping problematic ligand: %s" % resname)    
    
    return passed
    
def renumber_chain(pdb,chain):
    with open(pdb) as infile:
        pdblines = infile.readlines()
    to_pdb = ['NA']
    isopdb = ''
    pdbtext = os.path.basename(os.path.splitext(pdb)[0])
    isopdbname = "%s_%s.pdb" % (pdbtext,chain)
    atomnum = 1
    resnum = 1
    for line in pdblines:
        if len(line)>20 and line[:4] =='ATOM' and line[21]==chain:
            try:
                resnum = int(line[22:26])
                aa = AA[line[17:20]]
            except ValueError:
                sys.exit("Error parsing %s" % pdb)
            new_atomnum = "%5d" % atomnum
            if to_pdb[-1] != resnum:
                to_pdb.append(resnum)
            new_atomnum = "%5d" % atomnum
            new_resnum = len(to_pdb)-1
            new_line = line[:6]+new_atomnum+line[11:22]+"%4d" % new_resnum + line[26:]
            isopdb += new_line
            atomnum += 1
    with open(isopdbname, 'w') as outfile:
        outfile.write(isopdb)
        outfile.write("TER\n")
    return isopdbname

def query_modres(aa):
    global MODRES
    try:
        caa = MODRES[aa]
        return caa
    except KeyError:
        return None
        
def check_mutation(muts,pdb,chain):
    seenmuts = set()
    seen_chain = False
    runner = ['NA']
    temp = []
    for line in pdb:
        if len(line)>20 and line[:4] =='ATOM' and line[21]==chain:
            try:
                resnum = int(line[22:26])
                aa = AA[line[17:20]]
            except ValueError:
                sys.exit("Fatal: Error parsing %s" % rosetta_pdb)
            seen_chain = True
            if runner[-1] != resnum:
                runner.append(resnum)
                for x in range(len(muts)):
                    if resnum == int(muts[x][1]):
                        if resnum in seenmuts:
                            print("Warning, resnum %d listed more than once for mutations, only keeping first instance" % resnum)
                            continue
                        if aa == muts[x][0]:
                            temp.append(muts[x])
                            seenmuts.add(resnum)
                        else:
                            print("Warning, AA for mutation %s does not match pdb: pos %d is %s, skipping this mutation" % ("".join(str(y) for y in muts[x]), resnum,aa))
                            seenmuts.add(resnum)
    if not seen_chain:
        return None
    else:
        for x in muts:
            if x[1] not in seenmuts:
                print("Warning: Resnum %d does not appear in pdb. Dropping mutations at this position!" % x[1])        
        return list(temp)
        
def parse_mutation(mutation,pos_only=False):
    if mutation is None:
        return None
    if pos_only:
        try:
            pos = int(mutation)
            return pos
        except:
            pass
    dix = re.search("\d", mutation)
    try:
        assert dix
    except AssertionError:
        return None
    startaa = mutation[:dix.start()].upper()
    if len(startaa) == 3:
        startaa = AA[startaa]
    mutation = mutation[dix.start():]
    if not pos_only:
        try:
            assert startaa.isalpha() and (len(startaa)==1 or len(startaa)==3)
        except AssertionError:
            return None
    six = re.search("[a-zA-Z]",mutation)
    try:
        assert six
        if pos_only:
            return int(mutation[:six.start()])
    except AssertionError:
        if pos_only:
            return int(mutation)
        else:
            return None
    resnum = int(mutation[:six.start()])
    endaa = mutation[six.start():].upper()
    if len(endaa) == 3:
        endaa = AA[endaa]
    try:
        assert endaa.isalpha() and (len(endaa)==1 or len(endaa)==3)
    except AssertionError:
        return None
    mutation = [startaa,resnum,endaa]
    return mutation        
    
class cleaner(object):
    def __init__(self,inpdb,verbose=True,human=True,outpath = None,nmr=False,hardhuman=True):
        self._verbose = verbose
        self._human = human
        self._inpdbfile = inpdb
        self._ignorechain = 1
        self.outfile = ""
        self.outfile_rev = ""
        self._pdbfile = ""
        self._pdbfile_rev = ""
        self._gz = False
        self.nmr = nmr
        self._keepdna = False
        self.hardhuman = hardhuman #Keep only human residues no matter what. If false, keep nonhuman residues that make up their own chain
        if outpath is not None:
            self.outfile = "%s/" % outpath.rstrip('/')
            self.outfile_rev = "%s/" % outpath.rstrip('/')
        if inpdb[-2:] == 'gz':
            self._gz = True
            p1 = os.path.basename(os.path.splitext(inpdb[:-3])[0])
        else:
            p1 = os.path.basename(os.path.splitext(inpdb)[0])        

        if not self._gz and inpdb[-1].isdigit():
            self._biomodel=inpdb[-1]
            self.outfile += "%s_rosetta.pdb%s" % (p1, self._biomodel)
            self.outfile_rev += "%s_rev.pdb%s" % (p1, self._biomodel)
        elif self._gz and inpdb[-4].isdigit():
            self._biomodel=inpdb[-4]
            self.outfile += "%s_rosetta.pdb%s" % (p1, self._biomodel)
            self.outfile_rev += "%s_rev.pdb%s" % (p1, self._biomodel)
        else:
            self._biomodel= "0"
            self.outfile += "%s_rosetta.pdb" % p1
            self.outfile_rev += "%s_rev.pdb" % p1
            
        try:
            assert os.path.isfile(self._inpdbfile)
        except AssertionError:
            sys.exit("Error: %s not found" % self._inpdbfile)
        if self._gz:
            with gzip.open(self._inpdbfile,'rt') as infile:
                self._lines = infile.readlines()
        else:
            with open(self._inpdbfile,'rt') as infile:
                self._lines = infile.readlines()
        try:
            assert len(self._lines)>0
        except AssertionError:
            sys.exit("Error unable to parse %s" % self._inpdbfile)

    def check_and_print_pdb(self, residue_buffer, bbcheck):
        hasCA = False
        hasN = False
        hasC = False
        resnum = ""
        restype = ""
        for line in residue_buffer:
            resnum = int(line[22:26])
            atomname = line[12:16]
            occupancy = float(line[55:60])
            restype = line[17:20]
            if atomname == " CA ":# and occupancy > 0.0:
                hasCA = True
            if atomname == " N  ":# and occupancy > 0.0:
                hasN = True
            if atomname == " C  ":# and occupancy > 0.0:
                hasC = True
            if line[:6]=='HETATM' or (self._keepdna and restype in [' DT',' DG',' DC',' DA','DT','DG','DC','DA']):
                hasCA = True
                hasN = True
                hasC = True
        if line[:6] == 'HETATM' and bbcheck:
            hasCA = check_ligand(residue_buffer,restype)
        # If cleaning for use in Rosetta, each residue must have 3 bb atoms present
        # or Rosetta will skip the residue.
        if (bbcheck and hasCA and hasN and hasC) or not bbcheck:
            for line in residue_buffer:
                if line[21:26] in self._seen_residues and not self.nmr:
                    print("Warning: residue %s seen more than once!" % line[21:26])
                aa = line[17:20]
                newnum = line[22:27]
                chain = line[21]
                if chain == ' ':
                    if self._verbose:
                        print("Reverting blank chain to A:")
                        print(line)
                    chain = 'A'
                line_edit = line[0:17] + aa + line[20:21] + chain + newnum + line[27:]
                # write the residue line to rev
                if not bbcheck:
                    self._pdbfile_rev += line_edit
                else:
                    self._pdbfile += line_edit            
            self._seen_residues.append(line[21:26])
            return True
        else:
            if self._verbose:
                print("skipped res with missing BB atoms (rosetta file only):")
                print(line)
            return False

    def write(self, bbcheck=False,chain=None,keepligands=False,keepdna=False):
        oldresnum = '    '
        current_atoms = []
        residue_buffer = []
        residue_invalid = False 
        modifiedres = ''       
        self._seen_residues = []
        self._non_human_residues = {}
        self._human_residues = {}
        self._seen_residues = []
        self._skipped = 0
        self._skipped_alts = 0
        self._skipped_nonhuman = 0
        self._keepdna = keepdna
        if bbcheck:
            self._pdbfile = ""
        else:
            self._pdbfile_rev = ""

        for i in range(len(self._lines)):
            line = self._lines[i]
            if line.startswith("DBREF") and self._ignorechain:
                if len(line.split())<8: continue
                if not line.split()[7].endswith("HUMAN"):
                    if line.split()[2] not in self._non_human_residues:
                        try:
                            startpos = int(line.split()[3])
                        #In case DBREF range starts or ends on alternative position
                        except ValueError:
                            try:
                                startpos = int(line.split()[3][:-1])
                            except ValueError:
                                LOGGER.info("Warning: Unable to parse non-human DBREF range, remove manually!")
                                continue
                        try:
                            endpos = int(line.split()[4])
                        except ValueError:
                            try:
                                endpos = int(line.split()[4][:-1])
                            except:
                                LOGGER.info("Warning: Unable to parse non-human DBREF range, remove manually!")
                                continue
                        self._non_human_residues[line.split()[2]] = list(range(startpos,endpos+1))
                    else:
                        self._non_human_residues[line.split()[2]] += list(range(startpos,endpos+1))           
                elif line.split()[2] not in self._human_residues:
                    try:
                        startpos = int(line.split()[3])
                    except ValueError:
                        try:
                            startpos = int(line.split()[3][:-1])
                        except ValueError:
                            LOGGER.info("Warning! Unable to parse human DBREF range!! clean manually!")
                            continue
                    try:
                        endpos = int(line.split()[4])
                    except ValueError:
                        try:
                            endpos = int(line.split()[4][:-1])
                        except:
                            LOGGER.info("Warning! Unable to parse human DBREF range!! clean manually!")
                            continue
                    self._human_residues[line.split()[2]] = list(range(startpos,endpos+1))
                else:
                    try:
                        startpos = int(line.split()[3])
                    except ValueError:
                        try:
                            startpos = int(line.split()[3][:-1])
                        except ValueError:
                            LOGGER.info("Warning! Unable to parse human DBREF range!! clean manually!")
                            continue
                    try:
                        endpos = int(line.split()[4])
                    except ValueError:
                        try:
                            endpos = int(line.split()[4][:-1])
                        except:
                            LOGGER.info("Warning! Unable to parse human DBREF range!! clean manually!")
                            continue
                    self._human_residues[line.split()[2]] += list(range(startpos,endpos+1))
            if len(line) > 5 and line[:6] == 'ENDMDL': break  # Its an NMR model.
            if line[0:3] == 'TER':
                continue
            if self._ignorechain:
                line_edit = line
                #Fix modified residues which are stored as HETATMs
                if line[0:6] == 'HETATM' and line[17:20] in BAD_LIGANDS:
                    ok = False
                    if keepligands and line[17:20] not in BAD_LIGANDS:
                        ok = True
                    # Is it a modified residue ?
                    if line[17:20] in MODRES:
                      # if so replace it with its canonical equivalent !
                        line_edit = 'ATOM  '+line[6:17]+MODRES[line[17:20]] + line[20:]
                        modifiedres = modifiedres + line[17:20] + ',  '
                        ok = True
                    # other substitution (of atoms mainly)
                    if (line[17:20] == 'MSE'):  # Selenomethionine
                        if (line_edit[12:14] == 'SE'):
                            line_edit = line_edit[0:12]+' S'+line_edit[14:]
                        if len(line_edit) > 75:
                            if (line_edit[76:78] == 'SE'):
                                line_edit = line_edit[0:76]+' S'+line_edit[78:]
                    if not ok:
                        continue  # skip this atom if we havnt found a conversion
        
                if line_edit[0:4] == 'ATOM' or (line_edit[0:6]=='HETATM' and keepligands):
                    curatom = line_edit[12:16].strip()            
                    altcode = line[16]
                    resnum = line_edit[22:27]
                    tag = line_edit[21:27]
                    curchain = line[21]
                    if chain is not None and chain!=curchain: continue
                    if tag != oldresnum:
                        if residue_buffer != []:  # is there a residue in the buffer ?
                            if not residue_invalid:
                                if not self.check_and_print_pdb(residue_buffer,bbcheck):
                                    self._skipped += 1
                        current_atoms = []
                        residue_buffer = []
                        residue_invalid = False
        
                    insres = line[26]
                    if insres != ' ' and not residue_invalid:
                        if self._verbose and not residue_invalid:
                            LOGGER.info("Skipping residue with insertion code:")
                            LOGGER.info(line)
                        residue_invalid = True
                        self._skipped += 1
                        
                    longname = line_edit[17:20]
                    if (longname not in LONGER_NAMES and not residue_invalid) and not (line_edit[:6]=='HETATM' and keepligands) and not (keepdna and longname.strip() in DNA):
                        if self._verbose:
                            LOGGER.info("Unable to reconcile residue name:")
                            LOGGER.info(line)
                        self._skipped += 1
                        residue_invalid = True
                    if self._human and (len(self._non_human_residues) >0 or len(self._human_residues) > 0) and ((line[21] in self._non_human_residues and int(line[22:26]) in self._non_human_residues[line[21]]) or line[21] not in self._human_residues or int(line[22:26]) not in self._human_residues[line[21]]) and not residue_invalid:
                        if self._verbose:
                            LOGGER.info("skipping non-human residue:")
                            LOGGER.info(line)
                        self._skipped_nonhuman += 1
                        residue_invalid = True
                    #Only take novel atom coords, keep only first alt pos listed.
                    if not residue_invalid:
                        if curatom not in current_atoms:
                            current_atoms.append(curatom)
                        else:
                            if self._verbose:
                                LOGGER.info("Skipping alternate atom coordinate:")
                                LOGGER.info(line)
                            self._skipped_alts += 1
                            continue
        
                    oldresnum = tag
        
                    # Get rid of alternate position designation
                    if line_edit[16:17] != ' ':
                        line_edit = line_edit[:16]+' '+line_edit[17:]
                    residue_buffer.append(line_edit)

        if not residue_invalid and residue_buffer !=[]:
            if not self.check_and_print_pdb(residue_buffer,bbcheck):
                self._skipped += 1
        
        if self._skipped > 0:
            LOGGER.info("%d bad residues have been skipped" % self._skipped)
        if self._skipped_nonhuman > 0:
            LOGGER.info("%d non-human residues have been skipped" % self._skipped_nonhuman)
        if self._skipped_alts > 0:
            LOGGER.info("%d alternative atom positions skipped" % self._skipped_alts)
        if bbcheck:
            if len(self._pdbfile)==0 or len([x for x in self._pdbfile.split('\n') if x.startswith('ATOM')])==0: raise Exception("No residues remain in pdb structure during processing")

            with open(self.outfile, 'w') as wout:
                wout.write(self._pdbfile)
                wout.write("TER\n")
                return self.outfile
        else:
            if len(self._pdbfile_rev)==0 or len([x for x in self._pdbfile_rev.split('\n') if x.startswith('ATOM')])==0: raise Exception("No residues remain in pdb structure during processing")

            with open(self.outfile_rev,'w') as wout:
                wout.write(self._pdbfile_rev)
                wout.write("TER\n")
                return self.outfile_rev
