#!/usr/bin/env python
# Chris Moth 2017-Nov 30
# Extract ONLY human entries from the ~150GB match_complete.xml
# The resulting match_human.xml file in the same directory (specify in -c configfile)
# is used to help generate domain graphics for isoform-specific uniprot ACs in the pipeline
# output reports

import argparse,configparser
import json
import traceback
import sys,os,csv,time,pdb,glob,gzip,shutil
import subprocess as sp
from multiprocessing import cpu_count
from collections import defaultdict

import logging
LOGGER = logging.getLogger(__name__)

import xml.etree.ElementTree as ET
import sys

# # Setup the Config File Parser
# conf_parser = argparse.ArgumentParser(add_help=False)
# conf_parser.add_argument("-c", "--conf_file",
# help="Specify config file", metavar="FILE", required=True)
# args, remaining_argv = conf_parser.parse_known_args()
# conf_file = args.conf_file
# print 'Reading config file: %s'%conf_file
# config = ConfigParser.SafeConfigParser()
# config.read([conf_file])
# defaults = {}
# defaults.update(dict(config.items("Genome_PDB_Mapper")))
# print 'Interpro directory: %s'%defaults['interpro_dir']
# 
# xml_human_file = defaults['interpro_dir'] + 'match_humanonly.xml'

class Cxml2json(object):
  __tree = None
  __root = None
  def __init__(self,xml_human_filename):
    # Init these structures one time only - and only if the class is instantiated once
    if not Cxml2json.__tree:
      try:
          LOGGER.info("Attempting to parse %s" % xml_human_filename)
          Cxml2json.__tree = ET.parse(xml_human_filename)
          Cxml2json.__root = Cxml2json.__tree.getroot()
      except Exception as ex:
          LOGGER.error("xml parse error of %s:\n%s" % (xml_human_filename, str(ex)))
      

  def xmlForBaseUnp(self,uniprot_ac):
    found_list = []
    base_unp = uniprot_ac.split('-')[0]
    for child in Cxml2json.__root:
      if child.tag == 'protein' and child.get('id').startswith(base_unp):
        print((child.get('id')))
        found_list.append(child)
    return found_list

  def isCanonical(self,uniprot_ac):
    # If init failed, then don't proceed
    if Cxml2json.__root is None:
        LOGGER.warn("For %s, Cxml2json__root not set by __init__" % uniprot_ac)
        return False
 
    unp_split = uniprot_ac.split('-')
    if len(unp_split) < 2:
      return True  # no dash means you gave a canonical unp
    # DO NOT CHECK IN THIS TERRIBLE PATCH BELOW
    # if uniprot_ac == 'Q8NDA2-5':
    #     return True
    #  DO NOT CHECK IN THIS TERRIBLE PATCH ^^^ ABOVE
    # When given ABC-n, you want to make sure the sequence CRC64 matches
    # between base-unp and specific isoform, and then you can say YES!
    base_unp = unp_split[0]
    isoform_node = Cxml2json.__root.findall("protein/[@id='%s']"%uniprot_ac)
    if len(isoform_node) != 1:
      failure_message = "%s was not found in xml file.  Return 'not canonical'"%uniprot_ac
      LOGGER.warn(failure_message)
      return False
    baseunp_node = Cxml2json.__root.findall("protein/[@id='%s']"%base_unp)
    if len(baseunp_node) != 1:
      failure_message = "%s was not found in xml file"%base_unp
      LOGGER.warn(failure_message)
      return False
    return isoform_node[0].get('crc64') == baseunp_node[0].get('crc64')

  def PFAMgraphicsJSONfromUnp(self,uniprot_ac) -> dict:
    if Cxml2json.__root is None:
        LOGGER.warn("For %s, Cxml2json__root not set by __init__" % uniprot_ac)
        return {} 

    if self.isCanonical(uniprot_ac):
      uniprot_ac = uniprot_ac.split('-')[0]
    ### WARNING HACK BELOW
    if uniprot_ac == 'Q9HCJ0':
       return {}
    ### WARNING HACK ABOVE
    # The final JSON graphics string has length, and lists of dictionaries called "regions" and "motifs"
    # See http://pfam.xfam.org/help#tabview=tab10 "A Guide to the PFam Domain Graohics"
    # which explains how to build up the JSON in a little detail
    # http://pfam.xfam.org/help/domain_graphics_example.html  is an invaluable resource to "play" with different
    # regions and mitofs.
    PFAMxml = Cxml2json.__root.findall("protein/[@id='%s']"%uniprot_ac)
    assert PFAMxml != None,"Could not find %s in interpro XML"%uniprot_ac
    # Sometimes we get an empty list of graphics.  In these cases, we can do no more
    if not PFAMxml:
        return {}
    color_disorder = "#cccccc"
    color_coiledcoil = "#32cd32"

    # Rotate through a few colors as familys/regions are identified on the grey line for the entire protein
    graphicsJSON={}
    graphicsJSON['length'] = PFAMxml[0].get('length')
    JSONregions = []
    JSONmotifs = []
    for child in PFAMxml:
      for match in child:
        if match.tag == 'match':
          matchdict = defaultdict(lambda: '')
          matchdict['id'] = match.get('id')
          matchdict['name'] = match.get('name')
          matchdict['dbname'] = match.get('dbname')
          matchdict['status'] = match.get('statis')
          matchdict['evd'] = match.get('evd')
          for lcn_or_ipr in match:
            if lcn_or_ipr.tag == 'lcn':
              # print "lcn_or_ipr=",ET.tostring(lcn_or_ipr)
              matchdict['start'] = int(lcn_or_ipr.get('start'))
              matchdict['end'] = int(lcn_or_ipr.get('end'))
              matchdict['score'] = lcn_or_ipr.get('score')
            elif lcn_or_ipr.tag == 'ipr':
              matchdict['description'] = lcn_or_ipr.get('name')
              # print "lcn_or_ipr=",ET.tostring(lcn_or_ipr)
    
          # print matchdict
    
          if matchdict['dbname'] == 'PFAM':
            JSONregions.append(
              {'start': matchdict['start'], 'end': matchdict['end'], 
               # 'alistart': matchdict['start'], 'aliend': matchdict['end'], 
               'startStyle': 'straight', 'endStyle': 'straight', 
               'type': "pfama", 
               'text': matchdict['name'],
               'href': "/family/%s"%matchdict['id'],
               # Assign colours at end, after sorting 'colour': color_domains[color_index % len(color_domains)],
               'metadata': {'type': 'Domain' if 'omain' in matchdict['description'] else 'Family', 'start': matchdict['start'], 'end': matchdict['end'], 
                            'database': 'pfam',
                            'description': matchdict['description'],
                             # 'aliStart': matchdict['start'], 'aliEnd': matchdict['end'], 
                             'score': matchdict['score']}
              })
            # color_index += 1
          elif matchdict['name'] == 'disorder_prediction':
            JSONmotifs.append(
              {'start': matchdict['start'], 'end': matchdict['end'],
               'type': "disorder", 
                'colour': "#cccccc",
               'metadata': {'type': 'disorder', 'start': matchdict['start'], 'end': matchdict['end']}
              })

    graphicsJSON['motifs'] = JSONmotifs

    sortedJSONregions = sorted(JSONregions, key=lambda d: int (d['start']))
    color_index = 0
    color_domains = ["#2dcf00","#ff5353","#5b5bff","#ebd61d","#ba21e0","#ff9c42"]
    for region in sortedJSONregions:
      region['colour'] = color_domains[color_index % len(color_domains)]
      color_index += 1

    graphicsJSON['regions'] = sortedJSONregions
    return graphicsJSON


# xml2json = Cxml2json(xml_human_file)
# temp = xml2json.xmlForBaseUnp('Q9H091')
# print len(temp)
# print xml2json.isCanonical('Q9H091')
# print xml2json.isCanonical('Q9H091-1')
# print xml2json.isCanonical('Q9H091-2')
# print xml2json.isCanonical('Q9H091-3')

# graphicsJSON = xml2json.PFAMgraphicsJSONfromUnp('Q9H091')
# print 'Final JSON\n%s'%json.dumps(graphicsJSON)


