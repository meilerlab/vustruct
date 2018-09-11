import cPickle as pickle
import networkx as nx
import logging

logger = logging.getLogger()

logger.info( 'Importing Protein-Protein Network...' )
G_ppi=nx.read_gpickle('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Network_Analysis/UCSC_ppi_network_new.gpickle')
print str(G_ppi)

if 'asdf' in G_ppi:
  print 'found'
else:
  print 'not found'

