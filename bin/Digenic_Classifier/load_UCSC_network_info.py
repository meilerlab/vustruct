import cPickle as pickle
import networkx as nx
import logging

logger = logging.getLogger()

logger.info( 'Importing Protein-Protein Network...' )
G_ppi=nx.read_gpickle('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Network_Analysis/UCSC_ppi_network_new.gpickle')

dists_ppi=pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Network_Analysis/PPI_network_allpairs_shortest_path_lengths_new_Feb_22_18.txt'))

#rev_dist_ppi=pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Network_Analysis/PPI_rev_dist_dict.txt'))

logger.info( 'Importing Protein-Protein network DONE!! (G_ppi contains %d graph edges)'%G_ppi.size() )

logger.info( 'Importing Pathway Network...' )

G_pwy=nx.read_gpickle('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Network_Analysis/UCSC_pwy_network_new.gpickle')

dists_pwy=pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Network_Analysis/PWY_network_allpairs_shortest_path_lengths_new_Feb_22_18.txt'))

#rev_dist_pwy=pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Network_Analysis/PWY_rev_dist_dict.txt'))

logger.info( 'Importing Pathway Network DONE!! (G_pwy contains %d graph edges)'%G_pwy.size() )

logger.info( 'Importing Literature mined Network...' )

G_txt=nx.read_gpickle('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Network_Analysis/UCSC_txt_network_new.gpickle')

dists_txt=pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Network_Analysis/Txt_network_allpairs_shortest_path_lengths_new_Feb_22_18.txt'))

#rev_dist_txt=pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Network_Analysis/Txt_rev_dist_dict.txt'))

logger.info( 'Importing Literature mined Network DONE!! (G_txt contains %d graph edges)'%G_txt.size() )

logger.info( 'Importing Full Network Connections details...' )

codes_rec=pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Network_Analysis/UCSC_int_code_details.txt'))

logger.info( 'Full Network Connection details importing.. DONE!! len(codes_rec)=%d'%len(codes_rec) )
