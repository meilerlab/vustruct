import cPickle as pickle
import os
# coexpressdb=pickle.load(open('/Users/souhrid/Coexpress_db_dict.txt'))
# coex_list=open('/Users/souhrid/Coexpressdb_list.txt').read().split('\n')[:-1]
# gene_id=pickle.load(open('/Users/souhrid/Entrez_id_gene_id_conversion.txt'))
# DIDA_list = open ('/Users/souhrid/Downloads/DIDA_pairs_list.csv').read().split('\n')[1:-1]
# all_genes=open('/Users/souhrid/All_genes.txt').read().split('\n')[:-1]


DATAROOT='/dors/capra_lab/projects/psb_collab/UDN/data_for_systems_biology_analysis/'

coexpress_dict=pickle.load(open(os.path.join(DATAROOT,'Coexpression_Analysis/Coexpress_db_dict_select.txt')))
