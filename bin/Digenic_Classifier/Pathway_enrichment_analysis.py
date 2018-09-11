import cPickle as pickle
from scipy.stats import hypergeom

reactome_gene_to_path_codes = pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Pathway_Enrichment/reactome_gene_to_path_codes.txt'))

reactome_path_to_genes = pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Pathway_Enrichment/reactome_path_to_genes.txt'))

reactome_path_code_to_name = pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Pathway_Enrichment/reactome_path_code_to_name.txt'))

reactome_tot_codes_in_path = pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Pathway_Enrichment/reactome_tot_codes_in_path.txt'))

kegg_list_genes = pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Pathway_Enrichment/kegg_list_genes.txt'))

kegg_full_gene_name = pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Pathway_Enrichment/kegg_full_gene_name.txt'))

kegg_gene_name = pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Pathway_Enrichment/kegg_gene_name.txt'))

kegg_gene_to_path_codes = pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Pathway_Enrichment/kegg_gene_to_path_codes.txt'))

kegg_path_to_genes = pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Pathway_Enrichment/kegg_path_to_genes.txt'))

kegg_path_code_to_name = pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Pathway_Enrichment/kegg_path_code_to_name.txt'))

def get_enriched_pathways(genes,cutoff,pv):
    list_paths = []
    paths_with_gene_reactome = {}
    for i in genes:
        if i in reactome_gene_to_path_codes:
            for j in reactome_gene_to_path_codes[i]:
                if j not in paths_with_gene_reactome:
                    paths_with_gene_reactome[j] = []
                if i not in paths_with_gene_reactome[j]:
                    paths_with_gene_reactome[j].append(i)

    list_genes_reactome = []
    for i in paths_with_gene_reactome:
        a = len(paths_with_gene_reactome[i]) / float(len(reactome_path_to_genes[i]))
        b = len(paths_with_gene_reactome[i]) / float(len(genes))
        list_genes_reactome.append(
            (a * b, a * 100, b * 100, len(paths_with_gene_reactome[i]), reactome_tot_codes_in_path[i], i))

    paths_with_gene_kegg = {}
    for i in genes:
        if i in kegg_gene_to_path_codes:
            for j in kegg_gene_to_path_codes[i]:
                if j not in paths_with_gene_kegg:
                    paths_with_gene_kegg[j] = []
                if i not in paths_with_gene_kegg[j]:
                    paths_with_gene_kegg[j].append(i)

    list_genes_kegg = []
    for i in paths_with_gene_kegg:
        a = len(paths_with_gene_kegg[i]) / float(len(kegg_path_to_genes[i]))
        b = len(paths_with_gene_kegg[i]) / float(len(genes))
        list_genes_kegg.append((a * b, a * 100, b * 100, len(paths_with_gene_kegg[i]), len(kegg_path_to_genes[i]), i))

    report = []
    record={}
    c = 't'
    for i in sorted(list_genes_reactome, key=lambda tup: tup[0], reverse=True):
        pvalue = hypergeom.sf(len(paths_with_gene_reactome[i[5]]) - 1, 20000, reactome_tot_codes_in_path[i[5]],
                              len(genes))
        if len(paths_with_gene_reactome[i[5]]) > cutoff and pvalue < pv:
            ans = '\nPathway {name} ({code}) is enriched {pathpercent} % ({inpath}/{totpath}) (p value= {pvalue})  by {genepercent} % ' \
                  '({inudn}/{totudn}) of the UDN genes : {listgenes}'.format(name=reactome_path_code_to_name[i[5]],
                                                                             code=i[5],
                                                                             pathpercent=round(i[1], 2),
                                                                             inpath=len(paths_with_gene_reactome[i[5]]),
                                                                             totpath=reactome_tot_codes_in_path[i[5]],
                                                                             genepercent=round(i[2], 2),
                                                                             inudn=len(paths_with_gene_reactome[i[5]]),
                                                                             totudn=len(genes),
                                                                             listgenes=paths_with_gene_reactome[i[5]],
                                                                             pvalue=round(pvalue, 4))
            list_paths.append((reactome_path_code_to_name[i[5]],i[5],paths_with_gene_reactome[i[5]],round(pvalue,3)))
            report.append((1 / i[0] * pvalue, ans))
            record[i[5]]=paths_with_gene_reactome[i[5]]
            c = 'f'
    # if c == 't':
    #     report.append((0, 'No pathways Enriched!!'))

    c = 't'
    for i in sorted(list_genes_kegg, key=lambda tup: tup[0], reverse=True):
        pvalue = hypergeom.sf(len(paths_with_gene_kegg[i[5]]) - 1, 20000, len(kegg_path_to_genes[i[5]]), len(genes))
        if len(paths_with_gene_kegg[i[5]]) > cutoff and pvalue < pv:
            ans = '\nPathway {name} ({code}) is enriched {pathpercent} % ({inpath}/{totpath}) (p value = {pvalue}) by {genepercent} % ' \
                  '({inudn}/{totudn}) of the UDN genes : {listgenes}'.format(name=kegg_path_code_to_name[i[5]],
                                                                             code=i[5],
                                                                             pathpercent=round(i[1], 2),
                                                                             inpath=len(paths_with_gene_kegg[i[5]]),
                                                                             totpath=len(kegg_path_to_genes[i[5]]),
                                                                             genepercent=round(i[2], 2),
                                                                             inudn=len(paths_with_gene_kegg[i[5]]),
                                                                             totudn=len(genes),
                                                                             listgenes=paths_with_gene_kegg[i[5]],
                                                                             pvalue=round(pvalue, 4))
            list_paths.append((kegg_path_code_to_name[i[5]],i[5],paths_with_gene_kegg[i[5]],round(pvalue,3)))
            report.append((1 / i[0] * pvalue, ans))
            record[i[5]] = paths_with_gene_kegg[i[5]]
            c = 'f'
    # if c == 't':
    #     report.append((0, 'No pathways Enriched!!'))



    return list_paths,report,record

