import cPickle as pickle
import logging
logger = logging.getLogger()

hpo_full_name_to_codes=pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Phenotype_Enrichment/hpo_full_name_to_codes.txt'))

hpo_code_to_all_names=pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Phenotype_Enrichment/hpo_code_to_all_names.txt'))

similar_codes=pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Phenotype_Enrichment/similar_codes.txt'))

hpo_gene_to_code=pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Phenotype_Enrichment/hpo_gene_to_code.txt'))

hpo_code_to_gene=pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Phenotype_Enrichment/hpo_code_to_gene.txt'))

hpo_name_to_code=pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Phenotype_Enrichment/hpo_processed_name_to_code.txt'))

hpo_code_to_name=pickle.load(open('/dors/capra_lab/users/mukhes1/Projects/Digenic_Classifier_Project/DATA_Processed/Phenotype_Enrichment/hpo_code_to_processed_name.txt'))

from scipy.stats import hypergeom

def process_code(code,keyword):
    new_code=code.replace(keyword+':','').strip().rstrip()
    return new_code

def process_name(full_name):
    new_name=full_name.lower().replace(' ','')
    return new_name

def get_full_name_from_short_name(name):
    full_name=hpo_code_to_name[hpo_name_to_code[name]]
    return full_name

def correct_tense(name):
    if name.lower()[-3:] == 'ies':
        name = name[:-3] + 'y'
    elif name.lower()[-1] == 's':
        name = name[:-1]
    return name

def clean_synonym(term):
    term=term.replace('synonym:','').strip().rstrip()
    import re
    a=[i.start() for i in re.finditer(r'"',term)]
    s=a[0]
    e=a[1]
    clean=term[s+1:e]
    return clean

def read_patient_phenotype(file=file):
    phenotype_list=open(file).read().split('\n')[:-1]
    return phenotype_list



def map_patient_phenotype_terms_to_HPO_codes(phenotype_list):
    # converting HPO terms provided for patient to codes in a dict

    #print '\nPhenotype Terms:'
    report='\nPhenotype Terms:'
    hpo_udn = []
    y = []

    for i in phenotype_list:
        p = process_name(i)
        if p in hpo_name_to_code:
            if i not in y:
                y.append(i)

            if hpo_name_to_code[p] not in hpo_udn:
                hpo_udn.append(hpo_name_to_code[p])
                #print get_full_name_from_short_name(p),' ('+hpo_name_to_code[p]+')'
                report=report+str(get_full_name_from_short_name(p))+'\t'+' ('+hpo_name_to_code[p]+')'+'\n'

        else:
            #if text does not match, deleting 's' or 'ies' at the end
            x=correct_tense(i)
            p=process_name(x)
            if p in hpo_name_to_code:
                if i not in y:
                    y.append(i)
                if hpo_name_to_code[p] not in hpo_udn:
                    hpo_udn.append(hpo_name_to_code[p])
                    #print get_full_name_from_short_name(p),'('+hpo_name_to_code[p]+')'
                    report=report+str(get_full_name_from_short_name(p))+'\t'+'('+hpo_name_to_code[p]+')'+'\n'

            else:
        # matching the patient term without 's' or 'ies' at the end is a substring of an existing HPO term
                for j in hpo_name_to_code:
                    if p in j:
                        if p in hpo_name_to_code:
                            if i not in y:
                                y.append(i)
                            if hpo_name_to_code[j] not in hpo_udn:
                                hpo_udn.append(hpo_name_to_code[j])
                                #print get_full_name_from_short_name(j),'('+hpo_name_to_code[j]+')'
                                report=report+str(get_full_name_from_short_name(j))+'\t'+'('+hpo_name_to_code[j]+')'+'\n'

    #if all fails, doing a fuzzy string matching, using all the words in the term, and choosing the best fit
            if i not in y:
                from fuzzywuzzy import process,fuzz
                l=[]

                for k in x.split(' '):
                    for j in hpo_name_to_code:
                        if k in j:
                            l.append(j)
                #get max 3 matches total
                for q in process.extract(p,l,scorer=fuzz.token_sort_ratio,limit=1):
                    if i not in y:
                        y.append(i)
                    if hpo_name_to_code[q[0]] not in hpo_udn:
                        hpo_udn.append(hpo_name_to_code[q[0]])
                        #print get_full_name_from_short_name(q[0]),'('+hpo_name_to_code[q[0]]+')'
                        report=report+str(get_full_name_from_short_name(q[0]))+'\t'+'('+hpo_name_to_code[q[0]]+')'+'\n'
    #print
    report=report+'\n'

    if [i for i in phenotype_list if i not in y] ==[]:
        #print 'All phenotypes mapped to HPO codes!'
        report=report+'\nAll phenotypes mapped to HPO codes!'+'\n'
    else:
        #print 'Not in database', [i for i in phenotype_list if i not in y]
        report=report+'Not in database'+'\t'+ str([i for i in phenotype_list if i not in y])

    return hpo_udn,report



def run_HPO_analysis(hpolist,genelist):

    report=''
    genes_to_code={}
    list_codes=[]
    genes_for_phenotype={}
    phenotypes_for_gene={}
    for i in genelist:
        c=[]
        if i in hpo_gene_to_code:



            for j in hpo_gene_to_code[i]:
                if j in hpolist:

                    if j not in c:
                        c.append(j)
                    if j not in genes_for_phenotype:
                        genes_for_phenotype[j] = []
                    if i not in genes_for_phenotype[j]:
                        genes_for_phenotype[j].append(i)
                if c==[]:
                    if j in similar_codes:
                        for k in similar_codes[j]:
                            if k in hpolist:
                                if k not in c:
                                    c.append(k)
                                if k not in genes_for_phenotype:
                                    genes_for_phenotype[k]=[]
                                if i not in genes_for_phenotype[k]:
                                    genes_for_phenotype[k].append(i)

            if c!=[]:
                if i not in phenotypes_for_gene:
                    phenotypes_for_gene[i]=[]
            for j in c:
                if j not in phenotypes_for_gene[i]:
                    phenotypes_for_gene[i].append(j)
                if j not in list_codes:
                    list_codes.append(j)
    #
    list_genes=[]
    for i in hpolist:
        g=[]
        if i in hpo_code_to_gene.keys():
            for j in hpo_code_to_gene[i]:
                if j in genelist:
                    if j not in g:
                        g.append(j)
                    if j not in phenotypes_for_gene:
                        phenotypes_for_gene[j]=[]
                    if i not in phenotypes_for_gene[j]:
                        phenotypes_for_gene[j].append(i)
            if g ==[]:
                if i in similar_codes:
                    for j in similar_codes[i]:
                        if j in hpo_code_to_gene:
                            for k in hpo_code_to_gene[j]:
                                if k in genelist:
                                    if k not in g:
                                       g.append(k)
                                    if k not in phenotypes_for_gene:
                                        phenotypes_for_gene[k] = []
                                    if i not in phenotypes_for_gene[k]:
                                        phenotypes_for_gene[k].append(i)

        if g!=[]:
            if i not in genes_for_phenotype:
                genes_for_phenotype[i]=[]
            for j in g:
                if j not in list_genes:
                    list_genes.append(j)
                if j not in genes_for_phenotype[i]:
                    genes_for_phenotype[i].append(j)

    for i in phenotypes_for_gene:
        if i not in list_genes:
            list_genes.append(i)

    from scipy.stats import hypergeom
    for i in genes_for_phenotype:
        if i not in list_codes:
            list_codes.append(i)
        pvalue=hypergeom.sf(len(genes_for_phenotype[i])-1,20000,len(hpo_code_to_gene[i]),len(genelist))
        if pvalue < 0.05:
            print 'Enriched phenotypes:\n',i,hpo_code_to_name[i],'(p value =',round(pvalue,4),')'
            gene_code_pairs=[]

    print
    print '\n{numgenes} ({percentgenes}%) UDN genes could describe {percent_phenotype} %  ({numcodes}/{totcodes}) ' \
          'of the phenotype \nGenes: {genes} '.format(genes=list_genes,numgenes=len(list_genes),totgenes=len(genelist),
                percentgenes=round(float(len(list_genes)) / len(genelist) * 100.,2),
                    percent_phenotype=round(float(len(list_codes)/float(len(hpolist)))*100.,2),
                            numcodes=len(list_codes),totcodes=len(hpolist))
    report=report+'{numgenes} ({percentgenes}%) UDN genes could describe {percent_phenotype} %  ({numcodes}/{totcodes}) ' \
          'of the phenotype :{genes} '.format(genes=list_genes,numgenes=len(list_genes),totgenes=len(genelist),
                percentgenes=round(float(len(list_genes)) / len(genelist) * 100.,2),
                    percent_phenotype=round(float(len(list_codes)/float(len(hpolist)))*100.,2),
                            numcodes=len(list_codes),totcodes=len(hpolist))+'\n'

    udn_gene_number_hpo=[]
    for i in list_genes:
        udn_gene_number_hpo.append(((len(phenotypes_for_gene[i]),phenotypes_for_gene[i]),i))
    print
    for i in sorted(udn_gene_number_hpo,key=lambda tup:tup[0],reverse=True):
        if i[0][0] > 1:
            print '{gene} covers {percent} % ({numgene}/{totgene}) of the HPO terms!\nPhenotypes:'.format(gene=i[1],
                percent=round(float(i[0][0])/len(hpolist)*100,2),numgene=i[0][0],totgene=len(hpolist))
            report=report+'\n{gene} covers {percent} % ({numgene}/{totgene}) of the HPO terms!'.format(gene=i[1],
                percent=round(float(i[0][0])/len(hpolist)*100,2),numgene=i[0][0],totgene=len(hpolist))+'\n'
    print
    for i in sorted(genes_for_phenotype,key=lambda i:len(genes_for_phenotype[i]),reverse=True):
        if len(genes_for_phenotype[i]) > 1:
            print 'Phenotype "{phenotype}" is covered by {percentgenes} % ({numgenes}/{totgenes}) of the total genes in UDN report!'\
            .format(phenotype=hpo_code_to_name[i],numgenes=len(genes_for_phenotype[i]),totgenes=len(genelist),
                    percentgenes=round(float(len(genes_for_phenotype[i]))/len(genelist)*100,2))
            report=report+'\nPhenotype "{phenotype}" is covered by {percentgenes} % ({numgenes}/{totgenes}) of the total genes in UDN report!'\
            .format(phenotype=hpo_code_to_name[i],numgenes=len(genes_for_phenotype[i]),totgenes=len(genelist),
                    percentgenes=round(float(len(genes_for_phenotype[i]))/len(genelist)*100,2))+'\n'

    codes_tot = hpolist
    smallest_subset_genes = []
    codes_uncovered=[]
    for i in sorted(udn_gene_number_hpo, key=lambda tup: tup[0], reverse=True):
        codes_uncovered = [j for j in codes_tot if j not in i[0][1]]
        if len(codes_uncovered) < len(codes_tot):
            smallest_subset_genes.append(i[1])
        codes_tot = codes_uncovered
    codes_covered = [j for j in hpolist if j not in codes_uncovered]

    print
    print '{numgenes} ({percentgenes}%) genes describe {percent_phenotype} % ' \
          'of the phenotype :{genes} '.format(genes=list_genes,numgenes=len(list_genes),totgenes=len(genelist),percentgenes=round(float(len(list_genes)) / len(genelist) * 100.,2),
                                    percent_phenotype=round(float(len(list_codes)/float(len(hpolist)))*100.,2))
    report=report+'\n{numgenes} ({percentgenes}%) genes describe {percent_phenotype} % ' \
          'of the phenotype :{genes} '.format(genes=list_genes,numgenes=len(list_genes),totgenes=len(genelist),percentgenes=round(float(len(list_genes)) / len(genelist) * 100.,2),
                                    percent_phenotype=round(float(len(list_codes)/float(len(hpolist)))*100.,2))+'\n'
    print
    print '{num_smallest_genes} ({percent_smallest_genes}%) genes cover ' \
          '{percent_phenotype} % ({numcodes}/{totcodes}) of phenotype : {name_genes} '.format(name_genes=smallest_subset_genes,
                num_smallest_genes=len(smallest_subset_genes),tot_genes=len(genelist),
                    percent_smallest_genes=round(float(len(smallest_subset_genes))/len(genelist)*100,2),
                        percent_phenotype=round(float(len(codes_covered)) / len(hpolist) * 100,2),
                                        numcodes=len(codes_covered),totcodes=len(hpolist))
    report=report+'\n{num_smallest_genes} ({percent_smallest_genes}%) genes cover ' \
          '{percent_phenotype} % ({numcodes}/{totcodes}) of phenotype : {name_genes} '.format(name_genes=smallest_subset_genes,
                num_smallest_genes=len(smallest_subset_genes),tot_genes=len(genelist),
                    percent_smallest_genes=round(float(len(smallest_subset_genes))/len(genelist)*100,2),
                        percent_phenotype=round(float(len(codes_covered)) / len(hpolist) * 100,2),
                                        numcodes=len(codes_covered),totcodes=len(hpolist))+'\n'
    print
    print 'Phenotype terms covered by the {} genes : '.format(len(smallest_subset_genes)), \
            sorted([hpo_code_to_name[i]for i in codes_covered ])
    report=report+'\nPhenotype terms covered by the {} genes : '.format(len(smallest_subset_genes))+'\t'+ \
            str(sorted([hpo_code_to_name[i]for i in codes_covered]))+'\n'

    gene_code_pairs=[]
    for i in genes_for_phenotype:
        for j in genes_for_phenotype[i]:
            value=len(hpo_code_to_gene[i])*len(hpo_gene_to_code[j])
            if (i,j,value) not in gene_code_pairs:
                gene_code_pairs.append((i,j,value))

    for i in phenotypes_for_gene:
        for j in phenotypes_for_gene[i]:
            value=len(hpo_gene_to_code[i])*len(hpo_code_to_gene[j])
            if (j,i,value) not in gene_code_pairs:
                gene_code_pairs.append((j,i,value))

    return report,list_genes,list_codes,gene_code_pairs


def run_best_match_analysis(list_genes,list_codes,gene_code_pairs,genelist,hpolist):
    list_pairs=''
    report=''
    print '\nBest Match for gene:'
    list_pairs=list_pairs+'\nBest Match for gene:'+'\n'
    best_match_code_gene_pair=[]
    for i in list_genes:
        x=[]
        for j in gene_code_pairs:
            if i == j[1]:
                x.append(j)
        best_match=sorted(x,key=lambda tup:tup[2])[0][0]
        value=sorted(x,key=lambda tup:tup[2])[0][2]
        print i,'-->',hpo_code_to_name[best_match]+' ('+best_match+')'
        list_pairs=list_pairs+i+'\t-->\t'+hpo_code_to_name[best_match]+'\t('+best_match+')'+'\n'
        if (best_match,i,value) not in best_match_code_gene_pair:
            best_match_code_gene_pair.append((best_match,i,value))


    print '\nBest Match for Phenotype:'
    list_pairs=list_pairs+'\nBest Match for Phenotype:'+'\n'
    best_fit_genes=[]
    for i in list_codes:
        x=[]
        for j in gene_code_pairs:
            if i == j[0]:
                x.append(j)
        best_match = sorted(x, key=lambda tup: tup[2])[0][1]
        value = sorted(x, key=lambda tup: tup[2])[0][2]
        print hpo_code_to_name[i]+'('+i+')','-->', best_match
        list_pairs=list_pairs+hpo_code_to_name[i]+'\t('+i+')'+'\t-->\t'+ best_match+'\n'
        if (i, best_match, value) not in best_match_code_gene_pair:
            best_match_code_gene_pair.append((i, best_match, value))
        if best_match not in best_fit_genes:
            best_fit_genes.append(best_match)

    best_score_code={}
    best_score_gene={}
    for i in best_match_code_gene_pair:
        if i[0] not in best_score_code:
            best_score_code[i[0]]=[]
        if i[1] not in best_score_code[i[0]]:
            best_score_code[i[0]].append(i[1])
        if i[1] not in best_score_gene:
            best_score_gene[i[1]]=[]
        if i[0] not in best_score_gene[i[1]]:
            best_score_gene[i[1]].append(i[0])

    print '\nRecalculating using BEST MATCH Criteria.....\n'
    report=report+'\nRecalculating using BEST MATCH Criteria.....\n'+'Phenotype matches, using best match:'+'\n'
    print 'Phenotype matches, using best match:'
    for i in sorted(best_score_code,key=lambda i:len(best_score_code[i]),reverse=True):
        if len(best_score_code[i])>1:
            print 'Phenotype "{name}" ({code}) is described by {percentgenes} % ({numgenes}/{totgenes}) of the UDN genes : {genes}'.format(
            name=hpo_code_to_name[i],code=i,percentgenes=round(float(len(best_score_code[i]))/len(genelist)*100,2),numgenes=len(best_score_code[i]),
            totgenes=len(genelist),genes=best_score_code[i])
            report=report+'Phenotype "{name}" ({code}) is described by {percentgenes} % ({numgenes}/{totgenes}) of the UDN genes : {genes}'.format(
            name=hpo_code_to_name[i],code=i,percentgenes=round(float(len(best_score_code[i]))/len(genelist)*100,2),numgenes=len(best_score_code[i]),
            totgenes=len(genelist),genes=best_score_code[i])+'\n'

    print '\nGene matches using Best match:'
    report=report+'\nGene matches using Best match:'+'\n'

    for i in sorted(best_score_gene,key=lambda i:len(best_score_gene[i]),reverse=True):
        if len(best_score_gene[i])>1:
            print 'Gene {name} covers {percentterms} % ({numterms}/{totterms}) of the phenotype {phenotype}'.format(name=i,
                percentterms=round(float(len(best_score_gene[i]))/len(hpolist)*100,2),numterms=len(best_score_gene[i]),totterms=len(hpolist),
                    phenotype=[hpo_code_to_name[j] for j in best_score_gene[i]])
            report=report+'Gene {name} covers {percentterms} % ({numterms}/{totterms}) of the phenotype {phenotype}'.format(name=i,
                percentterms=round(float(len(best_score_gene[i]))/len(hpolist)*100,2),numterms=len(best_score_gene[i]),totterms=len(hpolist),
                    phenotype=[hpo_code_to_name[j] for j in best_score_gene[i]])+'\n'
            if i not in best_fit_genes:
                best_fit_genes.append(i)

    print '\nSmallest set of Best Matched genes:'
    report=report+'\nSmallest set of Best Matched genes:'+'\n'
    print '{numgenes} or {percentgenes} % ({numgenes}/{totgenes}) describe {percentterms} % ({numterms}/{totterms}) of the Phenotype:\nGenes: {genes}\nPhenotypes: {phenotypes}'\
        .format(numgenes=len(best_fit_genes),percentgenes=round(float(len(best_fit_genes))/len(genelist)*100,2),totgenes=len(genelist),
                percentterms=round(float(len(best_score_code))/len(hpolist)*100,2),numterms=len(best_score_code),
                totterms=len(hpolist),genes=best_fit_genes,phenotypes=[hpo_code_to_name[i] for i in best_score_code])
    report=report+'{numgenes} or {percentgenes} % ({numgenes}/{totgenes}) describe {percentterms} % ({numterms}/{totterms}) of the Phenotype:\nGenes: {genes}\nPhenotypes: {phenotypes}'\
        .format(numgenes=len(best_fit_genes),percentgenes=round(float(len(best_fit_genes))/len(genelist)*100,2),totgenes=len(genelist),
                percentterms=round(float(len(best_score_code))/len(hpolist)*100,2),numterms=len(best_score_code),
                totterms=len(hpolist),genes=best_fit_genes,phenotypes=[hpo_code_to_name[i] for i in best_score_code])+'\n'
    return report,list_pairs


def find_genes_using_codes(hpolist):
    report=''
    phenotypes_for_gene={}
    genes_for_phenotype={}
    list_genes=[]
    for i in hpolist:
        g=[]
        if i in hpo_code_to_gene.keys():
            for j in hpo_code_to_gene[i]:
                if j not in g:
                    g.append(j)
                if j not in phenotypes_for_gene:
                    phenotypes_for_gene[j]=[]
                if i not in phenotypes_for_gene[j]:
                    phenotypes_for_gene[j].append(i)
            # if g==[]:
            if i in similar_codes:
                for j in similar_codes[i]:
                    if j in hpo_code_to_gene:
                        for k in hpo_code_to_gene[j]:
                            if k not in g:
                               g.append(k)
                            if k not in phenotypes_for_gene:
                                phenotypes_for_gene[k] = []
                            if i not in phenotypes_for_gene[k]:
                                phenotypes_for_gene[k].append(i)

        if g!=[]:
            if i not in genes_for_phenotype:
                genes_for_phenotype[i]=[]
            for j in g:
                if j not in list_genes:
                    list_genes.append(j)
                if j not in genes_for_phenotype[i]:
                    genes_for_phenotype[i].append(j)

    for i in phenotypes_for_gene:
        if i not in list_genes:
            list_genes.append(i)

    gene_code_pairs=[]
    for i in genes_for_phenotype:
        for j in genes_for_phenotype[i]:
            value=len(hpo_code_to_gene[i])*len(hpo_gene_to_code[j])
            if (i,j,value) not in gene_code_pairs:
                gene_code_pairs.append((i,j,value))

    for i in phenotypes_for_gene:
        for j in phenotypes_for_gene[i]:
            value=len(hpo_gene_to_code[i])*len(hpo_code_to_gene[j])
            if (j,i,value) not in gene_code_pairs:
                gene_code_pairs.append((j,i,value))
    best_match_code_gene_pair=[]
    list_pairs=''
    list_pairs = list_pairs + '\nBest Match for Phenotype:' + '\n'
    best_fit_genes = []
    for i in hpolist:
        x = []
        for j in gene_code_pairs:
            if i == j[0]:
                x.append(j)
        if x!=[]:
            best_match = sorted(x, key=lambda tup: tup[2])[0][1]
            value = sorted(x, key=lambda tup: tup[2])[0][2]
            print hpo_code_to_name[i] + '(' + i + ')', '-->', best_match
            list_pairs = list_pairs + hpo_code_to_name[i] + '\t(' + i + ')' + '\t-->\t' + best_match + '\n'
            if (i, best_match, value) not in best_match_code_gene_pair:
                best_match_code_gene_pair.append((i, best_match, value))
            if best_match not in best_fit_genes:
                best_fit_genes.append(best_match)
        else:
            print 'No matches for',i
    #
    # best_score_code = {}
    # best_score_gene = {}
    # for i in best_match_code_gene_pair:
    #     if i[0] not in best_score_code:
    #         best_score_code[i[0]] = []
    #     if i[1] not in best_score_code[i[0]]:
    #         best_score_code[i[0]].append(i[1])
    #     if i[1] not in best_score_gene:
    #         best_score_gene[i[1]] = []
    #     if i[0] not in best_score_gene[i[1]]:
    #         best_score_gene[i[1]].append(i[0])
    #
    # print '\nRecalculating using BEST MATCH Criteria.....\n'
    # report = report + '\nRecalculating using BEST MATCH Criteria.....\n' + 'Phenotype matches, using best match:' + '\n'
    # print 'Phenotype matches, using best match:'
    # for i in sorted(best_score_code, key=lambda i: len(best_score_code[i]), reverse=True):
    #     if len(best_score_code[i]) > 1:
    #         print 'Phenotype "{name}" ({code}) is described by {numgenes} : {genes}'.format(
    #             name=hpo_code_to_name[i], code=i,numgenes=len(best_score_code[i]),genes=best_score_code[i])
    #
    #         report = report + 'Phenotype "{name}" ({code}) is described by {numgenes}: {genes}'.format(
    #             name=hpo_code_to_name[i], code=i,numgenes=len(best_score_code[i]),genes=best_score_code[i]) + '\n'


    return report,list_genes,list_pairs


def find_codes_using_genes(genelist):
    phenotypes_for_gene={}
    genes_for_phenotype={}
    list_codes=[]
    for i in genelist:
        c=[]
        if i in hpo_gene_to_code:
            for j in hpo_gene_to_code[i]:
                if j not in c:
                    c.append(j)
                if j not in genes_for_phenotype:
                    genes_for_phenotype[j] = []
                if i not in genes_for_phenotype[j]:
                    genes_for_phenotype[j].append(i)

                if j in similar_codes:
                    for k in similar_codes[j]:
                        if k not in c:
                            c.append(k)
                        if k not in genes_for_phenotype:
                            genes_for_phenotype[k]=[]
                        if i not in genes_for_phenotype[k]:
                            genes_for_phenotype[k].append(i)

            if c!=[]:
                if i not in phenotypes_for_gene:
                    phenotypes_for_gene[i]=[]
            for j in c:
                if j not in phenotypes_for_gene[i]:
                    phenotypes_for_gene[i].append(j)
                if j not in list_codes:
                    list_codes.append(j)



    for i in genes_for_phenotype:
        if i not in list_codes:
            list_codes.append(i)


    gene_code_pairs=[]
    for i in genes_for_phenotype:
        for j in genes_for_phenotype[i]:
            value=len(hpo_code_to_gene[i])*len(hpo_gene_to_code[j])
            if (i,j,value) not in gene_code_pairs:
                gene_code_pairs.append((i,j,value))

    for i in phenotypes_for_gene:
        for j in phenotypes_for_gene[i]:
            value=len(hpo_gene_to_code[i])*len(hpo_code_to_gene[j])
            if (j,i,value) not in gene_code_pairs:
                gene_code_pairs.append((j,i,value))

    list_pairs=''

    print '\nBest Match for gene:'
    list_pairs=list_pairs+'\nBest Match for gene:'+'\n'
    best_match_pairs=[]
    for i in list_codes:
        x=[]
        for j in gene_code_pairs:
            if i == j[0]:
                x.append(j)
        best_match = sorted(x, key=lambda tup: tup[2])[0][1]
        value = sorted(x, key=lambda tup: tup[2])[0][2]
        print hpo_code_to_name[i]+'('+i+')','-->', best_match
        list_pairs=list_pairs+hpo_code_to_name[i]+'\t('+i+')'+'\t-->\t'+ best_match+'\n'
        if (i, best_match, value) not in best_match_pairs:
            best_match_pairs.append((i, best_match, value))
    return list_codes,list_pairs,best_match_pairs

def get_all_phenotype_linked_genes(hpo_udn):
    hit_genes = []
    sig_hits={}
    for i in hpo_udn:
        if i in hpo_code_to_gene:
            hit_genes = list(set(hit_genes).union(hpo_code_to_gene[i]))

    tot_codes = len(hpo_code_to_gene)
    all_udn_codes = len(hpo_udn)

    for i in hit_genes:
        if i in hpo_gene_to_code:
            udn_codes_for_gene=len(list(set(hpo_gene_to_code[i]).intersection(hpo_udn)))-1
            all_codes_for_gene=len(hpo_gene_to_code[i])

            pvalue=hypergeom.sf(udn_codes_for_gene,tot_codes,all_codes_for_gene,all_udn_codes)
            if pvalue <=0.01:
                sig_hits[i]=[0,[]]
                sig_hits[i][0]=pvalue
                sig_hits[i][1]=[hpo_code_to_name[x] for x in list(set(hpo_gene_to_code[i]).intersection(hpo_udn))]
    return hit_genes,sig_hits

def get_related_codes(hpo_udn):
    codes={}
    for i in hpo_udn:
        if i in similar_codes:
            codes[i]=[i] + similar_codes[i]
            # codes=list(set(codes).union(list(set(similar_codes[i]).union([i]))))
    return codes
