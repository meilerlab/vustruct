import cPickle as pickle
from scipy.stats import hypergeom
import os

DATAROOT='/dors/capra_lab/projects/psb_collab/UDN/data_for_systems_biology_analysis/'

upstream_regulators=pickle.load(open(os.path.join(DATAROOT,'Regulation_Analysis/Regulators_of_genes.txt')))
downstream_targets=pickle.load(open(os.path.join(DATAROOT,'Regulation_Analysis/Targets_of_genes.txt')))

def get_enriched_regulators(genes):
    count=[]
    all_reg={}
    sig_regs={}
    for i in genes:
        if i in upstream_regulators:
            count.append(i)
            for j in upstream_regulators[i]:
                if j not in all_reg:
                    all_reg[j]=[]
                if i not in all_reg[j]:
                    all_reg[j]=list(set(all_reg[j]).union([i]))

    for i in all_reg:
        if len(all_reg[i]) > 1:
            target_udn_genes=len(all_reg[i])
            all_targets=len(downstream_targets[i])
            udn_genes=len(list(set(upstream_regulators.keys()).intersection(genes)))
            tot_genes=len(upstream_regulators)
            pvalue=hypergeom.sf(target_udn_genes-1,tot_genes,udn_genes,all_targets)
            if pvalue <=0.05:
                sig_regs[i]=[all_reg[i],pvalue]

    return sig_regs


def get_enriched_targets(genes):

    count = []
    all_tgts = {}
    sig_tgts={}
    for i in sorted(genes):
        if i in downstream_targets:
            count.append(i)
            for j in downstream_targets[i]:
                if j not in all_tgts:
                    all_tgts[j] = []
                if i not in all_tgts[j]:
                    all_tgts[j] = list(set(all_tgts[j]).union([i]))

    for i in all_tgts:
        if len(all_tgts[i]) > 1 :
            reg_udn_genes=len(all_tgts[i])
            all_regs=len(upstream_regulators[i])
            udn_genes=len(list(set(downstream_targets.keys()).intersection(genes)))
            tot_genes=len(downstream_targets)
            pvalue=hypergeom.sf(reg_udn_genes-1,tot_genes,udn_genes,all_regs)
            if pvalue <=0.05:
                sig_tgts[i]=[all_tgts[i],pvalue]

    return sig_tgts
