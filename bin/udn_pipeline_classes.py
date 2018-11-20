import datetime
import glob
import gzip
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import re
import shutil
import subprocess
import sys
import pprint
from collections import OrderedDict
from numpy import array, linalg, mean
from shutil import copyfile
from timeit import default_timer
from urllib import urlopen

import logging

# In case we are called by a client that has not configured logging,
# NullHandler ensures that does not fail
class NullHandler(logging.Handler):
  def emit(self, record):
    pass

LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(NullHandler())

sys.path.insert(0, '/dors/meilerlab/home/sliwosgr/gregit/')
from prepare import cleaner,query_modres

AA = {"R":"ARG","ARG":"R","H":"HIS","HIS":"H","K":"LYS","LYS":"K",
      "D":"ASP","ASP":"D","E":"GLU","GLU":"E","S":"SER","SER":"S",
      "T":"THR","THR":"T","N":"ASN","ASN":"N","Q":"GLN","GLN":"Q",
      "C":"CYS","CYS":"C","U":"SEC","SEC":"U","G":"GLY","GLY":"G",
      "P":"PRO","PRO":"P","A":"ALA","ALA":"A","V":"VAL","VAL":"V",
      "I":"ILE","ILE":"I","L":"LEU","LEU":"L","M":"MET","MET":"M",
      "F":"PHE","PHE":"F","Y":"TYR","TYR":"Y","W":"TRP","TRP":"W"
      }

def avg(atoms):
    coords = [0,0,0]
    for item in atoms:
        coords[0] += item[0]
        coords[1] += item[1]
        coords[2] += item[2]
    return [x/len(atoms) for x in coords]

def distance(a, b):
    aa = array(a)
    bb = array(b)
    return linalg.norm(aa-bb)    

#Returns [True,isopdbname,to_pdb indexing] if successful
def isolate_chain(pdb,chain):
    with open(pdb) as infile:
        pdblines = infile.readlines()
    to_pdb = ['NA']
    isopdb = ''
    pdbtext = os.path.basename(os.path.splitext(pdb)[0])
    isopdbname = "%s_%s_iso.pdb" % (pdbtext,chain)
    atomnum = 1
    for line in pdblines:
        if len(line)>20 and line[:4] =='ATOM' and line[21]==chain:
            try:
                resnum = int(line[22:26])
                aa = AA[line[17:20]]
            except ValueError:
                return[False, "NA","isolate_chain: Error parsing %s" % pdb]
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
    return [True, isopdbname, to_pdb]

def renum_pdb(pdb,chain_check=None,min_chains=None,write=True):
    with open(pdb) as infile:
        pdblines = infile.readlines()
    to_pdb = ['NA']
    renumpdb = ""
    pdbtext = "%s_%s" % (os.path.basename(os.path.splitext(pdb)[0]),os.path.splitext(pdb)[1][-1])    
    renumpdbname = "%s_renum.pdb" % pdbtext
    atomnum = 1
    allchains = set()
    if chain_check is not None:
        chain_seen = False
    else:
        chain_seen = True
    for line in pdblines:
        if len(line)>20 and line[:4] =='ATOM':
            chain = line[21]
            resnum = int(line[22:26])
            fullres = chain+str(resnum)
            if to_pdb[-1] != fullres:
                to_pdb.append(fullres)
            new_atomnum = "%5d" % atomnum
            new_resnum = len(to_pdb)-1
            new_line = line[:6]+new_atomnum+line[11:22]+"%4d" % new_resnum + line[26:]
            allchains.add(chain)
            renumpdb += new_line
            atomnum += 1
            if chain==chain_check:
                chain_seen = True
    if not chain_seen:
        return [False,"NA","Chain %s not seen in file" % chain_check,"NA"]
    elif min_chains is not None and len(allchains)<min_chains:
        return [False,"NA","File has less than minimum %s chains" % min_chains,"NA"]
    elif write:
        with open(renumpdbname,'w') as outfile:
            outfile.write(renumpdb)
            outfile.write("TER\n")
    return [True,renumpdbname,to_pdb,allchains]

class predict_ss(object):
    def __init__(self,fastafile,mutations,scripts_path='/dors/meilerlab/home/sliwosgr/udn_pipeline/predict_ss/runss_scripts'):
        self._scripts = scripts_path.rstrip("/")
        self._fastafile = fastafile
        self._mutation = mutations[0]
        self._residue = self._mutation[1]
        self._predictors = {'PP':['path_to_psipred2','psipred_ss2'],'Jufo':['pathtojufo','jufo9d_ss']}
        self._color_map = {'H':'yellow','E':'orange','C':'cyan','unknown':'green'}
        self._visbuffer = 20
        self._vis = []
        self._final_preds = []
        with open(self._fastafile) as infile:
            self._fasta = "".join(x.strip() for x in infile.readlines() if x[0]!='>')
        self._mutfasta = self._fasta[:self._mutation[1]-1]+self._mutation[2]+self._fasta[self._mutation[1]:]

        if len(mutations)>1:
            print "Warning, predict_ss currently only supports 1 mutation, using only first mutation (%s)" % "".join(str(x) for x in self._mutation)
        LOGGER.info("predict_ss.__init__  completed after setting:\n%s"%pprint.pformat(self.__dict__))
        
    def failure(self, message):
        print message
        return [False,message]


    def make_pic(self, visoutfile):
        if len(self._final_preds) == 0:
            "No prediction data to generate picture"
            return False

        residue = self._residue
        mutation = self._mutation
        disprange = range(max(1,residue-self._visbuffer),min(len(self._fasta),residue+self._visbuffer)+1)
        npreds = len(self._predictors)
        ixrow = ['res','aa','mut_aa']       
        for mit in sorted(self._predictors.keys()):
            ixrow += [mit,"%s_M" % mit]
        self._vis = [ixrow]
        self._vis += [x[:(npreds*2+3)] for x in self._final_preds[1:] if int(x[0]) in disprange]            
        vis = self._vis
        #Each x value is residue number and each y value is 1
        #Colors of each bar is determined by COLOR_MAP
        #Except for variant, which is white except for residue that corresponds to selection
        x = disprange
        y = [1]*len(x)
        if self._mutation is None:
            nplts = len(self._predictors)+2
        else:
            nplts = len(self._predictors)*3+1
            mutcol = {}
            for item in self._predictors:
                curix = vis[0].index("%s_M" % item)
                mutcol[item] = [self._color_map[it[curix]] for it in vis[1:]]
        h = (nplts+1)*.4
        w = max(len(disprange)/2.5,1)
    
        #Leave space for xlabel if buffer too small %TODO: Improve this nonsense
        if len(disprange) < 5:
            lb = .5
            rb = None
        else:
            rb = lb = None

        varcol = ['white']*len(disprange)
        #Color selection red, is index 19 unless res number is less than VIS_BUFFER.
        if residue<self._visbuffer+1:
            varcol[residue-1] = 'red'
        else:
            varcol[self._visbuffer] = 'red'
        wtcol = {}
        for item in self._predictors:
            curix = vis[0].index(item)
            wtcol[item] = [self._color_map[it[curix]] for it in vis[1:]]
    
        #Adjust pic size based on # plots and # of residues to display
        plt.rcParams["figure.figsize"] = [w,h]
        fig, axs = plt.subplots(nplts,1,sharex=True,sharey=True)
        if mutation is None:
            labels = ['']+sorted(self._predictors.keys(),reverse=True)+['']
            preds = [varcol]+[wtcol[it] for it in sorted(wtcol.keys(),reverse=True)]+[varcol]
        else:
            labels = ['']
            for it in sorted(self._predictors.keys(),reverse=True):       
                labels += ["%s_Mut" % it,"%s_WT" % it,'']
            labels += ['']
            preds = [varcol]
            for it in sorted(self._predictors.keys(),reverse=True):
                preds += [mutcol[it],wtcol[it],varcol]
        
        plt.subplots_adjust(hspace=None,bottom=.2,left=lb,right=rb)
        axs[-1].axis([min(x)-.5,max(x)+.5,0,1])
        axs[-1].set_xticks(x,1)
        axs[-1].set_yticks([])
        if disprange[-1]>999 or self._visbuffer<2:
            axs[-1].set_xticklabels(disprange,fontsize=8)
    
        for item in xrange(len(axs)):
            axs[item].set_ylabel(labels[item],rotation=0,va='center',ha='right')
            axs[item].tick_params(axis=u'both', which=u'both',length=0)
            axs[item].bar(disprange,y,1,color=preds[item],align='center')
    
        plt.xticks(x)
        plt.savefig(visoutfile)
        if os.path.isfile(visoutfile):
            return True
        else:
            return False        

    def ss(self,sequencefile):
        p1 = os.path.splitext(os.path.basename(sequencefile))[0]
        predfiles = {}
        predraw = {}
        for it in self._predictors:
            predfiles[it] = "%s.%s" % (os.path.splitext(sequencefile)[0],self._predictors[it][-1])
        subprocess.call("%s/runss_local %s" % (self._scripts,sequencefile),shell=True)
        for it in predfiles:
            if not os.path.isfile(predfiles[it]): [False,"SS fail: necessary output file %s not found." % predfiles[it]]
            with open(predfiles[it]) as infile:
                predraw[it] = [x.strip().split() for x in infile.readlines() if len(x.strip().split())>1 and x.strip().split()[0].isdigit()]
        return [True,predraw]

    def test_ss(self,sequencefile):
        p1 = os.path.splitext(os.path.basename(sequencefile))[0]
        predfiles = {}
        predraw = {}
        for it in self._predictors:
            predfiles[it] = "%s.%s" % (os.path.splitext(sequencefile)[0],self._predictors[it][-1])
#        subprocess.call("%s/runss_local %s" % (self._scripts,sequencefile),shell=True)
        for it in predfiles:
            if not os.path.isfile(predfiles[it]): [False,"SS fail: necessary output file %s not found." % predfiles[it]]
            with open(predfiles[it]) as infile:
                predraw[it] = [x.strip().split() for x in infile.readlines() if len(x.strip().split())>1 and x.strip().split()[0].isdigit()]
        return [True,predraw]                                           

    def run(self,testing = False):
        mutationtext = "".join(str(x) for x in self._mutation)
        p1 = os.path.splitext(os.path.basename(self._fastafile))[0]
        mutfilename = "%s_%s.fasta" % (p1,mutationtext)       
        with open(mutfilename,'w') as fout:
            fout.write(">%s_%s\n" % (p1,mutationtext))
            fout.write(self._mutfasta)
            fout.write("\n")
        if not testing:
            wtoutcome, wt_raw = self.ss(self._fastafile)
        else:
            wtoutcome,wt_raw = self.test_ss(self._fastafile)
        if not wtoutcome:
            return self.failure(wt_raw)
        if not testing:
            mutoutcome, mut_raw = self.ss(mutfilename)
        else:
            mutoutcome,mut_raw = self.test_ss(mutfilename)
        if not mutoutcome:
            return self.failure(mut_raw)
            
        label = ['RESIDUE','WT_AA','MUT_AA']
        for it in sorted(self._predictors.keys()):
            label += ["WT_%s_SS" % it, "MUT_%s_SS" % it]
        for it in sorted(self._predictors.keys()):
            label += ["WT_%s_CONF" % it, "MUT_%s_CONF" % it]
        for it in sorted(self._predictors.keys()):
            label += ["%s_change" % it]
        final_preds = [label]

        for x in xrange(len(self._fasta)):
            ref = sorted(self._predictors.keys())[0]
            res = int(wt_raw[ref][x][0])
            aa = wt_raw[ref][x][1]
            wtp = []
            wtc = []
            for it in sorted(self._predictors.keys()):
                wtp.append(wt_raw[it][x][2])
                wtc.append(max(float(wt_raw[it][x][3]),float(wt_raw[it][x][4]),float(wt_raw[it][x][5])))
    
            mutp = []
            mutc = []
            for it in sorted(self._predictors.keys()):
                mutp.append(mut_raw[it][x][2])
                mutc.append(max(float(mut_raw[it][x][3]),float(mut_raw[it][x][4]),float(mut_raw[it][x][5])))
            mut_change = []
            for it in xrange(len(wtp)):
                mut_aa = mut_raw[ref][x][1]
                if wtp[it]==mutp[it]:
                    mut_change.append('-')
                else:
                    mut_change.append("%s->%s" % (wtp[it],mutp[it]))               
            wtmut = [wmit for sublist in zip(wtp,mutp) for wmit in sublist]
            wtmutcon = [wmit for sublist in zip(wtc,mutc) for wmit in sublist]
            final_preds.append([str(res),aa,mut_aa]+wtmut+["%.3f" % xp for xp in wtmutcon]+mut_change)
        self._final_preds = final_preds
        return [True, final_preds]      
        
class ddg_monomer(object):
    def __init__(self,pdb,chain,mutations,rosetta_path='/dors/meilerlab/apps/rosetta/rosetta-3.7/main/source/bin'):
        self._rosetta = rosetta_path.rstrip("/")
        self._pdb = pdb
        self._chain = chain
        self._mutations = mutations
#        self._root = root_path
        self._break = "".join('*'*20)
        self._success = "/dors/meilerlab/home/sliwosgr/udn_pipeline/ddg_calculator/history/success_log.tab"
        self._fail = "/dors/meilerlab/home/sliwosgr/udn_pipeline/ddg_calculator/history/fail_log.tab"
        self._applications = ['minimize_with_cst.default.linuxgccrelease','per_residue_energies.linuxgccrelease','ddg_monomer.linuxgccrelease','score_jd2.linuxgccrelease']
        self._timers = ['minimize','rescore','ddg_monomer']
        self._to_pdb = []
        LOGGER.info("ddg_mmonomer.__init__ completed after setting:\n%s"%pprint.pformat(self.__dict__))
#        with open(self._pdb) as infile:
#            self._pdblines = infile.readlines()
                
    def log_success(self,times):
        mutationtext = ",".join("".join(str(x) for x in y) for y in self._mutations)
        nres = str(len(self._to_pdb)-1)
        curdate = datetime.datetime.now().strftime("%d_%m_%y")       
        curtime = datetime.datetime.now().strftime("%H_%M")
        wd = os.getcwd()
        if not os.path.isfile(self._success):
              with open(self._success,'w') as outfile:
                header = ['date','time','pdbfile','chain','nres','mutation']+self._timers+['path']
                outfile.write("\t".join(header))
                outfile.write("\n")
        with open(self._success,'a') as outfile:
            line = [curdate,curtime,self._pdb,self._chain,nres,mutationtext]+\
            ["%.2f" % (times[x]/60.0) if x in times else 'NaN' for x in self._timers] + [wd]
            outfile.write("\t".join(line))
            outfile.write("\n")

    def log_failure(self,message):
        wd = os.getcwd()
        mutationtext = ",".join("".join(str(x) for x in y) for y in self._mutations)
        curdate = datetime.datetime.now().strftime("%d_%m_%y")       
        curtime = datetime.datetime.now().strftime("%H_%M")
        if not os.path.isfile(self._fail):
            with open(self._fail,'w') as outfile:
                header = ['date','time','pdbfile','chain','mutation','message','path']
                outfile.write("\t".join(header))
                outfile.write("\n")
        with open(self._fail,'a') as outfile:
            line = [curdate,curtime,self._pdb,self._chain,mutationtext,message,wd]
            outfile.write("\t".join(line))
            outfile.write("\n")
            
    def failure(self, dumps, message):
        LOGGER.warning(message)
        if dumps is not None:
            for item in dumps:
                LOGGER.warning("Dumping %s to %s",item[0],item[1])
                with open(item[1],'w') as outfile:
                    outfile.write(item[2])
        #self.log_failure(message)
        return [False,message]


#    def protocols(self, selection=16):
        #Protocols defined by kellogg et al. 2010
#        prots = {
#        16: ["-in:file:fullatom\n-ddg:weight_file soft_rep_design\n-in:file:fullatom\n-ddg:weight_file soft_rep_design\n-ddg:local_opt_only false\n-ddg:min_cst true\n-ddg:mean false\n-ddg:min true\n-ddg:sc_min_only false\n-ddg:ramp_repulsive true\n-ddg:minimization_scorefunction talaris2014\n-ignore_zero_occupancy false\n")
#        6: ["-ddg:weight_file soft_rep_design\n-ddg:local_opt_only false\n-ddg:min_cst false\n-ddg:mean true
    def get_ddg(self,quality,silent=False):
        suffix = datetime.datetime.now().strftime("%d%m%y%H%M%S")       

        with open("ddg_predictions.out") as infile:
            muts = [x.strip().split()[1] for x in infile.readlines()[1:] if len(x.strip().split()) > 1]
        if not os.path.exists('top_models'):
            os.makedirs('top_models')
        top_three = {'wt': []}            

        if not silent:
            wt_files = glob.glob("repacked*.pdb")
            mut_files = glob.glob("mut*.pdb")
            modellist = "models_%s.ls" % suffix
            rescorefile = "ddg_%s.sc" % suffix
            with open(modellist,'w') as outfile:
                for x in wt_files:
                    outfile.write("\n".join(wt_files+mut_files))
                    outfile.write("\n")
            command = "score_jd2.linuxgccrelease -l %s -score:weights talaris2014 -out:file:scorefile %s" % (modellist,rescorefile)
            LOGGER.info("self._break = %s",self._break)
            LOGGER.info("Scoring all ddg models with command:\n%s/%s" , self._rosetta,command)
            runscore = subprocess.Popen("%s/%s" % (self._rosetta, command), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (scout,scerr) = runscore.communicate()
            try:
                assert os.path.isfile(rescorefile)
            except AssertionError:
                scoutfile = "modscorefail%s.out" % suffix
                scerrfile = "modscorefail%s.err" % suffix
                return self.failure([["rescore output",scoutfile,scout],["rescore stderr",scerrfile,scerr]],"Error scoring ddg models, scorefile is missing.")

            with open(rescorefile) as infile:
                raw_rescore = [x.strip().split() for x in infile.readlines() if x.strip()!='']

            try:
                assert len(raw_rescore) > 0
            except AssertionError:
                scoutfile = "modscorefail%s.out" % suffix
                scerrfile = "modscorefail%s.err" % suffix
                return self.failure([["rescore output",scoutfile,scout],["rescore stderr",scerrfile,scerr]],"Error scoring ddg models, scorefile is empty.")
        else:
            raw_rescore = []
            with open("wt_.out") as infile:
                raw_rescore += [x.strip().split() for x in infile.readlines() if x.strip().split()[-1]!='description' and x.strip().split()[0]=='SCORE:']
            for item in muts:
                with open("mut_%s.out" % item) as infile:
                    raw_rescore += [x.strip().split() for x in infile.readlines() if x.strip().split()[-1]!='description' and x.strip().split()[0]=='SCORE:']

        if not silent:
            try:
                score_ix = raw_rescore[0].index('total_score')
                name_ix = raw_rescore[0].index('description')
                headend = 0
            except ValueError:
                score_ix = raw_rescore[1].index('total_score')
                name_ix = raw_rescore[1].index('description')            
                headend = 1
        else:
            score_ix = 1
            name_ix = -1
            headend = -1

        for line in raw_rescore[headend+1:]:
            if silent:
                curmodel = line[name_ix]
            else:
                curmodel = "%s.pdb" % line[name_ix][:-5]
            curscore = float(line[score_ix])
            curkey = None
            if curmodel[:11] == 'repacked_wt':
                curkey = 'wt'
            else:
                for mut in muts:
                    if curmodel.startswith("mut_%s" % mut):
                        curkey = mut
            if curkey is None:
                continue
            if curkey not in top_three:
                top_three[curkey] = [[curscore,curmodel]]
            elif len(top_three[curkey])<3:
                top_three[curkey].append([curscore,curmodel])
            elif curscore<sorted(top_three[curkey])[-1][0]:
                print sorted(top_three[curkey])[-1][0]
                top_three[curkey][-1] == [curscore,curmodel]
        allmeans = {}
        for keys in top_three:
            topmod = sorted(top_three[keys])[0][1]
            if silent:
                if keys == 'wt':
                    outfile = "wt_.out"
                else:
                    outfile = "mut_%s.out" % keys
                command = "extract_pdbs.linuxgccrelease -in:file:silent %s -in:file:tags %s" % (outfile, topmod)
                topmod+= ".pdb"
                extcom = subprocess.Popen("%s/%s" % (self._rosetta, command), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                (cstout,csterr) = extcom.communicate()                            
            copyfile(topmod, 'top_models/%s' % topmod)
            curmean = mean([x[0] for x in top_three[keys]])
            allmeans[keys] = curmean
        scores = []            
        for keys in allmeans:
            if keys=='wt': continue
            scores.append("%s\t%.3f" % (keys, (allmeans[keys]-allmeans['wt'])))
        return [True,scores]                        
                
    def run(self, iterations=50, standard_ddg = True,squality=True,silent=True):
        commandfile = "commands.txt"
        #Always use silent file/standard ddg for low quality models
        if not squality:
            silent = True
            standard_ddg = True            
        for x in self._applications:
            try:
                assert os.path.isfile("%s/%s" % (self._rosetta,x))
            except AssertionError:
                return (False, "Required Rosetta application %s not found in path %s" % (x,self._rosetta))
        with open(commandfile,'a') as outfile:
            outfile.write("silent: %s\nstandard_ddg: %s\nRow16: %s\n" % (silent,standard_ddg,squality))
        mutationtext = ",".join("".join(str(y) for y in x) for x in self._mutations)
        suffix = datetime.datetime.now().strftime("%d%m%y%H%M%S")       
        isoresult, isopdbname, to_pdb = isolate_chain(self._pdb, self._chain)
        if not isoresult:
            return self.failure(None,to_pdb)
        else:
            self._to_pdb = to_pdb
        #Convert mutation to pose numbering for internal consistency
        mutation_pose = []
        pose_to_raw = {}

        for mut in self._mutations:
            rawmut = "".join(str(y) for y in mut)
            mutation_pose.append(mut[:])
            mutation_pose[-1][1] = to_pdb.index(mut[1])
            posemut = "".join(str(y) for y in mutation_pose[-1])
            pose_to_raw[posemut] = rawmut
            print "Mutation residue has been converted from %d to pose #%d (internally)." % (mut[1],mutation_pose[-1][1])
            with open(commandfile,'a') as outfile:
                outfile.write("Mutation residue has been converted from %d to pose #%d (internally).\n" % (mut[1],mutation_pose[-1][1]))
        print self._break
  
        ##Preminimize and gather restraints
        #Premin only takes a list of structs even if you only have 1
        templistfile = "temp%s.ls" % suffix
        with open(templistfile, 'w') as outfile:
            outfile.write(isopdbname)
            outfile.write("\n")
        command = """minimize_with_cst.default.linuxgccrelease -in:file:fullatom -fa_max_dis 9.0 -ddg:harmonic_ca_tether 0.5 -ddg::constraint_weight 1.0 -ddg::out_pdb_prefix minimized -ddg::sc_min_only false -score:weights talaris2014 -in:file:l %s -overwrite -ignore_zero_occupancy false""" % (templistfile)
        #Add to minimize protocol    -in:file:l %s 
        #Add to ddg -ddg:iterations %d\n" % iterations
        with open(commandfile,'a') as outfile:
            outfile.write("---\n")
            outfile.write("Preminimize\n")
            outfile.write("%s/%s" % (self._rosetta,command))
        minstart = default_timer()
        print "Generating starting constraints with command:"
        print "%s/%s" % (self._rosetta, command)
        runcst = subprocess.Popen("%s/%s" % (self._rosetta, command), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (cstout,csterr) = runcst.communicate()

        csts = []
        for line in cstout.split("\n"):
            if len(line)>7 and line[:7]=="c-alpha":
                cur = line.split()
                csts.append("AtomPair CA %s CA %s HARMONIC %s %s" % (cur[5],cur[7],cur[9],cur[12]))
        try:
            assert len(csts)>0
        except AssertionError:
            cstoutfile = "pmfail%s.out" % suffix
            csterrfile = "pmfail%s.err" % suffix
            return self.failure([["relax output",cstoutfile,cstout],["relax stderr",csterrfile,csterr]],"Error running preminimization generation of csts (cst list is empty).")

        minstop = default_timer()
        mintime = minstop-minstart
        os.remove(templistfile)
        minimized_pdb = "minimized.%s_0001.pdb" % (os.path.basename(isopdbname).split(".")[0])
        print minimized_pdb
        try:
            assert os.path.isfile(minimized_pdb)
        except AssertionError:
            cstoutfile = "pmfail%s.out" % suffix
            csterrfile = "pmfail%s.err" % suffix
            return self.failure([["relax output",cstoutfile,cstout],["relax stderr",csterrfile,csterr]],"Error running preminimization generation of csts (minimized pdb not found).")  

        ##Rescore model and get residue score
        command = """per_residue_energies.linuxgccrelease -s %s -score:weights talaris2014 \
         -out:file:silent min_residues_%s.sc""" % (minimized_pdb, suffix)
        print self._break
        print "Rescoring minimized model with command:"
        print "%s/%s" % (self._rosetta,command)
        with open(commandfile,'a') as outfile:
            outfile.write("---\nRescore:\n")
            outfile.write("%s/%s" % (self._rosetta,command))
        scorestart = default_timer()
        runscore = subprocess.Popen("%s/%s" % (self._rosetta, command), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (scoreout,scoreerr) = runscore.communicate()
        with open("min_residues_%s.sc" % suffix) as infile:
            raw_scores = [x.strip().split()[-2:] for x in infile.readlines() if len(x.strip().split())>0 and x.strip().split()[-1] != "description"]
        try:
            assert len(raw_scores) == len(to_pdb)-1
        except AssertionError:
            scoreoutfile = "scorefail%s.out" % suffix
            scoreerrfile = "scorefail%s.err" % suffix
            return self.failure([["score stdout",scoreoutfile,scoreout],["score stderr",scoreerrfile,scoreerr]],\
        "Error rescoring minimized model, scored residues %d != expected %d" % (len(raw_scores),len(to_pdb)-1))
        scorestop = default_timer()
        scoretime = scorestop-scorestart
        mutation_score = {}
        for x in self._mutations:
            z = "".join(str(y) for y in x)
            mutation_score[z] = None
        for x in xrange(len(raw_scores)):
            resnum = to_pdb[x+1]
            if self._mutations is not None:
                for z in self._mutations:
                    if resnum == z[1]:
                        zz = "".join(str(y) for y in z)
                        mutation_score[zz] = raw_scores[x][0]
        for ms in mutation_score:
            try:
                assert mutation_score[ms] is not None
            except AssertionError:
                scoreoutfile = "scorefail%s.out" % suffix
                scoreerrfile = "scorefail%s.err" % suffix
                failure([["score stdout",screoutfile,scoreout],["score stderr",scoreerrfile,scoreerr]],\
                    "Mutation residue was not found in score file. Something went wrong.")

        ##Set up ddg_monomer
        #Generate cst file
        cstfile = "%s_%s_%s.cst" % (self._pdb, self._chain, mutationtext)
        with open(cstfile,'w') as outfile:
            outfile.write("\n".join(csts))
            outfile.write("\n")

        #Generate mut and options files
        mutfile = "%s_%s_%s.mut" % (self._pdb,self._chain, mutationtext)
        with open(mutfile,'w') as outfile:
            outfile.write("total %d\n" %(len(mutation_pose)+1 if len(mutation_pose)>1 else 1))
            for mut in mutation_pose:
                outfile.write("1\n")
                outfile.write(" ".join(str(x) for x in mut))
                outfile.write("\n")
            if len(self._mutations)>1:
                outfile.write("%d\n" % len(mutation_pose))
                for mut in mutation_pose:
                    outfile.write(" ".join(str(x) for x in mut))
                    outfile.write("\n")

        optionsfile = "%s_%s_%s.options" % (self._pdb,self._chain,mutationtext)
        if squality:
            protocol = "Row16"
            #Row 16 options
            with open(optionsfile, 'w') as outfile:
                outfile.write("#Based on kellogg et al 2011 table 1 row 16\n\n")
                outfile.write("-in:file:fullatom\n-ddg:weight_file soft_rep_design\n")
                outfile.write("-fa_max_dis 9.0\n-ddg:iterations %d\n" % iterations)
                outfile.write("-ddg:dump_pdbs true\n-ddg:local_opt_only false\n-ddg:min_cst true\n")
                outfile.write("-ddg:suppress_checkpointing true\n-ddg:mean false\n-ddg:min true\n-ddg:sc_min_only false\n")
                outfile.write("-ddg:ramp_repulsive true\n-ddg:minimization_scorefunction talaris2014\n")
                outfile.write("-ignore_zero_occupancy false\n")
                outfile.write("-ddg:output_silent true\n")
        else:
            protocol = "Row3"
            #Row 3 options
            with open(optionsfile, 'w') as outfile:
                outfile.write("#Based on kellogg et al 2011 table 1 row 3\n\n")
                outfile.write("-in:file:fullatom\n-ddg:weight_file soft_rep_design\n")
                outfile.write("-fa_max_dis 9.0\n-ddg:iterations %d\n" % iterations)
                outfile.write("-ddg:dump_pdbs true\n-ddg:local_opt_only true\n")
                outfile.write("-ddg:suppress_checkpointing true\n-ddg:mean true\n-ddg:min false\n")
                outfile.write("-ddg:minimization_scorefunction talaris2014\n")
                outfile.write("-ignore_zero_occupancy false\n")            
                outfile.write("-ddg:output_silent true\n")

        ##Run actual ddg_monomer application
        if os.path.isfile("ddg_predictions.out"):
            shutil.move("ddg_predictions.out","ddg_predictions.old_%s" % suffix)
        command = "ddg_monomer.linuxgccrelease @%s -constraints::cst_file %s -s %s -ddg:mut_file %s" % (optionsfile,cstfile,minimized_pdb,mutfile)
        print self._break
        if squality:
            print "Running high quality (row 16) ddG monomer protocol"
        else:
            print "Running low quality (row 3) ddG monomer protocol"
        print "Running ddg_monomer application with %d iterations (this may take a while)" % iterations
        print "Command:"
        print "%s/%s" % (self._rosetta,command)
        with open(commandfile,'a') as outfile:
            outfile.write("---\nddG:\n")
            outfile.write("%s/%s" % (self._rosetta,command))
            outfile.write("\n")
        ddgstart = default_timer()
        runddg = subprocess.Popen("%s/%s" % (self._rosetta, command), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (ddgout,ddgerr) = runddg.communicate()
        ddgstop = default_timer()
        ddgtime = ddgstop-ddgstart
        
        #Make sure all expected files have been generated
        if not silent:
            expected_files = ["repacked_wt_round_%d.pdb" % (x+1) for x in xrange(iterations)]+["ddg_predictions.out"]
            for mut in mutation_pose:
                expected_files += ["mut_%s_round_%d.pdb" % ("".join(str(x) for x in mut),x+1) for x in xrange(iterations)]
        else:
            expected_files = ["ddg_predictions.out"]
        for x in expected_files:
            try:
                assert os.path.isfile(x)
            except AssertionError:
                ddgoutfile = "ddgfail%s.out" % suffix
                ddgerrfile = "ddgfail%s.err" %suffix
                return self.failure([["ddg stdout",ddgoutfile,ddgout],["ddg stderr",ddgerrfile,ddgerr]],\
                "Something went wrong with ddg_monomer, expected file %s is missing!" % x)

        #These files are always empty
        if os.path.isfile("wt_traj"):
            os.remove("wt_traj")
        for mut in mutation_pose:
            if os.path.isfile("mutant_traj%s" % "".join(str(x) for x in mut)):
                os.remove("mutant_traj%s" % "".join(str(x) for x in mut))
    
        with open("ddg_predictions.out") as infile:
            final_ddg = [x.strip().split() for x in infile.readlines() if len(x.strip().split())>2 and x.strip().split()[1]!='description']

        expected = len(self._mutations)+1 if len(self._mutations)>1 else 1
        try:
            assert len(final_ddg)==expected
        except AssertionError:
            ddgoutfile = "ddgfail%s.out" % suffix
            ddgerrfile = "ddgfail%s.err" % suffix
            with open(ddgoutfile, 'w') as outfile:
                outfile.write(ddgout)
            with open(ddgerrfile,'w') as outfile:
                outfile.write(ddgerr)    
            return self.failure([["ddg stdout",ddgoutfile,ddgout],["ddg stderr",ddgerrfile,ddgerr]],\
            "Something went wrong with ddg_monomer, expected %d prediction line, got %d" % (expected,len(final_ddg)))

        mutation_ddg = {}
        temp = mutation_pose[:]
        temp.sort(key=lambda x:x[1])
        tempraw = "".join("".join(str(x) for x in y) for y in sorted(self._mutations,key=lambda z: z[1]))
        #Standard ddg uses the output ddG score as provided by application
        #Nonstandard ddg uses the difference in pose scores between the mean score of the top 3 models per condition.       
        if standard_ddg:
            for line in final_ddg:
                mut = line[1]
                dg = line[2]
                try:
                    curmut = pose_to_raw[mut]
                except KeyError:
                    curmut = tempraw               
                mutation_ddg[curmut] = dg         
        else:
            scout, scres = self.get_ddg(squality,silent)      
            if not scout:
                return [False,scres]
            for line in scres:           
                mut,dg = line.split("\t")
                try:
                    curmut = pose_to_raw[mut]
                except KeyError:
                    curmut = tempraw
                mutation_ddg[curmut] = dg
        results = ["%s\t%s\t%s\t%s\t%s\t%s\t%s" % ("File","Chain","Residue","WT_Res_Score","Mutation","ddG","Protocol")]

        for x in sorted(self._mutations,key=lambda zz: zz[1]):
            xstr = "".join(str(y) for y in x)
            results += ["%s\t%s\t%s\t%s\t%s\t%s\t%s" % (os.path.basename(self._pdb),self._chain,str(x[1]),mutation_score[xstr],xstr,mutation_ddg[xstr],protocol)]

        if len(self._mutations)>1:
            results += ["%s\t%s\t%s\t%s\t%s\t%s\t%s" % (os.path.basename(self._pdb),self._chain,"All","All","All",mutation_ddg[tempraw],protocol)]          
            #self.log_success({'minimize': mintime,'rescore': scoretime,'ddg_monomer':ddgtime})
        return [True, results]                                                  
        
class interface_analyzer(object):
    def __init__(self,pdb,chain,mutations,rosetta_path='/dors/meilerlab/apps/rosetta/rosetta-3.7/main/source/bin'):
        self._pdb = pdb
        self._chain = chain
        self._mutations = mutations
        self._rosetta = rosetta_path.rstrip('/')
        self._break = "".join('*'*20)
        self._success = "/dors/meilerlab/home/sliwosgr/udn_pipeline/interface_energy_calculator/history/success_log.tab"
        self._fail = "/dors/meilerlab/home/sliwosgr/udn_pipeline/interface_energy_calculator/history/fail_log.tab"
        self._timers = ['minimize','interface_analyzer']
        self._applications = ['minimize_with_cst.default.linuxgccrelease','residue_energy_breakdown.default.linuxgccrelease']
        self._index = ['resi1','resi2','total']
        self._topdb = []
        LOGGER.info("interface_analyzer.__init__ completed after setting:\n%s"%pprint.pformat(self.__dict__))

    def log_success(self,times):
        mutationtext = ",".join("".join(str(x) for x in y) for y in self._mutations)
        nres = len(self._topdb)-1
        curdate = datetime.datetime.now().strftime("%d_%m_%y")       
        curtime = datetime.datetime.now().strftime("%H_%M")
        wd = os.getcwd()  
        if not os.path.isfile(self._success):
            with open(self._success,'w') as outfile:
                header = ['date','time','pdbfile','chain','nres','residue']+self._timers+['path']
                outfile.write("\t".join(header))
                outfile.write("\n")
            with open(self._success,'a') as outfile:
                line = [curdate,curtime,self._pdb,self._chain,nres,mutationtext]+\
                ["%.2f" % (times[x]/60.0) if x in times else 'NaN' for x in self._timers]+[wd]
                outfile.write("\t".join(line))
                outfile.write("\n")

    def log_failure(self,message):
        wd = os.getcwd()
        mutationtext = ",".join("".join(str(x) for x in y) for y in self._mutations)
        curdate = datetime.datetime.now().strftime("%d_%m_%y")       
        curtime = datetime.datetime.now().strftime("%H_%M")
        if not os.path.isfile(self._fail):
            with open(self._fail,'w') as outfile:
                header = ['date','time','pdbfile','chain','residue','message','path']
                outfile.write("\t".join(header))
                outfile.write("\n")
        with open(self._fail,'a') as outfile:
            line = [curdate,curtime,self._pdb,self._chain,mutationtext,message,wd]
            outfile.write("\t".join(line))
            outfile.write("\n")
            
    def failure(self, dumps, message):
        print message
        if dumps is not None:
            for item in dumps:
                print "Dumping %s to %s" % (item[0],item[1])
                with open(item[1],'w') as outfile:
                    outfile.write(item[2])
        #self.log_failure(message)
        return [False,message]
           
    def run(self):
        for x in self._applications:
            try:
                assert os.path.isfile("%s/%s" % (self._rosetta,x))
            except AssertionError:
                return [False, "Error: required Rosetta application %s not found in provided path" % x]
        mutationtext = ",".join("".join(str(y) for y in x) for x in self._mutations)
        suffix = datetime.datetime.now().strftime("%d%m%y%H%M%S")       

        #Renum pdb for consistency and check to make sure chain is present + at least 1 other chain
        #If failed chain checks, still write results to files for simple post processing
        renumresult, renumpdbname, to_pdb, allchains = renum_pdb(self._pdb, self._chain,2)      
        if not renumresult:
            print "Unable to calculate interface energy: %s" % to_pdb
            results = ["\t".join(["File","Chain","Residue","Interface_Chain","Total_Score","Abs_Total_Score"])]
            for x in self._mutations:
                results += ["\t".join([self._pdb,self._chain,str(x[1]),'NA','NA','NA'])]
            return [True,results]
        self._to_pdb = to_pdb
        residue_pose = []

        for mut in self._mutations:
            res = self._chain+str(mut[1])
            try:
                residue_pose.append(self._to_pdb.index(res))
            except ValueError:
                print "%s not found in biounit, skipping" % res
                residue_pose.append(None)
                continue
            print "Selected residue number has been converted from %s to pose #%d (internally)." % (res[1:],residue_pose[-1])
        print self._break

        ##Preminimize to equilibriate to Rosetta energy environment
        #Premin only takes a list of structs even if you only have 1
        templistfile = "temp%s.ls" % suffix
        with open(templistfile, 'w') as outfile:
            outfile.write(renumpdbname)
            outfile.write("\n")

        pmlogfile = "%s_premin.log" % os.path.splitext(renumpdbname)[0]
        command = "minimize_with_cst.default.linuxgccrelease -in:file:fullatom -fa_max_dis 9.0 -ddg:harmonic_ca_tether 0.5 -ddg::constraint_weight 1.0 -ddg::out_pdb_prefix  minimized -ddg::sc_min_only false -score:weights talaris2014 -in:file:l %s -overwrite -ignore_zero_occupancy false" % templistfile
        minstart = default_timer()
        print "Equilibriating structure to Rosetta energy function with command:"
        print "%s/%s" % (self._rosetta, command)
        runmin = subprocess.Popen("%s/%s" % (self._rosetta, command), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (minout,minerr) = runmin.communicate()
        minstop = default_timer()
        mintime = minstop-minstart
        os.remove(templistfile)

        minimized_pdb = "minimized.%s_0001.pdb" % os.path.basename(os.path.splitext(renumpdbname)[0])
        try:
            assert os.path.isfile(minimized_pdb)
        except AssertionError:
            minoutfile = "pmfail%s.out" % suffix
            minerrfile = "pmfail%s.err" % suffix
            return(self.failure([["relax output",minoutfile,minout],["relax stderr",minerrfile,minerr]],"Error running preminimization."))

        ##Generate residue-pair energy score raw file
        command = "residue_energy_breakdown.default.linuxgccrelease -s %s -score:weights talaris2014 -out:file:silent residue_pairs_%s.raw -overwrite" % (minimized_pdb,suffix)
        print self._break
        print "Generating residue-pair energies from minimized model with command:"
        print "%s/%s" % (self._rosetta,command)
        scorestart = default_timer()
        runscore = subprocess.Popen("%s/%s" % (self._rosetta, command), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (scoreout,scoreerr) = runscore.communicate()
        rawscorefile = "residue_pairs_%s.raw" % suffix
        try:
            assert os.path.isfile(rawscorefile)
        except AssertionError:
            scoreoutfile = "scorefail%s.out" % suffix
            scoreerrfile = "scorefail%s.err" % suffix
            return(self.failure([["score stdout",scoreoutfile,scoreout],["score stderr",scoreerrfile,scoreerr]],\
            "Residue-pair energy extraction failed."))
   
        with open(rawscorefile) as infile:
            raw = [x.strip().split() for x in infile.readlines()]
            try:
                r1ix = raw[0].index(self._index[0])
                r2ix = raw[0].index(self._index[1])
                tix = raw[0].index(self._index[2])
            except ValueError:
                scoreoutfile = "scorefail%s.out" % s
                scoreerrfile = "scorefail%s.err" % s
                return(self.failure([["score stdout",scoreoutfile,scoreout],["score stderr",scoreerrfile,scoreerr]],\
                "Unable to locate required header entries in raw energy file %s, may need to update constants in script." % rawscorefile))
            #Note: disregarding onebody scores. No good reason so may need to change in future.
            raw_scores = [[int(x[r1ix]),int(x[r2ix]),float(x[tix])] for x in raw[1:] if x[r2ix]!='--']

        processed_scores = {}
        processed_scores_abs = {}   
        for item in raw_scores:
            c1, r1 = to_pdb[item[0]][0],to_pdb[item[0]][1:]
            c2, r2 = to_pdb[item[1]][0],to_pdb[item[1]][1:]
            curscore = item[2]
            if c1==c2: continue
            if c1 not in processed_scores:
                processed_scores[c1] = {}
            if r1 not in processed_scores[c1]:
                processed_scores[c1][r1] = {}
            if c2 not in processed_scores[c1][r1]:
                processed_scores[c1][r1][c2] = curscore
            else:
                processed_scores[c1][r1][c2] += curscore
            if c2 not in processed_scores:
                processed_scores[c2] = {}
            if r2 not in processed_scores[c2]:
                processed_scores[c2][r2] = {}
            if c1 not in processed_scores[c2][r2]:
                processed_scores[c2][r2][c1] = curscore
            else:
                processed_scores[c2][r2][c1] += curscore
            if c1 not in processed_scores_abs:
                processed_scores_abs[c1] = {}
            if r1 not in processed_scores_abs[c1]:
                processed_scores_abs[c1][r1] = {}
            if c2 not in processed_scores_abs[c1][r1]:
                processed_scores_abs[c1][r1][c2] = abs(curscore)
            else:
                processed_scores_abs[c1][r1][c2] += abs(curscore)
            if c2 not in processed_scores_abs:
                processed_scores_abs[c2] = {}
            if r2 not in processed_scores_abs[c2]:
                processed_scores_abs[c2][r2] = {}
            if c1 not in processed_scores_abs[c2][r2]:
                processed_scores_abs[c2][r2][c1] = abs(curscore)
            else:
                processed_scores_abs[c2][r2][c1] += abs(curscore)
        scorestop = default_timer()
        scoretime = scorestop-scorestart

        results = ["\t".join(["File","Chain","Residue","Interface_Chain","Total_Score","Abs_Total_Score"])]
        for x in sorted(allchains):
            if x==self._chain: continue
            for y in xrange(len(self._mutations)):
                if residue_pose[y] is None:
                    cursc = 0.0
                    curabs = 0.0
                else:
                    try:
                        cursc = processed_scores[self._chain][str(self._mutations[y][1])][x]
                        curabs = processed_scores_abs[self._chain][str(self._mutations[y][1])][x]
                    except KeyError:
                        cursc = 0.0
                        curabs = 0.0
                results += ["\t".join([self._pdb,self._chain,str(self._mutations[y][1]),x,"%.3f" % cursc,"%.3f" % curabs])]
        #self.log_success({'minimize': mintime,'interface_analyzer': scoretime})        
        return [True, results]
        
class dssp(object):
    def __init__(self,pdb,chain,mutations,dssp_path='/dors/meilerlab/home/sliwosgr/udn_pipeline/dssp/dssp_local.exe'):
        self._pdb = pdb
        self._chain = chain
        self._mutations = mutations
        self._dssp = dssp_path.rstrip('/')
        self._break = "".join('*'*20)
        self._TOTAL_SSA = {"A": 113, "R":241, "N":158, "D":151, "C":140, "Q":189,\
                           "E":183, "G":85, "H":194, "I":182, "L":180, "K":211, "M":204, "F":218,\
                           "P":143, "S":122, "T":146, "W":259, "Y":229, "V":160}
        self._success = "/dors/meilerlab/home/sliwosgr/udn_pipeline/dssp/history/success_log.tab"
        self._fail = "/dors/meilerlab/home/sliwosgr/udn_pipeline/dssp/history/fail_log.tab"
        LOGGER.info("dssp.__init__  completed after setting:\n%s"%pprint.pformat(self.__dict__))
        
    def log_success(self):
        mutationtext = ",".join("".join(str(x) for x in y) for y in self._mutations)
        curdate = datetime.datetime.now().strftime("%d_%m_%y")       
        curtime = datetime.datetime.now().strftime("%H_%M")
        wd = os.getcwd()  
        if not os.path.isfile(self._success):
            with open(self._success,'w') as outfile:
                header = ['date','time','pdbfile','chain','residue','path']
                outfile.write("\t".join(header))
                outfile.write("\n")
            with open(self._success,'a') as outfile:
                line = [curdate,curtime,self._pdb,self._chain,mutationtext,wd]
                outfile.write("\t".join(line))
                outfile.write("\n")

    def log_failure(self,message):
        wd = os.getcwd()
        mutationtext = ",".join("".join(str(x) for x in y) for y in self._mutations)
        curdate = datetime.datetime.now().strftime("%d_%m_%y")       
        curtime = datetime.datetime.now().strftime("%H_%M")
        if not os.path.isfile(self._fail):
            with open(self._fail,'w') as outfile:
                header = ['date','time','pdbfile','chain','residue','message','path']
                outfile.write("\t".join(header))
                outfile.write("\n")
        with open(self._fail,'a') as outfile:
            line = [curdate,curtime,self._pdb,self._chain,mutationtext,message,wd]
            outfile.write("\t".join(line))
            outfile.write("\n")
            
    def failure(self, dumps, message):
        print message
        if dumps is not None:
            for item in dumps:
                print "Dumping %s to %s" % (item[0],item[1])
                with open(item[1],'w') as outfile:
                    outfile.write(item[2])
        #self.log_failure(message)
        return [False,message]
        
    def run(self):
        results = []
        centroids = {}
        atoms = []
        chainpdbs = {}
        rawisolated = {}
        isolated_sas = {}       
        with open(self._pdb) as infile:
            pdblines = infile.readlines()
        mutationtext = ",".join("".join(str(y) for y in x) for x in self._mutations)
        suffix = datetime.datetime.now().strftime("%d%m%y%H%M%S")
        try:
            assert os.path.isfile(self._dssp)
        except AssertionError:
            return(self.failure(False, "Error: %s not found" % self._dssp))
        
        command = "%s -i %s" % (self._dssp,self._pdb)
        print "Running DSSP for whole biounit with command:"
        print command        
        rundssp = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (out,err) = rundssp.communicate()        
        raw = out.split("\n")
        try:
            assert len(raw)>0
        except AssertionError:
            return(self.failure(None,"Failed to run dssp on whole biomodel"))

        prev = None
        for line in pdblines:
            if line[:4] != "ATOM": continue
            res = int(line[22:26])
            chain = line[21]
            runner = chain+str(res)
            coord = [float(line[30:38]),float(line[38:46]),float(line[46:54])]
            if not prev or runner==prev:
                atoms.append(coord)
                prev = runner
            else:
                avgc = avg(atoms)
                prev = runner
                atoms = [coord]
                centroids[runner] = avgc
            if len(atoms) > 0:
                avgc = avg(atoms)
                centroids[runner] = avgc
            if chain not in chainpdbs:
                chainpdbs[chain] = [line]
            else:
                chainpdbs[chain].append(line)
        print self._break
        
        for chain in chainpdbs:
            if chain!=self._chain: continue
            tempfile = 'temp_%s_%s.temp.' % (chain, suffix)
            with open(tempfile, 'w') as outtemp:
                outtemp.write("".join(chainpdbs[chain]))  
            command = "%s -i %s" % (self._dssp,tempfile)
            print "Running DSSP on isolated chain %s with command:" % self._chain
            print command
            rundssp = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (out,err) = rundssp.communicate()
            if len(err)>0:
                return(self.failure(None, "dssp failed with error: %s" % " ".join(err.split("\n"))))
            rawisolated[chain] = out.split("\n")
            if len(rawisolated[chain])==0:
                return(self.failure(None,"Failed DSSP with isolated chain %s" % chain))
            os.remove('temp_%s_%s.temp.' % (chain, suffix))
        print self._break
        
        for chain in chainpdbs:
            if chain != self._chain: continue
            for line in rawisolated[chain]:
                if len(line)==0: continue
                try:
                    res = int(line[5:11])
                except ValueError:
                    continue
                try:
                    sas = int(line[35:38])
                except ValueError:
                    continue
                chainds = line[11]
                aa = line[13]
                #DSSP assigns lowercase partners to bridged cysteins
                if aa.islower():
                    aa = 'C'
                if chainds==" " or aa==" " or chainds!=chain: continue
                rsa = float(sas)/self._TOTAL_SSA[aa]
                if rsa>1:
                    rsa = 1.0
                dist_surf = 0 if rsa>=.1 else None
                isolated_sas[chain+str(res)] = [sas,rsa,dist_surf]

        for line in raw:
            if len(line)==0: continue
            try:
                res = int(line[5:11])
            except ValueError:
                continue
            try:
                sas = int(line[35:38])
            except ValueError:
                continue
            aa = line[13]
            #DSSP assigns lowercase partners to bridged cysteins
            if aa.islower():
                aa = 'C'
            chain = line[11]
            if chain==" " or aa==" " or chain!=self._chain: continue
            ss = "X" if line[16]==" " else line[16]
            rsa = float(sas)/self._TOTAL_SSA[aa]
            if rsa>1:
                rsa = 1.0
            dist_surf = 0 if rsa>=.1 else None
            centroid = centroids[chain+str(res)]
            iso_sas,iso_rsa,iso_ds = isolated_sas[chain+str(res)]
            dsas = iso_sas-sas
            dras = iso_rsa-rsa
            if dist_surf is None and iso_ds == 0:
                interface = "Y"
            else:
                interface= "N"
            closest = 'SELF'
            isoclosest = 'SELF'
            results.append([str(res),chain,aa,ss,str(sas),"%.2f" % rsa,dist_surf,str(iso_sas),"%.2f" % iso_rsa,iso_ds,str(dsas),"%.2f" % dras,closest,isoclosest,centroid])

        ixsas = 4
        ixrsa = 5
        ixds = 6
        ixisosas = 7
        ixisorsa = 8
        ixisods = 9            

        for line in results:
            #If both isolated and non-isolated are surface residues, continue
            if line[ixds] == 0 and line[ixisods] == 0: continue
            mindist = line[ixds]
            isomindist = line[ixisods]
            cur = line[-1]
            for other in results:
                if other[ixds] != 0 and other[ixisods] != 0: continue  
                c1 = array((float(other[-1][0]),float(other[-1][1]),float(other[-1][2])))
                c2 = array((float(cur[0]),float(cur[1]),float(cur[2])))
                dist = linalg.norm(c2-c1)
                if other[ixds] == 0:
                    if mindist is None or dist < mindist:
                        mindist = dist
                        line[-3] = other[1]+str(other[0])
                if other[ixisods] == 0 and line[1]==other[1]:
                    if isomindist is None or dist < isomindist:
                        isomindist = dist
                        line[-2] = other[1]+str(other[0])
            line[ixds] = mindist
            line[ixisods] = isomindist

        header =  "\t".join(["structure","resnum","chain","aa","SSE","SAS","RAS","dist_surface_res","SAS_isolated","RAS_isolated","dist_surface_res_isolated","deltaSAS","deltaRAS","clos_surf_whole","clos_surf_iso"])
        final_results = [header]
        seen = []
        for line in results:
            for mut in self._mutations:
                if line[0] == str(mut[1]) and line[1]==self._chain:
                    line[ixds] = "%.2f" % line[ixds]
                    line[ixisods] = "%.2f" % line[ixisods]
                    final_results += ["\t".join([self._pdb] + line[:-1])]
                    seen.append(mut)
        for mut in self._mutations:
            if mut not in seen:
                final_results += ["\t".join([self._pdb, str(mut[1]),self._chain] + ['NA']*10)]
        #self.log_success()
        return [True,final_results]
        
class ligands(object):
    def __init__(self,pdb,chain,mutations):
        self._pdb = pdb
        self._chain = chain
        self._mutations = mutations
        self._break = "".join('*'*20)
        self._SKIPPED = ['HOH']
        self._success = '/dors/meilerlab/home/sliwosgr/udn_pipeline/ligand_distances/history/success_log.tab'
        self._fail = '/dors/meilerlab/home/sliwosgr/udn_pipeline/ligand_distances/history/fail_log.tab'
        self._timers = ['total_time']
        self.residues = {}
        self.ligands = {}
        LOGGER.info("ligands.__init__  completed after setting:\n%s"%pprint.pformat(self.__dict__))
        
    def log_success(self,times):
        mutationtext = ",".join("".join(str(x) for x in y) for y in self._mutations)
        curdate = datetime.datetime.now().strftime("%d_%m_%y")       
        curtime = datetime.datetime.now().strftime("%H_%M")
        wd = os.getcwd()  
        if not os.path.isfile(self._success):
            with open(self._success,'w') as outfile:
                header = ['date','time','pdbfile','chain','nres','nlig','residue']+self._timers+['path']
                outfile.write("\t".join(header))
                outfile.write("\n")
            with open(self._success,'a') as outfile:
                line = [curdate,curtime,self._pdb,self._chain,len(self.residues),len(self.ligands),mutationtext]+\
                ["%.2f" % (times[x]/60.0) if x in times else 'NaN' for x in self._timers]
                outfile.write("\t".join(line))
                outfile.write("\n")

    def log_failure(self,message):
        wd = os.getcwd()
        mutationtext = ",".join("".join(str(x) for x in y) for y in self._mutations)
        curdate = datetime.datetime.now().strftime("%d_%m_%y")       
        curtime = datetime.datetime.now().strftime("%H_%M")
        if not os.path.isfile(self._fail):
            with open(self._fail,'w') as outfile:
                header = ['date','time','pdbfile','chain','residue','message','path']
                outfile.write("\t".join(header))
                outfile.write("\n")
        with open(self._fail,'a') as outfile:
            line = [curdate,curtime,self._pdb,self._chain,mutationtext,message,wd]
            outfile.write("\t".join(line))
            outfile.write("\n")
            
    def failure(self, dumps, message):
        print message
        if dumps is not None:
            for item in dumps:
                print "Dumping %s to %s" % (item[0],item[1])
                with open(item[1],'w') as outfile:
                    outfile.write(item[2])
        #self.log_failure(message)
        return [False,message]              

    def run(self):
        times = {}
        starttime = default_timer() 
        res_distances = {}      
        respos = [x[1] for x in self._mutations]
        mutationtext = ",".join("".join(str(y) for y in x) for x in self._mutations)
        suffix = datetime.datetime.now().strftime("%d%m%y%H%M%S")
        header = "\t".join(['filename','chain','residue','ligand','distance'])
        result = [header]
        if self._pdb[-3:] == '.gz':
            with gzip.open(self._pdb) as infile:
                pdblines = infile.readlines()
        else:
            with open(self._pdb) as infile:
                pdblines = infile.readlines()
                
        for line in pdblines:
            if len(line)<54 or (line[:4]!='ATOM' and line[:6]!='HETATM'): continue
            if len(line)>6 and line[:6]=='ENDMDL':break #Only take first model of NMR structures
            curchain = line[21]
            #Skip alternate positions and residues with insertion codes
            if line[16] !='A' and line[16] != ' ': continue
            label = line[17:20].strip()
            if line[26] !=' ' or (curchain!=self._chain and label not in ['DT','DC','DG','DA']): continue
            resnum = int(line[22:26])
            if line[:4]=='ATOM' and (resnum not in respos and label not in ['DT','DC','DG','DA']): continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])  

            #Skip anything in skip list
            if label in self._SKIPPED: continue
            #Treat modified residues as normal residues
            #Otherwise, store hetatom as ligand
            if line[:6]=='HETATM' and query_modres(line[17:20]) is None and line[17:20]!='MSE':
                full_label = "_".join([str(resnum),label])
                curatom = [resnum,full_label,x,y,z]
                if full_label not in self.ligands:
                    self.ligands[full_label] = [curatom]
                else:
                    self.ligands[full_label].append(curatom)
            #DNA as a ligand
            elif label in ['DT','DC','DG','DA']:
                full_label = 'DNA'
                curatom = [resnum,full_label,x,y,z]
                if full_label not in self.ligands:
                    self.ligands[full_label] = [curatom]
                else:
                    self.ligands[full_label].append(curatom)                                
            else:
                curatom = [resnum,label,x,y,z]
                if resnum not in self.residues:
                    self.residues[resnum] = [curatom]
                else:
                    self.residues[resnum].append(curatom)                  

            
        try:
            assert len(self.ligands)>0
        except AssertionError:
            print "No ligands in %s" % self._pdb
            for x in respos:
                result += ["\t".join([self._pdb,self._chain,str(x),'NA','NA'])]
            return [True, result]
            
        for item in self.residues:
            for het in self.ligands:
                mindist = None
                for atom1 in self.residues[item]:
                    for atom2 in self.ligands[het]:
                        a = atom1[-3:]
                        b = atom2[-3:]
                        dist = distance(a,b)
                        if mindist is None or dist<mindist:
                            mindist = dist
                if item not in res_distances:
                    res_distances[item] = [[het,mindist]]
                else:
                    res_distances[item].append([het,mindist])

        for item in res_distances:
            for lig in res_distances[item]:
                result += ["\t".join([self._pdb,self._chain,str(item),lig[0],"%.3f" % lig[1]])]

        stoptime = default_timer()
        times['total_time'] = starttime-stoptime       
        #self.log_success(times)
        return [True,result]

class uniprot(object):
    def __init__(self,gene,mutations):
        self._ANNOTATIONS = {
        'FT':[
        'ACT_SITE','BINDING','CA_BIND','CARBOHYD','CHAIN','COILED','COMPBIAS','CROSSLINK','DISULFID','DNA_BIND','DOMAIN',
        'HELIX','INTRAMEM','LIPID','METAL','MOD_RES','MOTIF','MUTAGEN','NP_BIND','PEPTIDE','PROPEP',
        'REGION','REPEAT','SIGNAL','SITE','STRAND','TOPO_DOM','TRANSIT','TRANSMEM','TURN','VARIANT','ZN_FING']}
        self._PAIRS = ['DISULFID']             
        self._DESCRIPTIONS = {
        'ACT_SITE': 'Amino acid(s) directly involved in the activity of an enzyme',
        'BINDING': 'Binding site for any chemical group (co-enzyme, prosthetic group, etc.)',
        'CA_BIND': 'Position(s) of calcium binding region(s) within the protein',
        'CARBOHYD': 'Covalently attached glycan group(s)',
        'CHAIN': 'Extent of a polypeptide chain in the mature protein',
        'COILED': 'Positions of regions of coiled coil within the protein',
        'COMPBIAS': 'Region of compositional bias in the protein',
        'CROSSLINK':'Residues participating in covalent linkage(s) between proteins',
        'DISULFID':'Cysteine residues participating in disulfide bonds',
        'DNA_BIND':'Position and type of a DNA-binding domain',
        'DOMAIN':'Position and type of each modular protein domain',
        'HELIX':'Helical regions within the experimentally determined protein structure',
        'INTRAMEM':'Extent of a region located in a membrane without crossing it',
        'LIPID':'Covalently attached lipid group(s)',
        'METAL':'Binding site for a metal ion',
        'MOD_RES':'Modified residues excluding lipids, glycans and protein cross-links',
        'MOTIF':'Short (up to 20 amino acids) sequence motif of biological interest',
        'MUTAGEN':'Site which has been experimentally altered by mutagenesis',
        'NP_BIND':'Nucleotide phosphate binding region',
        'PEPTIDE':'Extent of an active peptide in the mature protein',
        'PROPEP':'Part of a protein that is cleaved during maturation or activation',
        'REGION':'Region of interest in the sequence',
        'REPEAT':'Positions of repeated sequence motifs or repeated domains',
        'SIGNAL':'Sequence targeting proteins to the secretory pathway or periplasmic space.',
        'SITE':'Any interesting single amino acid site on the sequence',
        'STRAND':'Beta strand regions within the experimentally determined protein structure',
        'TOPO_DOM':'Location of non-membrane regions of membrane-spanning proteins',
        'TRANSIT':'Extent of a transit peptide for organelle targeting',
        'TRANSMEM':'Extent of a membrane-spanning region',
        'TURN':'Turns within the experimentally determined protein structure',
        'VARIANT':'Description of a natural variant of the protein',
        'ZN_FING':'Position(s) and type(s) of zinc fingers within the protein'}
        self._UNP_MAP = '/dors/meilerlab/home/sliwosgr/udn_pipeline/uniprot_gene_id.tab'
        self._gene = gene
        self._mutations = mutations
        self._fasta = {}
        LOGGER.info("uniprot.__init__  completed after setting:\n%s"%pprint.pformat(self.__dict__))
    
    def check_aa(self,unp):
        warnings = []
        mutpos = [x[1] for x in self._mutations]
        seenmuts = set()
        if unp not in self._fasta:
            return [False,"Warning: failed to retrieve fasta for %s; %s" % (self._gene, unp)]
        for x in xrange(1,len(self._fasta[unp])+1):
            if x in mutpos:
                for y in self._mutations:
                    if y[1]==x:
                        seenmuts.add(y[1])
                        if y[0]!=self._fasta[unp][x-1]:
                            warnings.append(["Warning: AA at pos %d (%s) in %s fasta does not match WT AA in mutation: %s" % (x,self._fasta[unp][x-1],unp,"".join(str(z) for z in y))])
        for x in mutpos:
            if x>(len(self._fasta[unp])):
                warnings.append(["Warning: mutation position %d falls outside the %s fasta length of %d" % (x,unp,len(self._fasta[unp]))])
        return warnings

    def run(self,ignorewhole = True):
        unp = []
        try:
            assert os.path.isfile(self._UNP_MAP)
        except AssertionError:
            return [False,"%s does not exist" % self._UNP_MAP]
        with open(self._UNP_MAP) as infile:
            for line in infile:
                cur = line.strip().split()
                if cur[0]==self._gene.upper():
                    unp.append(cur[1])
        try:
            assert len(unp) > 0
        except AssertionError:
            return [False,"Gene %s not found in uniprot map %s" % (self._gene,self._UNP_MAP)]
        final_set = [self.get_annotations(unp[x],ignorewhole)[1]+"\n------\n" if x<len(unp)-1 else self.get_annotations(unp[x],ignorewhole)[1] for x in xrange(len(unp))]
#        if not outcome:
#            return [False,final_set]
        try:
            assert len(final_set) > 0
        except AssertionError:
            return [False,"Failed to get anything at all"]
        return [True,final_set]
    
    def get_annotations(self,unp,ignorewhole):
        annotations = ["UniprotID\tAnnotation\tStart_Res\tEnd_Res\tDescription"]
        respos = [x[1] for x in self._mutations]
        final_set = []
        fastafile = "%s.fasta" % (self._gene)
        seen_annotations = set()
        urladd = "http://www.uniprot.org/uniprot/%s.txt" % unp
        url = urlopen(urladd)
        raw_data = [x.rstrip() for x in url.readlines() if len(x.strip())>0]
        url.close()
        try:
            assert len(raw_data)>0
        except AssertionError:
            return [False,"Error with %s: Failed to load anything from url" % unp]
        try:
            length = int(raw_data[0].split()[3])
        except:
            return [False,"Error with %s: expected protein length at ix 3, got %s instead" % (unp,raw_data[0].split()[3])]

        within = False
        runner = []
        infasta = False
        fasta = ''
        for line in raw_data[1:]:
            if line[:2] == 'SQ':
                infasta = True
                continue
            elif infasta:
                if line[0]==' ':
                    fasta += "".join(line.strip().split())
                    continue
                else:
                    infasta = False
                    continue
            current = line.split()
                   
            if len(line)>4 and line[5] == ' ':
                if within:
                    cont = " ".join(current[1:])
                    runner[-1]+=(" "+cont)
                continue           
               
            if len(current)>=4:
                if within and len(runner)>0:
                    final_set.append(runner)
                    runner=[]
                    within = False
                category = current[0]
                subcat = current[1]
                if category not in self._ANNOTATIONS.keys() or subcat not in self._ANNOTATIONS[category]:
                    continue
                try:               
                    if current[2][0] == '<' or current[2][0]=='?':
                        start = int(current[2][1:])
                        posdesc = "uncertain start: %s; " % current[2]
                    else:
                        start = int(current[2])
                        posdesc = ""
                    if current[3][0] == '>' or current[3][0]=='?':
                        end = int(current[3][1:])
                        posdesc += "uncertain end: %s; " % current[2]
                    else:
                        end = int(current[3])
                except:
                    print "Warning: expected residue numbers at pos 2 and 3, got %s and %s instead, skipping this entry" % (current[2],current[3])
                    within = False
                    continue
                if start==1 and end==length and subcat not in self._PAIRS and ignorewhole:
                    within = False
                    continue
                if subcat not in self._PAIRS:
                    keep = False
                    for item in respos:
                        if start<=item and end>=item:
                            keep = True
                    if keep:
                        try:
                            desc = posdesc+" ".join(current[4:])
                        except:
                            desc = posdesc
                        runner = [subcat,start,end,desc]
                        seen_annotations.add(subcat)
                        within = True
                else:
                    keep = False
                    for item in respos:
                        if start==item or end==item:
                            keep = True
                    if keep:
                        try:
                            desc = posdesc+" ".join(current[4:])
                        except:
                            desc = posdesc
                        runner = [subcat,start,end,desc]
                        seen_annotations.add(subcat)
                        within = True                   
        self._fasta[unp] = fasta
        warnings = self.check_aa(unp)
        if len(final_set)>0:
            annotations += ["\t".join(str(y) for y in [unp]+x) for x in final_set]
        else:
            annotations = ["No annotations for residue(s): %s" % ",".join(str(x) for x in respos)]
        #    sys.exit("No annotations for residue %s in gene %s; uniprot %s" % (pos,gene,unp))
        
        final = ["Gene: %s\tUniprot: %s\tLength %d AA\n" % (self._gene,unp,length)]+annotations
        if len(warnings)>0:
            final += ["*********Warnings*********"]
            for item in warnings:
                final+=item
        if len(annotations)>1:
            final += ["\nDescriptions:"]
            for item in sorted(seen_annotations):
                final += ["%s: %s" % (item,self._DESCRIPTIONS[item])]
                
        if len(fasta)>0:
            with open(fastafile,'w') as fout:
                fout.write(">%s; %s\n" % (self._gene,unp))
                fout.write(fasta)
                fout.write("\n")
        else:
            print "Warning: failed to retrieve fasta for gene %s" % self._gene
        return [True,"\n".join(final)]

