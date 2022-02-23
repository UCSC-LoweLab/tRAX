#!/usr/bin/env python

import pysam
import sys
import argparse
import string
import itertools
from collections import defaultdict
import os.path
from trnasequtils import *
import time
import numpy as np
from scipy import stats
import random
from scipy.spatial.distance import euclidean



'''
~/pythonsource/TRAX/modelmismatch.py --covfile=mousetgirtplusalkbmism/mousetgirtplusalkbmism-coverage.txt --trnafile=~/pythonsource/trnatest/trnadbs/mm10_new/mm10-trnatable.txt

~/pythonsource/TRAX/modelmismatch.py --covfile=../mousetgirtminus/mousetgirtminus-coverage.txt --trnafile=~/pythonsource/trnatest/trnadbs/mm10_new/mm10-trnatable.txt  --samplefile=../mousetgirtminus.txt --bytrna  >trnacompare.txt

~/pythonsource/TRAX/modelmismatch.py --covfile=thioseq/thioseq-coverage.txt --trnafile=~/pythonsource/trnatest/trnadbs/hg19/hg19-trnatable.txt  --samplefile=thiosamples.txt --bytrna  >trnacompare.txt


 ~/pythonsource/TRAX/modelmismatch.py --covfile=thioseq/thioseq-coverage.txt --trnafile=~/pythonsource/trnatest/trnadbs/hg19/hg19-trnatable.txt  --samplefile=thiosamples.txt  --runname=thioseq >trnacompare.txt
 
~/pythonsource/TRAX/modelmismatch.py --covfile=angko/angko-coverage.txt --trnafile=~/pythonsource/trnatest/trnadbs/hg19/hg19-trnatable.txt  --samplefile=angkosamples2.txt  --runname=angko2 

~/pythonsource/TRAX/modelmismatch.py --covfile=sharma2018/sharma2018-coverage.txt --trnafile=~/pythonsource/trnatest/trnadbs/mm10/mm10-trnatable.txt  --samplefile=sharma2018samples.txt  --runname=sharma2018mism 

~/pythonsource/TRAX/modelmismatch.py --covfile=/projects/lowelab/users/abagi/data/trax_data/1_22_2020_c305_dmseq/trax_output_hg38_hek293ctrl/trax_output_hg38_hek293ctrl-coverage.txt --trnafile=~/pythonsource/trnatest/trnadbs/hg38/hg38-trnatable.txt  --samplefile=/projects/lowelab/users/abagi/data/trax_data/1_22_2020_c305_dmseq/samples.txt  --runname=brainletsamples

~/pythonsource/TRAX/modelmismatch.py --covfile=sharma2017samples/sharma2017samples-coverage.txt --trnafile=~/pythonsource/trnatest/trnadbs/mm10/mm10-trnatable.txt  --samplefile=sharma2017_samples.txt  --runname=sharma2017mism 

~/pythonsource/TRAX/modelmismatch.py --covfile=cruzsamples/cruzsamples-coverage.txt --trnafile=~/pythonsource/trnatest/trnadbs/mm10/mm10-trnatable.txt  --samplefile=cruzsamples.txt  --runname=cruzsamplesmism 


~/pythonsource/TRAX/modelmismatch.py --covfile=ottrfullmouse/ottrfullmouse-coverage.txt --trnafile=~/pythonsource/trnatest/trnadbs/mm10/mm10-trnatable.txt  --samplefile=mousevirusfull.txt  --runname=ottrfullmousec


~/pythonsource/TRAX/modelmismatch.py --covfile=humanexosome/humanexosome-coverage.txt --trnafile=~/pythonsource/trnatest/trnadbs/hg19/hg19-trnatable.txt  --samplefilehumanexosomesamples.txt  --runname=humanexosome

~/pythonsource/TRAX/modelmismatch.py --covfile=humanexosome/humanexosome-coverage.txt --trnafile=~/pythonsource/trnatest/trnadbs/hg19/hg19-trnatable.txt  --samplefilehumanexosomesamples.txt  --runname=humanexosome

~/pythonsource/TRAX/modelmismatch.py --covfile=epidtypesamples/epidtypesamples-coverage.txt --trnafile=~/pythonsource/trnatest/trnadbs/mm10/mm10-trnatable.txt  --samplefile=epidtypesamples.txt  --runname=epidtype

~/pythonsource/TRAX/modelmismatch.py --covfile=spermatazoatype/spermatazoatype-coverage.txt --trnafile=~/pythonsource/trnatest/trnadbs/mm10/mm10-trnatable.txt  --samplefile=spermatazoasamples.txt  --runname=spermatazoatype



~/pythonsource/TRAX/modelmismatch.py --covfile=miseqcond/miseqcond-coverage.txt --trnafile=~/pythonsource/trnatest/trnadbs/mm10/mm10-trnatable.txt  --samplefile=miseqsamples.txt  --runname=miseqsamples

'''
gapchars = set("-._~")
positions = list(str(curr) for curr in list([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,'17a',18,19,20,'20a','20b',21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14','e15','e16','e17','e18','e19',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76]))



def readclusterfile(filename):
    trnaclusts = defaultdict(set)
    for currline in open(filename):
        fields = currline.split()
        trnaclusts[int(fields[1])].add(fields[0].strip('"')) 
    return trnaclusts
def eucliddistance(firfreqs, secfreqs):
    
    return np.linalg.norm(np.array(firfreqs)-np.array(secfreqs))

def freqtolist(freqdict):
    for currbase in ["A","T","C","G", "-"]:
        yield freqdict.freqlists[currbase]
    

def getrandfreqs(freqsum, numfreqs):
    randlist = list(random.randrange(freqsum) for i in range(0,numfreqs - 1))  
    randlist.extend([0,freqsum])
    randlist.sort()
    results = list()
    for i in range(1, len(randlist)):
        results.append(randlist[i] - randlist[i - 1])
    return results
        
def drawcounts(fircounts, length, pcounts = .5):
    probtable = list((curr+pcounts)/(1.*(sum(fircounts)+ len(fircounts)*pcounts)) for curr in fircounts)
    #print >>sys.stderr, probtable
    #print >>sys.stderr, range(0, 3)
    #countlist = list(np.random.choice(4, length, p=probtable))
    #return countlist.count(0),countlist.count(1),countlist.count(2),countlist.count(3)
    
    results = np.random.multinomial(length, probtable)
    #print >>sys.stderr,fircounts
    #print >>sys.stderr,length
    #print >>sys.stderr,list(results)
    #print >>sys.stderr,sum(list(results))
    #print >>sys.stderr,"**"
    return list(results)
    
                 
        
class covpos:
    def __init__(self, Sample, Feature, position):
        self.Sample = Sample
        self.Feature = Feature
        self.position = position
    def __eq__(self, other):
       return self.Sample == other.Sample and self.Feature == other.Feature and self.position == other.position
    def __hash__(self):
       return hash(self.Sample) + hash(self.Feature)+hash(self.position)
           
class covline:
    def __init__(self, acounts, ccounts, gcounts,tcounts, deletions, actualbase, percentunique = None):
        self.acounts = int(acounts )
        self.ccounts = int(ccounts )
        self.gcounts = int(gcounts )
        self.tcounts = int(tcounts )
        self.deletions = int(deletions)
        self.actualbase = actualbase
        self.percentunique = percentunique
        if self.actualbase == "A":
            self.acounts
        elif self.actualbase == "T":
            self.tcounts    
        elif self.actualbase == "C":
            self.ccounts                      
        elif self.actualbase == "G":
            self.gcounts
        else:
            print >>sys.stderr, "No base "+self.actualbase 
            sys.exit(1)
        
        self.total = self.acounts +self.ccounts +self.gcounts +self.tcounts +self.deletions 
    def percentdict(self, pseudocounts = .1):
        return {"apercent":self.acounts/(self.total +pseudocounts)   ,"cpercent":self.ccounts/(self.total +pseudocounts)   ,"gpercent":self.gcounts/(self.total +pseudocounts)   ,"tpercent":self.tcounts/(self.total +pseudocounts), "delpercent":self.delcounts/(self.total +pseudocounts)   }
    def basecounts(self):
        return [self.acounts,self.ccounts,self.gcounts,self.tcounts, self.deletions]
    def percentcounts(self, pseudocounts = .01):
        return [self.acounts/(self.total +pseudocounts),self.ccounts/(self.total +pseudocounts),self.gcounts/(self.total +pseudocounts),self.tcounts/(self.total +pseudocounts)]
    def totalcounts(self):
        return sum([self.acounts,self.ccounts,self.gcounts,self.tcounts])
    def mismatchpercent(self, pseudocounts = .01):
        if self.actualbase == "A":
            return (self.total + pseudocounts - self.acounts)/(self.total +pseudocounts)     
        if self.actualbase == "T":
            return (self.total + pseudocounts -  self.tcounts)/(self.total+pseudocounts)     
        if self.actualbase == "C":
            return (self.total + pseudocounts - self.ccounts)/(self.total+pseudocounts )                        
        if self.actualbase == "G":
            return (self.total + pseudocounts - self.gcounts)/(self.total +pseudocounts)                
    def matchpercent(self, pseudocounts = .01):
        return 1 - self.mismatchpercent(pseudocounts = pseudocounts)
    def shufflecounts(self):
        shuffledcounts = getrandfreqs(self.totalcounts(), 4)
        #print >>sys.stderr, self.basecounts()
        #print >>sys.stderr, sum(self.basecounts())
        #print >>sys.stderr, shuffledcounts
        #print >>sys.stderr, sum(shuffledcounts)
        #print >>sys.stderr, "***"
        
        return covline(shuffledcounts[0],shuffledcounts[1],shuffledcounts[2],shuffledcounts[3],self.deletions, self.actualbase,percentunique =  self.percentunique)
    def drawcomparison(self, firdata):
        #print >>sys.stderr, self.basecounts()
        #print >>sys.stderr, sum(self.basecounts())
        #print >>sys.stderr, secdata.basecounts()
        #print >>sys.stderr, sum(secdata.basecounts())

        newcounts = drawcounts(firdata.basecounts(),self.totalcounts())
        #print >>sys.stderr, sum(self.basecounts())
        #print >>sys.stderr, sum(newcounts)
        #print >>sys.stderr, sum(firdata.basecounts())
        #print >>sys.stderr, "***"
        return covline(newcounts[0],newcounts[1],newcounts[2],newcounts[3],self.deletions, self.actualbase,percentunique =  self.percentunique)
                   
class positiondist:
    def __init__(self,nuclist):
        self.freqlists = nuclist
        #print >>sys.stderr, nuclist.keys()
        self.length = max(len(nuclist[curr]) for curr in nuclist.iterkeys())
        
    def getmeans(self):
        for currbase in freqlists:
            amean = sum(self.freqlists["A"]) / self.length
            cmean = sum(self.freqlists["C"]) / self.length
            tmean = sum(self.freqlists["T"]) / self.length
            gmean = sum(self.freqlists["G"]) / self.length
            delmean = sum(self.freqlists["0"]) / self.length
        return {"A":amean,"T":tmean,"C":cmean,"G":gmean,"-":delmean}
    def printtable(self, output = sys.stdout):
        print "\t".join(["A","T","C","G"])
        for i in range(length):
            print "\t".join([self.freqlists["A"][i],self.freqlists["T"][i],self.freqlists["C"][i],self.freqlists["G"][i]])


def scalecounts(secset, firset):
    scalingfactor = sum(firset)/(1.*sum(secset))
    newsecset = list(curr*scalingfactor for curr in secset)
    return newsecset

def chitest(firset, secset, pcounts = .01, scaled = False):
    fircounts = list(curr+pcounts for curr in firset)
    if scaled:
        seccounts = list(curr+pcounts for curr in secset)
        seccounts = scalecounts(seccounts,firset)
    else:
        seccounts = list(curr+pcounts for curr in secset)
    return stats.chisquare(fircounts, seccounts)



def countentropy(firset, secset, pcounts = .01):
    firset = list(curr+pcounts for curr in firset)
    secset = list(curr+pcounts for curr in secset)
    firprobs = list((curr)/(1.*sum(firset)) for curr in firset) #+ pcounts*len(firset)
    secprobs = list((curr)/(1.*sum(secset)) for curr in secset)
    return stats.entropy(firprobs, secprobs)
    
def log2(num):
    return math.log(num, 2)
    
def bhattacharyyadistance(firset, secset):    
    """Calculates Bhattacharyya distance (https://en.wikipedia.org/wiki/Bhattacharyya_distance)."""
    sim = - np.log(np.sum([np.sqrt(p*q) for (p, q) in zip(firset, secset)]))
    assert not np.isnan(sim), 'Error: Similarity is nan.'
    if np.isinf(sim):
        # the similarity is -inf if no term in the review is in the vocabulary
        return 0
    return sim 
def bhatcounts(firset, secset, pcounts = .01):
    firset = list(curr+pcounts for curr in firset)
    secset = list(curr+pcounts for curr in secset)
    firprobs = list((curr)/(1.*sum(firset)) for curr in firset) #+ pcounts*len(firset)
    secprobs = list((curr)/(1.*sum(secset)) for curr in secset)
    return bhattacharyyadistance(firprobs, secprobs)
    
sqrt2 = np.sqrt(2)
def hellingerdistance(firset, secset):    

    return euclidean(np.sqrt(firset), np.sqrt(secset))/sqrt2 
def hellingercounts(firset, secset, pcounts = .01):
    firset = list(curr+pcounts for curr in firset)
    secset = list(curr+pcounts for curr in secset)
    firprobs = list((curr)/(1.*sum(firset)) for curr in firset) #+ pcounts*len(firset)
    secprobs = list((curr)/(1.*sum(secset)) for curr in secset)
    return hellingerdistance(firprobs, secprobs)

class covcomparison:
    def __init__(self, firpos, secpos, firdata, secdata):
        self.firpos = firpos
        self.secpos = secpos
        self.firdata = firdata
        self.secdata = secdata
        
    def chisquare(self, pcounts = .01, mincount = 50):
        #scalingfactor = sum(self.firdata.basecounts())/(1.*sum(self.secdata.basecounts()))
        if self.hascounts(mincount = mincount):
            #print self.firdata.basecounts()
            #print sum(self.firdata.basecounts())
            #
            #print self.secdata.basecounts()
            #print sum(self.secdata.basecounts())
            return chitest(self.firdata.basecounts(),self.secdata.basecounts(), pcounts = pcounts)
        else:
            return "NA", 1
            
    def reversepair(self):
        return covcomparison(self.secpos, self.firpos, self.secdata, self.firdata)
    def eucliddist(self):
        return eucliddistance(self.firdata.basecounts(), self.secdata.basecounts())
    def hascounts(self, mincount = 50):
        return sum(self.firdata.basecounts()) > mincount and sum(self.secdata.basecounts()) > mincount
    def countentropy(self):
        return countentropy(self.firdata.basecounts(), self.secdata.basecounts())
    def bhatdistance(self):
        return bhatcounts(self.firdata.basecounts(), self.secdata.basecounts())        
    def hdistance(self):
        return hellingercounts(self.firdata.basecounts(), self.secdata.basecounts())        
        
    def sameorigbase(self):
        return self.firdata.actualbase == self.secdata.actualbase
    def bothhavebase(self, base):
        return self.firdata.actualbase == base and self.secdata.actualbase == base
    def bothunique(self, threshold = .9):
        if self.firdata.percentunique is None or self.secdata.percentunique is None:
            return False
        return self.firdata.percentunique > threshold and  self.secdata.percentunique > threshold
        
    def containsmismatches(self, threshold = .01):
        return self.firdata.mismatchpercent()  > threshold or self.secdata.mismatchpercent() > threshold
    def containsbothmismatches(self, threshold = .01):
        return self.firdata.mismatchpercent()  > threshold and self.secdata.mismatchpercent() > threshold
    def compareprint(self):
        
        pval, chiscore = self.chisquare()
        eucdist = self.eucliddist()
        entropy = self.countentropy()
        print self.firdata.basecounts()
        print self.secdata.basecounts()  
        #print list(scalecounts(self.secdata.basecounts(), self.firdata.basecounts()))
        print "\t".join([self.firpos.Sample, self.firpos.Feature,self.firpos.position,self.secpos.Sample, self.secpos.Feature,self.secpos.position,str(eucdist),str(entropy),str(pval),str(chiscore)])
    def shufflecomparison(self):
        newfirst = self.firdata.shufflecounts()
        newsec = self.secdata.shufflecounts()
        return covcomparison(self.firpos, self.secpos, newfirst, newsec)
        
    def redrawcomparison(self):
        newfirst = self.firdata
        newsec = self.secdata.drawcomparison(self.firdata)
        return covcomparison(self.firpos, self.secpos, newfirst, newsec)
class covdata:
    def __init__(self):
        self.covdict = dict()
        self.nucpos = defaultdict(dict)
    def getpos(self,currpos):
        return self.covdict[currpos]
    def addline(self, Sample, Feature, position,apercent, cpercent, gpercent,tpercent, deletions, actualbase, percentunique = None):
        #print "||"+ Sample+ " "+Feature +" "+ position
        self.covdict[covpos( Sample, Feature, position)] = covline(apercent, cpercent, gpercent,tpercent, deletions, actualbase, percentunique =  percentunique)
        self.nucpos[Feature][position] = actualbase
        
    def getmismatchpos(self, featlist, samplelist,poslist, minmismatch = .2):
        posset = set()
        for currpos in poslist:
            for currsample in samplelist:
                for currfeat in featlist:
                    currcovpos = covpos(currsample, currfeat, currpos)
                    if currcovpos in self.covdict and self.covdict[currcovpos].mismatchpercent() > minmismatch:
                        posset.add(currpos)
                        break
        return posset
    def comparetrnas(self, featlist, Sample, position):
        featpercents = dict()
        for currfeat in featlist:
            currpos = covpos( Sample, currfeat, position)
            if currpos in self.covdict:
                featpercents[currfeat] = self.covdict[currpos].percentdict()
            
        return featpercents
        

    def getposset(self, samplelist, featurelist, positionlist, currbase = None):
        covposlist = list()

        for currpos in positionlist:
            for currfeat in featurelist:
                 for currsample in samplelist:
                     #print currsample+" "+currfeat+" "+currpos
                     currcovpos = covpos( currsample, currfeat, currpos)
                     #print >>sys.stderr, currpos
                     if currpos == 56:
                         #print >>sys.stderr, "**||"
                         pass
                     if currcovpos in self.covdict and (currbase is None or self.covdict[currcovpos].actualbase == currbase):
                         
                         covposlist.append(currcovpos)
        return covposlist
        
    def compareposset(self, covposlist):
                
        for firpos, secpos in itertools.combinations(covposlist, 2):
            yield covcomparison(firpos, secpos, self.covdict[firpos], self.covdict[secpos] ) 

def readcovfile(covfile):
    headers = None
    headerdict = None
    totalcount = 0
    badcount = 0
    covcounts = covdata()
    for linenum, currline in enumerate(open(covfile)):
        fields = currline.split()
        if linenum == 0:
            headers = fields
            headerdict = {headers:i for i, headers in enumerate(headers)}
            #print headerdict["Sample"]
        elif len(fields) < 2:
            continue
        else:
            Sample = fields[headerdict["Sample"]]
            Feature = fields[headerdict["Feature"]]
            position = fields[headerdict["position"]]
            #print fields[headerdict["Sample"]]
            gpercent = fields[headerdict["guanines"]]
            cpercent = fields[headerdict["cytosines"]]			
            tpercent = fields[headerdict["thymines"]]
            apercent = fields[headerdict["adenines"]]
            
            actualbase= fields[headerdict["actualbase"]].upper()
            mismatchedbases= fields[headerdict["mismatchedbases"]]
            deletions = fields[headerdict["deletions"]]
            coverage= float(fields[headerdict["coverage"]])
            uniquereads = float(fields[headerdict["uniquecoverage"]])
            totalcount += 1
            if actualbase in gapchars:
                #print >>sys.stderr, currline
                continue
            if actualbase == "U":
                actualbase = "T"
            covcounts.addline(Sample, Feature, position,apercent, cpercent, gpercent,tpercent, deletions, actualbase, percentunique = (uniquereads)/(1.*coverage + .1))
            #print  str(float(coverage))+":"+ str(float(gpercent) + float(cpercent) + float(tpercent) + float(apercent)+ float(deletions))
            #if float(coverage) !=  :
                #print >>sys.stderr, currline
                             
                             
    return covcounts                           


nucbases = ["A","T","C","G"]

'''
data = read.table("comparetrnas.txt",header = TRUE,row.names = NULL, stringsAsFactors=FALSE)
sortdata = data[order(-data$hdist),]
head(sortdata[sortdata$firname == "M_dm_Heart_M6_minusAlkB" & !(sortdata$firfeat %in% excludelist) & !(sortdata$secfeat %in% excludelist),])


head(sortdata[sortdata$groupname == "M_dm_Heart_M6_minusAlkB_34pos_C_clust" & sortdata$firfeat %in% c("tRNA-Leu-CAG-2","tRNA-Leu-CAG-1") &  sortdata$secfeat %in% c("tRNA-Leu-CAG-2","tRNA-Leu-CAG-1"),])

unique(sortdata[sortdata$firpercent < .9 & sortdata$secpercent < .9 & sortdata$hdist > .5,"groupname"])  
unique(sortdata[sortdata$hdist > .5,"groupname"])  


data = read.table("comparesamples.txt",header = TRUE,row.names = NULL, stringsAsFactors=FALSE)
sortdata = data[order(-data$hdist),]
unique(sortdata[sortdata$hdist > .5,"groupname"])  

''' 

def getreplicates(positions, trnainfo,sampleinfo):
    for currpos in positions:
        for curramino in trnainfo.allaminos():
            for currtrna in trnainfo.getaminotranscripts(curramino):
                for currreplicate in sampleinfo.allreplicates():
                    clustname = currreplicate + "_"+currpos+"pos"
                    yield clustname, sampleinfo.getrepsamples(currreplicate),[currtrna], [currpos] 
                    
def gettrnasamples(positions, trnainfo,sampleinfo):
    for currpos in positions:
        #print >>sys.stderr, currpos
        for currsample in sampleinfo.getsamples():
        #for currsample in ["M_dm_Heart_M5_minusAlkB"]:
        #for currsample in ["M_dm_Liver_M5_minusAlkB"]:
            for curramino in trnainfo.allaminos():
                clustname = currsample +"_"+curramino+ "_"+currpos+"pos"
                yield clustname, [currsample],trnainfo.getaminotranscripts(curramino), [currpos]  
                
def getsamples(positions, trnainfo,sampleinfo):
    for currpos in positions:
        for curramino in trnainfo.allaminos():
            for currtrna in trnainfo.getaminotranscripts(curramino):
                clustname = currtrna + "_"+currpos+"pos"
                yield clustname, sampleinfo.getsamples(),[currtrna], [currpos] 
                
                
def twopos(number):
    return "{:.2f}".format(number)
    
    
def gettrnainfo(outfile, covcounts, positions, trnainfo,sampleinfo,mismatchthreshold = .1, mincount = 50):
    
    totalcounts = defaultdict(lambda: defaultdict(int))

    mismatchcounts = defaultdict(lambda: defaultdict(int))
    mismatchpercents = defaultdict(lambda: defaultdict(list))
    mismatchlists = defaultdict(lambda: defaultdict(list))
    for groupname, samplelist, trnalist, poslist in getsamples(positions, trnainfo,sampleinfo):
        
        for currbase in nucbases:
            #print >>sys.stderr, poslist
            array = 1
            posset = covcounts.getposset(samplelist,trnalist, poslist, currbase )
            
            for currpos in posset:
                if covcounts.getpos(currpos).totalcounts() < mincount:
                    continue
                totalcounts[currpos.Feature][currpos.position] += 1
                mismatchpercents[currpos.Feature][currpos.position].append(covcounts.getpos(currpos).mismatchpercent())
                mismatchlists[currpos.Feature][currpos.position].append(covcounts.getpos(currpos).totalcounts())
                
                if covcounts.getpos(currpos).mismatchpercent()  > mismatchthreshold:                    
                    #print >>sys.stderr, covcounts.getpos(currpos).mismatchpercent()
                    mismatchcounts[currpos.Feature][currpos.position] += 1
                    
                if currpos.Feature == "tRNA-Ile-TAT-1" and currpos.position == "45":
                    #print >>sys.stderr, covcounts.getpos(currpos).mismatchpercent()
                    #print >>sys.stderr, mismatchpercents[currpos.Feature][currpos.position]
                    pass
    print >>outfile, "\t".join(currposition for currposition in positions)
    positionmismatches = defaultdict(set)
    for currposition in positions:
        for currtrna in trnainfo.gettranscripts():
        
            if mismatchcounts[currtrna][currposition] > 1:
                positionmismatches[currposition].add(currtrna)
            if mismatchcounts[currtrna][currposition] != 0 and mismatchcounts[currtrna][currposition] != totalcounts[currtrna][currposition]:
                #print >>sys.stderr, currtrna
                #print >>sys.stderr, currposition
                #print >>sys.stderr, str(mismatchcounts[currtrna][currposition])+"/"+str(totalcounts[currtrna][currposition])
                #print >>sys.stderr, mismatchpercents[currtrna][currposition]
                #print >>sys.stderr, mismatchlists[currtrna][currposition]
                
                pass
        print >>outfile, currtrna +"\t"+"\t".join(str(mismatchcounts[currtrna][currposition])+"/"+str(totalcounts[currtrna][currposition]) for currposition in positions)
   
    for currposition in positions:
        if len(positionmismatches[currposition]) > 1:
            print >>sys.stderr, currposition+":"+str(positionmismatches[currposition])
        
        
def createtable(outfile, covcounts, pairgroup, minreads = 50, skipmatches = True, shufflemode = False, drawmode = False):
    allclusters = set()
    totalcounted = 0
    skipunique = 0
    entropies = list()


    #print mismatchlocs
    postotal = defaultdict(int)
    posmismatch = defaultdict(int)
    
    
    print >>outfile, "\t".join(["groupname","firname", "firfeat","firpos","firrefbase","firtotal","firpercent","fircounts","secname", "secfeat","secpos","secrefbase","sectotal","secpercent","seccounts","entropy","bdist","hdist","pval","chiscore"])
    
    for groupname, samplelist, trnalist, poslist in pairgroup:
        
        for currbase in nucbases:
            #print >>sys.stderr, poslist
            array = 1
            posset = covcounts.getposset(samplelist,trnalist, poslist, currbase )
            
            for currpair in covcounts.compareposset(posset):
               
               
               currgroupname = groupname + "_"+currbase+"_clust"
               testgroup = "tRNA-Ala-TGC-5_58pos_A_clust"  #
               allclusters.add(currgroupname)
               #print >>sys.stderr, ",".join(curr.Feature for curr in poslist)
               #print >>sys.stderr, currgroupname
               
               if currgroupname != testgroup:
                   pass
                   #continue
               
               if not currpair.bothhavebase(currbase):
                   #print >>sys.stderr, "bothhavebase"
                   continue
               
               '''
               if currpair.firpos.Feature == "tRNA-Met-CAT-6" and currpair.secpos.Feature != 'tRNA-Gly-GCC-3':
                   
                   print >>sys.stderr, currpair.secpos.Feature
                   print >>sys.stderr, currpair.firdata.actualbase
                   print >>sys.stderr, currpair.secdata.actualbase
                   #self.firdata.actualbase
               if currpair.secpos.Feature == "tRNA-Met-CAT-6" and currpair.firpos.Feature != 'tRNA-Gly-GCC-3':
                   print >>sys.stderr, currpair.firpos.Feature 
                   print >>sys.stderr, currpair.firdata.actualbase
                   print >>sys.stderr, currpair.secdata.actualbase
                   
               '''

               if not currpair.hascounts(mincount = minreads):
                   #print >>sys.stderr, "hascounts"
                   #print >>sys.stderr, minreads
                   continue
               

               if skipmatches and not currpair.containsbothmismatches():
                   continue
               
               if shufflemode:
                   currpair = currpair.shufflecomparison()
               elif drawmode:
                   currpair = currpair.redrawcomparison()
               
               totalcounted += 1    
               if not currpair.bothunique():
                   skipunique += 1
                   pass
               #currpair.compareprint() 
               #print >>sys.stderr, "["+",".join(str(curr) for curr in currpair.firdata.basecounts())+"]" "["+",".join(str(curr) for curr in currpair.secdata.basecounts())+"]" +str(currpair.firdata.mismatchpercent())+":"+str(str(currpair.secdata.mismatchpercent()))
               
               
               currentropy = currpair.countentropy()
               #reventropy = currpair.reversepair().countentropy()
               
               #reventropy = currpair.reversepair().countentropy()
               #print >>sys.stderr, str(currentropy)+":"+str(reventropy)
               entropies.append(currentropy)
               chiscore, pval = currpair.chisquare()
               revchiscore, revpval = currpair.reversepair().chisquare()
               #print >>sys.stderr, str(chiscore)+":"+str(revchiscore)
               bhatd = currpair.bhatdistance()
               #revbhatd = currpair.reversepair().bhatdistance()            
               #print >>sys.stderr, str(bhatd)+":"+str(revbhatd)
               hdistance = currpair.hdistance()
               revhdist = currpair.reversepair().hdistance()            
               #print >>sys.stderr, str(hdistance)+":"+str(revhdist)
               postotal[currpair.firpos.position] += 1
               '''
               if currpair.firpos.Feature == "tRNA-Met-CAT-6" or currpair.secpos.Feature == "tRNA-Met-CAT-6":
                   print >>sys.stderr, poslist
                   print >>sys.stderr, "\t".join([currgroupname,currpair.firpos.Sample, currpair.firpos.Feature,currpair.firpos.position,str(currpair.firdata.totalcounts()),str(currpair.firdata.matchpercent()),"["+",".join(str(curr) for curr in currpair.firdata.basecounts())+"]",currpair.secpos.Sample, currpair.secpos.Feature,currpair.secpos.position,str(currpair.secdata.totalcounts()),str(currpair.secdata.matchpercent()),"["+",".join(str(curr) for curr in currpair.secdata.basecounts())+"]",str(currentropy),str(bhatd),str(hdistance),str(pval),str(chiscore)])

               if currpair.firpos.Feature == "tRNA-Gly-GCC-3" or currpair.secpos.Feature == "tRNA-Gly-GCC-3":
                   print >>sys.stderr, "\t".join([currgroupname,currpair.firpos.Sample, currpair.firpos.Feature,currpair.firpos.position,str(currpair.firdata.totalcounts()),str(currpair.firdata.matchpercent()),"["+",".join(str(curr) for curr in currpair.firdata.basecounts())+"]",currpair.secpos.Sample, currpair.secpos.Feature,currpair.secpos.position,str(currpair.secdata.totalcounts()),str(currpair.secdata.matchpercent()),"["+",".join(str(curr) for curr in currpair.secdata.basecounts())+"]",str(currentropy),str(bhatd),str(hdistance),str(pval),str(chiscore)])

               '''
               print >>outfile, "\t".join([currgroupname,currpair.firpos.Sample, currpair.firpos.Feature,currpair.firpos.position,currpair.firdata.actualbase,str(currpair.firdata.totalcounts()),str(currpair.firdata.matchpercent()),""+",".join(str(curr) for curr in currpair.firdata.basecounts())+"",currpair.secpos.Sample, currpair.secpos.Feature,currpair.secpos.position,currpair.secdata.actualbase,str(currpair.secdata.totalcounts()),str(currpair.secdata.matchpercent()),""+",".join(str(curr) for curr in currpair.secdata.basecounts())+"",str(currentropy),str(bhatd),str(hdistance),str(pval),str(chiscore)])
               #currpair.compareprint() 
               if currentropy > 1:
                   posmismatch[currpair.firpos.position] += 1
                   pass
               if pval < .05:
                   #currpair.compareprint() 
                   pass
               if chiscore > 50000: 
                   #
                   pass
           #print freqcounts.freqlists
           #repfreqs[currreplicate] = covcounts
           #print freqcounts.freqlists

            
def createcombinedtable(outfile, covcounts, pairgroup, minreads = 50, skipmatches = True, shufflemode = False, drawmode = False):
    allclusters = set()
    totalcounted = 0
    skipunique = 0
    entropies = list()


    #print mismatchlocs
    postotal = defaultdict(int)
    posmismatch = defaultdict(int)
    
    
    print >>outfile, "\t".join(["groupname","firname", "firfeat","firpos","firrefbase","firtotal","firpercent","fircounts","secname", "secfeat","secpos","secrefbase","sectotal","secpercent","seccounts","entropy","bdist","hdist","pval","chiscore"])
    
    for groupname, samplelist, trnalist, poslist in pairgroup:
        
        for currbase in nucbases:
            #print >>sys.stderr, poslist
            array = 1
            posset = covcounts.getposset(samplelist,trnalist, poslist, currbase )
            
            for currpair in covcounts.compareposset(posset):
               
               
               currgroupname = groupname + "_"+currbase+"_clust"
               allclusters.add(currgroupname)
               #print >>sys.stderr, ",".join(curr.Feature for curr in poslist)
               #print >>sys.stderr, currgroupname
               
               if currgroupname != testgroup:
                   pass
                   #continue
               
               if not currpair.bothhavebase(currbase):
                   #print >>sys.stderr, "bothhavebase"
                   continue
               
               '''
               if currpair.firpos.Feature == "tRNA-Met-CAT-6" and currpair.secpos.Feature != 'tRNA-Gly-GCC-3':
                   
                   print >>sys.stderr, currpair.secpos.Feature
                   print >>sys.stderr, currpair.firdata.actualbase
                   print >>sys.stderr, currpair.secdata.actualbase
                   #self.firdata.actualbase
               if currpair.secpos.Feature == "tRNA-Met-CAT-6" and currpair.firpos.Feature != 'tRNA-Gly-GCC-3':
                   print >>sys.stderr, currpair.firpos.Feature 
                   print >>sys.stderr, currpair.firdata.actualbase
                   print >>sys.stderr, currpair.secdata.actualbase
                   
               '''

               if not currpair.hascounts(mincount = minreads):
                   #print >>sys.stderr, "hascounts"
                   #print >>sys.stderr, minreads
                   continue
               

               if skipmatches and not currpair.containsbothmismatches():
                   continue
               
               if shufflemode:
                   currpair = currpair.shufflecomparison()
               elif drawmode:
                   currpair = currpair.redrawcomparison()
               
               totalcounted += 1    
               if not currpair.bothunique():
                   skipunique += 1
                   pass
               #currpair.compareprint() 
               #print >>sys.stderr, "["+",".join(str(curr) for curr in currpair.firdata.basecounts())+"]" "["+",".join(str(curr) for curr in currpair.secdata.basecounts())+"]" +str(currpair.firdata.mismatchpercent())+":"+str(str(currpair.secdata.mismatchpercent()))
               
               
               currentropy = currpair.countentropy()
               #reventropy = currpair.reversepair().countentropy()
               
               #reventropy = currpair.reversepair().countentropy()
               #print >>sys.stderr, str(currentropy)+":"+str(reventropy)
               entropies.append(currentropy)
               chiscore, pval = currpair.chisquare()
               revchiscore, revpval = currpair.reversepair().chisquare()
               #print >>sys.stderr, str(chiscore)+":"+str(revchiscore)
               bhatd = currpair.bhatdistance()
               #revbhatd = currpair.reversepair().bhatdistance()            
               #print >>sys.stderr, str(bhatd)+":"+str(revbhatd)
               hdistance = currpair.hdistance()
               revhdist = currpair.reversepair().hdistance()            
               #print >>sys.stderr, str(hdistance)+":"+str(revhdist)
               postotal[currpair.firpos.position] += 1
               '''
               if currpair.firpos.Feature == "tRNA-Met-CAT-6" or currpair.secpos.Feature == "tRNA-Met-CAT-6":
                   print >>sys.stderr, poslist
                   print >>sys.stderr, "\t".join([currgroupname,currpair.firpos.Sample, currpair.firpos.Feature,currpair.firpos.position,str(currpair.firdata.totalcounts()),str(currpair.firdata.matchpercent()),"["+",".join(str(curr) for curr in currpair.firdata.basecounts())+"]",currpair.secpos.Sample, currpair.secpos.Feature,currpair.secpos.position,str(currpair.secdata.totalcounts()),str(currpair.secdata.matchpercent()),"["+",".join(str(curr) for curr in currpair.secdata.basecounts())+"]",str(currentropy),str(bhatd),str(hdistance),str(pval),str(chiscore)])

               if currpair.firpos.Feature == "tRNA-Gly-GCC-3" or currpair.secpos.Feature == "tRNA-Gly-GCC-3":
                   print >>sys.stderr, "\t".join([currgroupname,currpair.firpos.Sample, currpair.firpos.Feature,currpair.firpos.position,str(currpair.firdata.totalcounts()),str(currpair.firdata.matchpercent()),"["+",".join(str(curr) for curr in currpair.firdata.basecounts())+"]",currpair.secpos.Sample, currpair.secpos.Feature,currpair.secpos.position,str(currpair.secdata.totalcounts()),str(currpair.secdata.matchpercent()),"["+",".join(str(curr) for curr in currpair.secdata.basecounts())+"]",str(currentropy),str(bhatd),str(hdistance),str(pval),str(chiscore)])

               '''
               print >>outfile, "\t".join([currgroupname,currpair.firpos.Sample, currpair.firpos.Feature,currpair.firpos.position,currpair.firdata.actualbase,str(currpair.firdata.totalcounts()),str(currpair.firdata.matchpercent()),""+",".join(str(curr) for curr in currpair.firdata.basecounts())+"",currpair.secpos.Sample, currpair.secpos.Feature,currpair.secpos.position,currpair.secdata.actualbase,str(currpair.secdata.totalcounts()),str(currpair.secdata.matchpercent()),""+",".join(str(curr) for curr in currpair.secdata.basecounts())+"",str(currentropy),str(bhatd),str(hdistance),str(pval),str(chiscore)])
               #currpair.compareprint() 
               if currentropy > 1:
                   posmismatch[currpair.firpos.position] += 1
                   pass
               if pval < .05:
                   #currpair.compareprint() 
                   pass
               if chiscore > 50000: 
                   #
                   pass
           #print freqcounts.freqlists
           #repfreqs[currreplicate] = covcounts
           #print freqcounts.freqlists

            
    
def main(**argdict):
    covfile = argdict["covfile"]
    skipmatches = argdict["skipperfect"]
    runname = argdict["runname"]
    minreads = 50
    if argdict["minreads"] is not None:
        minreads = int(argdict["minreads"])
        
        
    trnainfo = transcriptfile(os.path.expanduser(argdict["trnafile"]))
    sampleinfo = samplefile(os.path.expanduser(argdict["samplefile"]))
    
    #sliver = .000001
    #print >>sys.stderr, sum([1 - 3*sliver,sliver,sliver,sliver])
    #print >>sys.stderr, bhattacharyyadistance([1 - 3*sliver,sliver,sliver,sliver],[sliver,sliver,sliver,1 - 3*sliver])
    #sys.exit()
    
    covcounts = readcovfile(covfile)
    
    #print >>sys.stderr, covcounts.covdict[covpos("M_dm_Heart_M6_minusAlkB","tRNA-Leu-CAG-1","34")].basecounts()
    #print >>sys.stderr, covcounts.covdict[covpos("M_dm_Heart_M6_minusAlkB","tRNA-Leu-CAG-1","34")].mismatchpercent()
    #print >>sys.stderr, covcounts.covdict[covpos("M_dm_Heart_M6_minusAlkB","tRNA-Leu-CAG-1","34")].actualbase
    #print >>sys.stderr, covcounts.covdict[covpos("M_dm_Heart_M6_minusAlkB","tRNA-Leu-CAG-1","34")].ccounts
    #
    #print >>sys.stderr, covcounts.covdict[covpos("Mouse_Brain_M5_minusAlkB","tRNA-Ala-TGC-6","58")].percentcounts()
    #print >>sys.stderr, covcounts.covdict[covpos("Mouse_Brain_M5_minusAlkB","tRNA-Ala-TGC-6","58")].basecounts()
    #print >>sys.stderr, covcounts.covdict[covpos("Mouse_Brain_M5_minusAlkB","tRNA-Ala-TGC-6","58")].mismatchpercent()
    #sys.exit()tRNA-Gly-GCC-3
    #print >>sys.stderr, "{:.2f},{:.2f},{:.2f},{:.2f}".format(*covcounts.covdict[covpos("M_dm_Heart_M5_minusAlkB","tRNA-Gly-GCC-3","58")].percentcounts())
    #print >>sys.stderr, "{:.2f},{:.2f},{:.2f},{:.2f}".format(*covcounts.covdict[covpos("M_dm_Heart_M5_minusAlkB","tRNA-Tyr-GTA-4","58")].percentcounts())
    #print >>sys.stderr, "{:.2f},{:.2f},{:.2f},{:.2f}".format(*covcounts.covdict[covpos("M_dm_Heart_M5_minusAlkB","tRNA-Arg-CCG-2","58")].percentcounts())
    #print >>sys.stderr, "{:.2f},{:.2f},{:.2f},{:.2f}".format(*covcounts.covdict[covpos("M_dm_Heart_M5_minusAlkB","tRNA-Val-CAC-3","58")].percentcounts())
    #print >>sys.stderr, "{:.2f},{:.2f},{:.2f},{:.2f}".format(*covcounts.covdict[covpos("M_dm_Heart_M5_minusAlkB","tRNA-Leu-AAG-3","58")].percentcounts())
    #
    #print >>sys.stderr,covcounts.covdict[covpos("M_dm_Heart_M5_minusAlkB","tRNA-Gly-GCC-3","57")].actualbase + ":"+covcounts.covdict[covpos("M_dm_Heart_M5_minusAlkB","tRNA-Gly-GCC-3","59")].actualbase
    #print >>sys.stderr,covcounts.covdict[covpos("M_dm_Heart_M5_minusAlkB","tRNA-Tyr-GTA-4","57")].actualbase + ":"+covcounts.covdict[covpos("M_dm_Heart_M5_minusAlkB","tRNA-Tyr-GTA-4","59")].actualbase
    #print >>sys.stderr,covcounts.covdict[covpos("M_dm_Heart_M5_minusAlkB","tRNA-Arg-CCG-2","57")].actualbase + ":"+covcounts.covdict[covpos("M_dm_Heart_M5_minusAlkB","tRNA-Arg-CCG-2","59")].actualbase
    #print >>sys.stderr,covcounts.covdict[covpos("M_dm_Heart_M5_minusAlkB","tRNA-Val-CAC-3","57")].actualbase + ":"+covcounts.covdict[covpos("M_dm_Heart_M5_minusAlkB","tRNA-Val-CAC-3","59")].actualbase
    #print >>sys.stderr,covcounts.covdict[covpos("M_dm_Heart_M5_minusAlkB","tRNA-Leu-AAG-3","57")].actualbase + ":"+covcounts.covdict[covpos("M_dm_Heart_M5_minusAlkB","tRNA-Leu-AAG-3","59")].actualbase
    
    #print scaledchitest([0, 4, 0, 232],[0, 0, 0, 105])
    
        
    clusttest = False
    if clusttest:
        #trnaclusts = readclusterfile("/projects/lowelab/users/holmes/pythonsource/trnatest/test/poscompare/M_dm_Heart_M5_minusAlkB_58pos_A_clust_groups.txt")
        trnaclusts = readclusterfile("/projects/lowelab/users/holmes/pythonsource/trnatest/test/poscompare/tRNA-Ala-TGC-6_48pos_A_clust_groups.txt")
        for currclust in trnaclusts.iterkeys():
            print  "cluster "+str(currclust)
            lista = list()
            listc = list()
            listg = list()
            listt = list()
            upbases = defaultdict(int)
            downbases = defaultdict(int)
            aminos = defaultdict(int)
            
            for currtrna in trnaclusts[currclust]:
                #print  currtrna
                #trnainfo = covcounts.covdict[covpos("M_dm_Heart_M5_minusAlkB",currtrna,"58")]
                trnainfo = covcounts.covdict[covpos(currtrna,"tRNA-Ala-TGC-6","58")]
                #print  "{:.2f},{:.2f},{:.2f},{:.2f}".format(*trnainfo.percentcounts())
                #print covcounts.covdict[covpos("M_dm_Heart_M5_minusAlkB",currtrna,"57")].actualbase + ":"+covcounts.covdict[covpos("M_dm_Heart_M5_minusAlkB",currtrna,"59")].actualbase
                perccounts = trnainfo.percentcounts()
                #upstreambase = covcounts.covdict[covpos("M_dm_Heart_M5_minusAlkB",currtrna,"57")]
                #downstreambase = covcounts.covdict[covpos("M_dm_Heart_M5_minusAlkB",currtrna,"59")]
                #upbases[upstreambase.actualbase] += 1
                #downbases[downstreambase.actualbase] += 1
                lista.append(perccounts[0])
                listc.append(perccounts[1])
                listg.append(perccounts[2])
                listt.append(perccounts[3])
                #aminos[currtrna.split("-")[1]] += 1
            print "A: "+twopos(np.average(lista)) +"+-"+ twopos(np.std(lista))
            print "C: "+twopos(np.average(listc)) +"+-"+ twopos(np.std(listc))
            print "G: "+twopos(np.average(listg)) +"+-"+ twopos(np.std(listg))
            print "T: "+twopos(np.average(listt)) +"+-"+ twopos(np.std(listt))
            for curramino in aminos.iterkeys():
                #print curramino+":" +str(aminos[curramino])
                pass
            for curr in nucbases:
                #print curr+" upstream:" +str(upbases[curr])
                pass
            for curr in nucbases:
                #print curr+" downstream:" +str(downbases[curr])
                pass
        sys.exit()

    
    #minreads = 50
    
    #sys.exit()
    
    mismatchpositions =  covcounts.getmismatchpos(trnainfo.gettranscripts(),sampleinfo.getsamples(),positions)
    pairgroup = None
    mismatchlocs = list(currpos for currpos in positions if currpos in mismatchpositions)

    '''
    if trnamode:
        pairgroup = gettrnasamples(mismatchlocs, trnainfo,sampleinfo)
    elif samplemode:
        pairgroup = getsamples(mismatchlocs, trnainfo,sampleinfo)
    else:
        pairgroup = getreplicates(mismatchlocs, trnainfo,sampleinfo)
    '''
    #for groupname, samplelist, trnalist, poslist in gettrnasamples(["58"], trnainfo,sampleinfo):
    #for groupname, samplelist, trnalist, poslist in getreplicates(mismatchlocs, trnainfo,sampleinfo):
    #for groupname, samplelist, trnalist, poslist in getsamples(mismatchlocs, trnainfo,sampleinfo):
    #no 56pos in results
    #outfile = open(runname+"-mismatchpos.txt","w")
    #gettrnainfo(outfile, covcounts, mismatchlocs, trnainfo,sampleinfo)
    #sys.exit()
    
    outfile = open(runname+"-repcompare.txt","w")
    pairgroup = getreplicates(mismatchlocs, trnainfo,sampleinfo)
    createtable(outfile, covcounts, pairgroup, minreads = minreads, skipmatches = False,shufflemode = False, drawmode = False)
    outfile.close()
    
    outfile = open(runname+"-repcomparedraw.txt","w")
    pairgroup = getreplicates(mismatchlocs, trnainfo,sampleinfo)
    createtable(outfile, covcounts, pairgroup, minreads = minreads, skipmatches = False,shufflemode = False, drawmode = True)
    outfile.close()
     
    
    outfile = open(runname+"-trnacompare.txt","w")
    pairgroup = gettrnasamples(mismatchlocs, trnainfo,sampleinfo)
    createtable(outfile, covcounts, pairgroup, minreads = minreads, skipmatches = skipmatches,shufflemode = False, drawmode = False)
    outfile.close()
        
    outfile = open(runname+"-trnacomparedraw.txt","w")
    pairgroup = gettrnasamples(mismatchlocs, trnainfo,sampleinfo)
    createtable(outfile, covcounts, pairgroup, minreads = minreads, skipmatches = skipmatches,shufflemode = False, drawmode = True)
    outfile.close()
    
        
    outfile = open(runname+"-samplecompare.txt","w")
    pairgroup = getsamples(mismatchlocs, trnainfo,sampleinfo)
    createtable(outfile, covcounts, pairgroup, minreads = minreads, skipmatches = False,shufflemode = False, drawmode = False)
    outfile.close()
    
        

        
    outfile = open(runname+"-samplecomparedraw.txt","w")
    pairgroup = getsamples(mismatchlocs, trnainfo,sampleinfo)
    createtable(outfile, covcounts, pairgroup, minreads = minreads, skipmatches = False,shufflemode = False, drawmode = True)
    outfile.close()
        

 
            
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
    parser.add_argument('--covfile',
                       help='coverage file')
    parser.add_argument('--trnafile',
                       help='trna file')
    parser.add_argument('--samplefile',
                       help='sample file')
    parser.add_argument('--runname',
                       help='run name')
    parser.add_argument('--minreads',
                       help='minimum reads')
    parser.add_argument('--skipperfect', action="store_true", default=False,
                       help='skip perfect matches to reference base')
    args = parser.parse_args()
    argdict = vars(args)
    main(**argdict)