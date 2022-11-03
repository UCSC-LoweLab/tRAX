#!/usr/bin/env python

import pysam
import sys
import argparse
import string
import itertools
from collections import defaultdict
import os.path
from trnasequtils import *
import numpy as np
from scipy import stats
import random
from scipy.spatial.distance import euclidean



'''


/projects/lowelab/users/holmes/pythonsource/trnatest/test/poscompare/Mouse_Brain_M4_minusAlkB_67pos_A_clust_groups.txt

/projects/lowelab/users/holmes/pythonsource/TRAX/comparetrnasets.py --trnafile=~/pythonsource/trnatest/trnadbs/mm10_new/mm10-trnatable.txt --stkfile=~/pythonsource/trnatest/trnadbs/mm10_new/mm10-trnaalign.stk

/projects/lowelab/users/holmes/pythonsource/TRAX/comparetrnasets.py --trnafile=~/pythonsource/trnatest/trnadbs/hg19/hg19-trnatable.txt --stkfile=~/pythonsource/trnatest/trnadbs/hg19/hg19-trnaalign.stk

'''

positions = list([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,'17a',18,19,20,'20a','20b',21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14','e15','e16','e17','e18','e19',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76])


testpositions = list([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76])

testpositions = set(str(curr) for curr in testpositions)

bases = ['A','U','C','G']        

def basetoprob(infreqs):
    total = float(sum(infreqs[curr] for curr in bases)) + .01
    return list(infreqs[curr]/total for curr in bases)
def getdist(fir, sec):
    #print basetoprob(fir)
    #print basetoprob(sec)
    #print hellingerdistance(basetoprob(fir),basetoprob(sec))
    return hellingerdistance(basetoprob(fir),basetoprob(sec))

sqrt2 = np.sqrt(2)
def hellingerdistance(firset, secset):    

    return euclidean(np.sqrt(firset), np.sqrt(secset))/sqrt2 


def multifasta(allseqs, filename):
    fafile = open(filename, "w")


    for seqname, seq in allseqs.items():
        fafile.write(">"+seqname+"\n")
        fafile.write(seq+"\n")
    fafile.flush()
#print >>sys.stderr, len(positions)
#this gets the tRNA numbers by the sprinzel numbering system
def gettnanums(trnaalign, margin = 0):
    trnanum = list()
    currcount = 0
    enum = 1
    gapnum = 1
    intronnum = 1
    for i in range(margin):
        trnanum.append('head'+str(margin - i))
    for i, struct in enumerate(trnaalign.consensus):
        if currcount >= len(positions):
            trnanum.append('gap'+str(gapnum))
            gapnum += 1
            currcount += 1
        elif struct in  set("+=*"):
            #special case to account for differences between loci/transcripts
            if currcount == 0 and struct == '=':
                currcount = 1
                gapnum = 1
            if positions[currcount] == 'e':
                trnanum.append('e'+str(enum))
                enum += 1
                currcount += 1
                gapnum = 1
            elif positions[currcount] == '-':
                trnanum.append(str(currcount)+'.gap'+str(gapnum))
                gapnum += 1
                currcount += 1
            else:
                trnanum.append(str(positions[currcount]))
                currcount += 1
                gapnum = 1
        else:
            #if intron

            if positions[currcount] == 38:
                trnanum.append('intron'+str(intronnum))
                intronnum += 1
            else:
                
                trnanum.append(str(currcount)+'.gap'+str(gapnum))
                gapnum += 1
    for i in range(margin):
        trnanum.append('tail'+str(i+1))
    return trnanum


def clusttrnafile(clustfile):
    clustdict = defaultdict(set)
    infile = open(clustfile)
    for currline in infile:
        fields = currline.split()
        trnaname = fields[0].strip('"')
        cluster = fields[1].strip()
        clustdict[cluster].add(trnaname)
    return clustdict



def traxtrnafile(logfoldchange, pvalfile, countfile, trnas, pvalcutoff = .05, minlogdiff = 0):
    logchange = defaultdict(lambda: defaultdict(dict))
    pvals = defaultdict(lambda: defaultdict(dict))
    counts = defaultdict(lambda: defaultdict(dict))
    maxlogdiff = 1
    headers = list()
    
    subtypes = set()
    for i, currline in enumerate(open(logfoldchange)):
        fields = currline.split()
        if i == 0:
            headers = list(curr.strip('"') for curr in fields)
            continue
        rowset = fields[0].strip('"').split("_")
        rowname = rowset[0]
        rowtype = ""
        if len(rowset) > 1:
            rowtype = rowset[1]
            subtypes.add(rowtype)
        for j, curr in enumerate(headers):
            if fields[j+1] == "NA":
                logchange[headers[j]][rowname][rowtype] = 0.0
            else:
                logchange[headers[j]][rowname][rowtype] = float(fields[j+1])
    for i, currline in enumerate(open(pvalfile)):
        fields = currline.split()
        if i == 0:
            headers = list(curr.strip('"') for curr in fields)
            continue
        rowset = fields[0].strip('"').split("_")
        rowname = rowset[0]
        rowtype = ""
        if len(rowset) > 1:
            rowtype = rowset[1]
            subtypes.add(rowtype)
        for j, curr in enumerate(headers):
            #print >>sys.stderr, headers[j]
            #print >>sys.stderr, rowname
            if fields[j+1] == "NA":
                pvals[headers[j]][rowname][rowtype] = 1.0
            else:
                pvals[headers[j]][rowname][rowtype] = float(fields[j+1])
    trnasums = defaultdict(dict)
    for i, currline in enumerate(open(countfile)):
        fields = currline.split()
        if i == 0:
            headers = list(curr.strip('"') for curr in fields)
            continue
        rowset = fields[0].strip('"').split("_")
        rowname = rowset[0]
        rowtype = ""
        #print >>sys.stderr, rowset
        if len(rowset) > 1:
            rowtype = rowset[1]
            subtypes.add(rowtype)
        #if rowname.startswith("tRNA"):
            #print >>sys.stderr, rowname
        for j, curr in enumerate(headers):
            counts[headers[j]][rowname][rowtype] = float(fields[j+1])
        trnasums[rowname][rowtype] = sum(float(curr) for curr in fields[1:])
    
    #print >>sys.stderr, trnas
    samplenames = list(pvals.keys())
    #print >>sys.stderr, samplenames
    

    
    
    for samplepair in samplenames:
        
        #print >>sys.stderr, pvals[samplepair].keys()
        for currtype in subtypes:
            clustdict = defaultdict(set)
            for currtrna in trnas:
                if currtrna in pvals[samplepair]:
                    #print >>sys.stderr, pvals[samplepair][currtrna]
                    pass
                
                if currtrna in pvals[samplepair] and currtype in pvals[samplepair][currtrna]:
                    if pvals[samplepair][currtrna][currtype] < pvalcutoff:
                        #print >>sys.stderr, "**"
                        if logchange[samplepair][currtrna][currtype] > minlogdiff:
                            clustdict["plus"].add(currtrna)
                        elif logchange[samplepair][currtrna][currtype] < -minlogdiff:
                            clustdict["minus"].add(currtrna)
                    elif trnasums[currtrna][currtype] > 60 and abs(logchange[samplepair][currtrna][currtype]) < maxlogdiff:
                        clustdict["neutral"].add(currtrna)
                        
                        #if currtype == 'SNORD92':
                        #print >>sys.stderr, pvals[samplepair][currtrna][currtype]
            if len(clustdict["plus"]) == 0 and len(clustdict["minus"]) == 0:
                continue
            yield samplepair+"_"+currtype, clustdict 
    




def traxgenefile(logfoldchange, pvalfile, countfile, genes, pvalcutoff = .05, minlogdiff = 0):
    logchange = defaultdict(lambda: defaultdict(dict))
    pvals = defaultdict(lambda: defaultdict(dict))
    counts = defaultdict(lambda: defaultdict(dict))
    
    headers = list()
    
    subtypes = set()
    for i, currline in enumerate(open(logfoldchange)):
        fields = currline.split()
        if i == 0:
            headers = list(curr.strip('"') for curr in fields)
            continue
        rowset = fields[0].strip('"').split("_")
        rowname = rowset[0]
        rowtype = ""
        if len(rowset) > 1:
            rowtype = rowset[1]
            subtypes.add(rowtype)
        for j, curr in enumerate(headers):
            if fields[j+1] == "NA":
                logchange[headers[j]][rowname] = 0.0
            else:
                logchange[headers[j]][rowname] = float(fields[j+1])
    for i, currline in enumerate(open(pvalfile)):
        fields = currline.split()
        if i == 0:
            headers = list(curr.strip('"') for curr in fields)
            continue
        rowset = fields[0].strip('"').split("_")
        rowname = rowset[0]
        rowtype = ""
        if len(rowset) > 1:
            rowtype = rowset[1]
            subtypes.add(rowtype)
        for j, curr in enumerate(headers):
            #print >>sys.stderr, headers[j]
            #print >>sys.stderr, rowname
            if fields[j+1] == "NA":
                pvals[headers[j]][rowname] = 1.0
            else:
                pvals[headers[j]][rowname] = float(fields[j+1])
    trnasums = defaultdict(dict)
    for i, currline in enumerate(open(countfile)):
        fields = currline.split()
        if i == 0:
            headers = list(curr.strip('"') for curr in fields)
            continue
        rowset = fields[0].strip('"').split("_")
        rowname = rowset[0]
        rowtype = ""
        #print >>sys.stderr, rowset
        if len(rowset) > 1:
            rowtype = rowset[1]
            subtypes.add(rowtype)
        #if rowname.startswith("tRNA"):
            #print >>sys.stderr, rowname
        for j, curr in enumerate(headers):
            counts[headers[j]][rowname][rowtype] = float(fields[j+1])
        trnasums[rowname][rowtype] = sum(float(curr) for curr in fields[1:])
    
    #print >>sys.stderr, trnas
    samplenames = list(pvals.keys())
    #print >>sys.stderr, samplenames
    

    
    
    for samplepair in samplenames:
        #print >>sys.stderr, pvals[samplepair].keys()
        #print >>sys.stderr, pvals[samplepair].keys()
        clustdict = defaultdict(set)
        for currgene in genes:
            #print >>sys.stderr, currgene.name
            #print >>sys.stderr, pvals[samplepair].keys()
            if currgene.name in pvals[samplepair]:
                #print >>sys.stderr, pvals[samplepair].keys()
                #print >>sys.stderr, currtype
                pass
            
            #print >>sys.stderr, pvals[samplepair].keys()
            if currgene.name in pvals[samplepair]:
                #print >>sys.stderr, "**||"
                if pvals[samplepair][currgene.name] < pvalcutoff:
                    
                    if logchange[samplepair][currgene.name]> minlogdiff:
                        clustdict["plus"].add(currgene.name)
                    elif logchange[samplepair][currgene.name] < -minlogdiff:
                        clustdict["minus"].add(currgene.name)
                elif trnasums[currgene.name] > 60:
                    clustdict["neutral"].add(currgene.name)
                    
                    #if currtype == 'SNORD92':
                    #print >>sys.stderr, pvals[samplepair][currtrna][currtype]
        if len(clustdict["plus"]) == 0 and len(clustdict["minus"]) == 0:
            continue
        yield samplepair, clustdict 
    

newpositions = list(str(curr) for curr in positions)        
def freqline(indict):
    bases = ['A', 'C', 'G', 'U', '-']
    return " ".join(curr+":"+str(indict[curr]) for curr in bases)
def comparesets(trnalistone, trnalisttwo,trnastk, positionnums):
    firfreqcounts = defaultdict(lambda: defaultdict(int))
    secfreqcounts = defaultdict(lambda: defaultdict(int))
    for currtrna in trnalistone:
        for i, posname in enumerate(positionnums):
            if posname in testpositions:
                firfreqcounts[posname][trnastk.aligns[currtrna][i]] += 1
    #print firfreqcounts["34"]
            
    for currtrna in trnalisttwo:
        for i, posname in enumerate(positionnums):
            if posname in testpositions:
                secfreqcounts[posname][trnastk.aligns[currtrna][i]] += 1
    #print secfreqcounts["34"]
    
    #print getdist(firfreqcounts["34"],secfreqcounts["34"])
    #sys.exit()
    baseprobs = dict()   
    for currpos in newpositions:
        currdist = getdist(firfreqcounts[currpos],secfreqcounts[currpos])
        baseprobs[currpos] = currdist
    #print baseprobs["34"]
    
    #sys.exit()
    #print ",".join(trnalistone)
    #print ",".join(trnalisttwo)
    basediffs = sorted(newpositions, key = lambda x: baseprobs[x], reverse = True)
    for i in range(0,5):
        print(basediffs[i])
        print("\t"+ freqline(firfreqcounts[basediffs[i]]))
        print("\t"+ freqline(secfreqcounts[basediffs[i]]))
        
        
def printseqs(pairname, positivelist, negativelist, neutrallist,genomefile, filterlabel = None):
    filtername = ""
    if filterlabel is not None:
        filtername = "_"+filterlabel
    
    genome = fastaindex( genomefile, genomefile +".fai") 
    posseqs = genome.getseqs(positivelist)
    negseqs = genome.getseqs(negativelist)
    neutseqs = genome.getseqs(neutrallist)
    multifasta(posseqs, pairname + filtername+ "_positive-seqs.fa")
    multifasta(negseqs, pairname + filtername+ "_negative-seqs.fa")
    multifasta(neutseqs,pairname + filtername+ "_neutral-seqs.fa")
        
def getgenetypes(genetypename):
    genetypes = dict()
    genetypefile = open(genetypename, "r")
    for currline in genetypefile:
        fields = currline.split("\t")
        genetypes[fields[0]] = fields[1].rstrip()
    return genetypes
        
def main(**argdict):
    edgemargin = 0
    if "stkfile" in argdict and argdict["stkfile"] is not None:
        trnastk = list(readrnastk(open(os.path.expanduser(argdict["stkfile"]), "r")))[0]
        positionnums = gettnanums(trnastk, margin = edgemargin)
        trnastk = trnastk.addmargin(edgemargin)
        
        #print "\t".join(trnastk.aligns['tRNA-Ala-CGC-1'])
        #print "\t".join(positionnums)
        
        trnainfo  = transcriptfile(os.path.expanduser(argdict["trnafile"]))
        baseinfo = defaultdict(dict)
        
        #print positionnums
        for currtrna in trnainfo.gettranscripts():
            for i, posname in enumerate(positionnums):
                baseinfo[currtrna][posname] = trnastk.aligns[currtrna][i]
        #print baseinfo[currtrna]

    
        trnalistone = trnainfo.getaminotranscripts("Arg")
        trnalisttwo = trnainfo.getaminotranscripts("Glu")
    
    '''
    clustfile = "/projects/lowelab/users/holmes/pythonsource/trnatest/test/poscompare/Mouse_Brain_M4_minusAlkB_67pos_A_clust_groups.txt"
    trnasets = clusttrnafile(clustfile)
    for currcluster in trnasets.keys():
        print currcluster
        print "**"
        noncluster = set(itertools.chain.from_iterable(trnasets[curr] for curr in trnasets.keys() if curr != currcluster))
        comparesets(trnasets[currcluster], noncluster,trnastk)
    '''
    logfoldchange = '/projects/sharma/users/holmes/sharmasra/sharma2017samples/sharma2017samples-logvals.txt'
    pvalfile = '/projects/sharma/users/holmes/sharmasra/sharma2017samples/sharma2017samples-padjs.txt'
    countfile = '/projects/sharma/users/holmes/sharmasra/sharma2017samples/sharma2017samples-readcounts.txt'
    genetypefile = '/projects/sharma/users/holmes/sharmasra/sharma2017samples/sharma2017samples-genetypes.txt'
    ensemblgtf = "/projects/lowelab/users/holmes/pythonsource/trnatest/trnadbs/mm10/mm10.gtf.gz"
    genomefile = "/projects/lowelab/users/holmes/pythonsource/trnatest/trnadbs/mm10/mm10-tRNAgenome.fa"
    
    
    #logfoldchange = '/projects/lowelab/users/holmes/pythonsource/trnatest/kitcomp/armseqhumannew2/armseqhumannew2-logvals.txt'
    #pvalfile =      '/projects/lowelab/users/holmes/pythonsource/trnatest/kitcomp/armseqhumannew2/armseqhumannew2-padjs.txt'
    #countfile =     '/projects/lowelab/users/holmes/pythonsource/trnatest/kitcomp/armseqhumannew2/armseqhumannew2-readcounts.txt'
     
    
    genetypes = getgenetypes(genetypefile)
    genefeats = list(readgtf(ensemblgtf, filterpsuedo = False, replacename = True))
    featnames =  {curr.name: curr for curr in genefeats} 
    
    #trnatranscripts = trnainfo.gettranscripts() 
    
    pvalcutoff = .005
    
    #print "**"
    filtertype = "miRNA"
    filtertype = None
    flank = 50
    for currpair, currtrnaset in traxgenefile(logfoldchange, pvalfile, countfile,genefeats, pvalcutoff = pvalcutoff):
        clustpair = list(currtrnaset.keys())
        #print "**"
        #print genetypes.values()
        #print str(len(currtrnaset["plus"]))+"+"+str(len(currtrnaset["minus"]))+"+"+str(len(currtrnaset["neutral"]))+"/"+str(len(genefeats))
        #print currtrnaset["neutral"]
        posgenes = {curr: featnames[curr].addmargin(0) for curr in currtrnaset["plus"] if filtertype is None or genetypes[curr] == filtertype }
        minusgenes = {curr: featnames[curr].addmargin(0) for curr in currtrnaset["minus"] if filtertype is None or genetypes[curr] == filtertype }
        neutralgenes = {curr: featnames[curr].addmargin(0) for curr in currtrnaset["neutral"] if filtertype is None or genetypes[curr] == filtertype }
        printseqs(currpair, posgenes,minusgenes, neutralgenes, genomefile, filterlabel = filtertype) 
        
        allgenes = list(currtrnaset["plus"] | currtrnaset["minus"] | currtrnaset["neutral"])
        allgenes = list(curr for curr in allgenes if genetypes[curr] == filtertype)
        random.shuffle(allgenes)
        
        poslength = len(list(posgenes.keys()))
        neglength = len(list(minusgenes.keys()))
        posrand = allgenes[0:poslength]
        negrand = allgenes[poslength: poslength + neglength]
        neutrand = allgenes[poslength + neglength:]
        #print >>sys.stderr, "***"
        #print >>sys.stderr, str(poslength)
        #print >>sys.stderr, str(len(posrand))
        
        posgenes = {curr: featnames[curr].addmargin(0) for curr in posrand if filtertype is None or genetypes[curr] == filtertype }
        minusgenes = {curr: featnames[curr].addmargin(0) for curr in negrand if filtertype is None or genetypes[curr] == filtertype }
        neutralgenes = {curr: featnames[curr].addmargin(0) for curr in neutrand if filtertype is None or genetypes[curr] == filtertype }
        printseqs(currpair+"_rand", posgenes,minusgenes, neutralgenes, genomefile, filterlabel = filtertype) 

    ''' 
    for currpair, currtrnaset in traxtrnafile(logfoldchange, pvalfile, countfile, trnatranscripts, pvalcutoff = pvalcutoff):
        clustpair = list(currtrnaset.keys())
        print currpair
        print str(len(currtrnaset["plus"]))+"+"+str(len(currtrnaset["minus"]))+"+"+str(len(currtrnaset["neutral"]))+"/"+str(len(trnatranscripts))
        #comparesets(currtrnaset["plus"], currtrnaset["minus"],trnastk, positionnums)
        print "Plus exclusive"
        comparesets(currtrnaset["plus"] , currtrnaset["minus"] | currtrnaset["neutral"],trnastk, positionnums)
        print "Minus exclusive"
        comparesets(currtrnaset["minus"], currtrnaset["plus"] | currtrnaset["neutral"],trnastk, positionnums)
    '''                                                                              
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
    parser.add_argument('--trnafile', 
                       help='trna file')
    parser.add_argument('--stkfile',
                       help='Stockholm file')

    
    
    '''
    Perform check on sizefactor file to ensure it has all samples
    '''
    args = parser.parse_args()
    argdict = vars(args)
    main(**argdict)