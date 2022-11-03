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
from multiprocessing import Process, Queue, Pool


bam_match = 0
bam_cins = 1
bam_cdel = 2

#pseudocount = 30


def cigarreflength(cigar):
    return sum(curr[1] for curr in cigar if curr[0] in set([0,bam_cdel]))
    
def cigarreadlength(cigar):
    return sum(curr[1] for curr in cigar if curr[0] in set([0,bam_cins]))
    
def cigarrefcoverage(cigar):
    nextsum = 1
    for curr in cigar:
        if curr[0] == bam_cins:
            nextsum += curr[1]
            pass
        elif curr[0] == bam_cdel:
            for i in range(curr[1]):
                yield 0
        elif curr[0] == bam_match:
            for i in range(curr[1]):
                yield nextsum
                nextsum = 1
                
                
    


gapchars = set("-._~")
class readcoverage:
    def __init__(self, region):
        self.region = region
        self.samplereads = 0
        self.coverage = list()
        self.length = region.length()
        self.totalreads = 0
        for i in range(0,region.length()):
            self.coverage.append(0)

                
    def coveragelist(self):
        return self.coverage
    def coveragealign(self, alignment, gapoutput = None,sizefactor = 1):          
        if len(self.coverage) != len(alignment.translate( str.maketrans('', '', str(gapchars)))):
            print("Alignment length does not match bed length "+str(len(self.coverage))+" "+str(len(string.translate(alignment, None, str(gapchars)))), file=sys.stderr) 
            print(self.coverage, file=sys.stderr)
            print(string.translate(alignment, None, str(gapchars)), file=sys.stderr)
            sys.exit(1)
        i = 0
        lastcoverage = None
        for curr in alignment:
            #print >>sys.stderr, curr
            if curr in gapchars:
                if lastcoverage is None:
                    yield gapoutput
                    
                else:
                    yield lastcoverage
            else:
                lastcoverage = self.coverage[i]/sizefactor
                yield self.coverage[i]/sizefactor
                i += 1
    def addread(self, read):
        self.totalreads += 1
        if self.region.strand == "+":
            start = max([0, read.start - self.region.start])
            end = min([self.length, read.end - self.region.start])
            
        else:
            start = max([0, self.region.end - read.end])
            end = min([self.length, self.region.end - read.start])
        
        for currpos in range(self.length):
            if start <= currpos <= end - 1:
                self.coverage[currpos] += 1
    def addbase(self, base):
        #self.totalreads += 1
 
        posbase = base - self.region.start
        if 0 <= posbase <= len(self.coverage) - 1:
            self.coverage[posbase] += 1
        else:
            pass
        
        
def nasum(operands, naval = "NA"):
    if sum("NA"== curr for curr in operands) == len(operands):
        return "NA"
    elif not any("NA"== curr for curr in operands):
        return sum(operands)
    else:
        print("Trying to add incompatible alignments", file=sys.stderr)
        sys.exit(1)
    #return ",".join(str(curr) for curr in operands)
def sumsamples(coverage,sampledata, repname, currfeat, sizefactors = defaultdict(lambda: 1)):
    return (nasum(curr) for curr in zip(*(allcoverages[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name],sizefactor = sizefactor[currsample]) for currsample in sampledata.getrepsamples(repname))))
    
count = 0







def readtrnanums(numfile, margin = 0):
    for i in range(margin):
        yield 'head'+str(margin - i)
    for currline in open(numfile):
        fields = currline.rstrip().split("\t")
        yield fields[3]
    for i in range(margin):
        yield 'tail'+str(i+1)


eukpositions = list([-1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,'17a',18,19,20,'20a','20b',21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14','e15','e16','e17','e18','e19',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76])
archpositions = list([-1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,'17a',18,19,20,'20a','20b',21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76])


bactpositions = list([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14','e15','e16','e17',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76])
mitopositions = list([-1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14','e15','e16','e17',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76])




#print >>sys.stderr, len(positions)
#this gets the tRNA numbers by the sprinzel numbering system


def gettnanums(trnaalign, margin = 0, orgtype = "euk"):
    trnanum = list()
    currcount = 0
    enum = 1
    gapnum = 1
    intronnum = 1
    positions = eukpositions 
    if orgtype == "arch":
    	positions = archpositions
    elif orgtype == "mito":
        positions = mitopositions
    elif orgtype == "bact":
        positions = bactpositions 
    for i in range(margin):
        trnanum.append('head'+str(margin - i))
    for i, struct in enumerate(trnaalign.consensus):
        if currcount >= len(positions):
            trnanum.append('gap'+str(gapnum))
            gapnum += 1
            currcount += 1
        elif struct in  set("+=*"):
            #special case to account for differences between loci/transcripts
            if currcount == 0 and struct == '=' and orgtype != "bact":
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
    
class coverageinfo:
    def __init__(self, readcounts, allcoverage,readstarts, readends, multaminocoverages, multaccoverages, multtrnacoverages,uniquecoverages, uniquegenomecoverages,multigenomecoverages, readmismatches,adeninemismatches,thyminemismatches,cytosinemismatches, guanosinemismatches, readskips, trimcoverage = None, trimmismatches = None  ):
        self.readcounts = readcounts
        self.allcoverages = allcoverage
        self.readstarts = readstarts
        self.readends = readends
        
        self.multaminocoverages = multaminocoverages
        self.multaccoverages = multaccoverages
        self.multtrnacoverages = multtrnacoverages
        self.uniquecoverages = uniquecoverages
        self.uniquegenomecoverages = uniquegenomecoverages
        self.multigenomecoverages = multigenomecoverages
        
        
        self.trimcoverage = trimcoverage
        self.trimmismatches = trimmismatches
        
        
        self.readmismatches = readmismatches
        self.adeninemismatches = adeninemismatches
        self.thyminemismatches = thyminemismatches
        self.cytosinemismatches = cytosinemismatches
        self.guanosinemismatches = guanosinemismatches    
        self.readskips = readskips
        
        
        
class locicoverageinfo:
    def __init__(self, readcounts, allcoverage):
        self.readcounts = readcounts
        self.allcoverages = allcoverage
def getlocicoverage(currsample, sampledata, trnaloci,maxmismatches = None, minextend = None): 
    currbam = sampledata.getbam(currsample)
    allcoverages = dict()
    readcounts = dict()
    
    try:
        #print >>sys.stderr, currbam
        if not os.path.isfile(currbam+".bai"):
            pysam.index(""+currbam)
        bamfile = pysam.Samfile(""+currbam, "rb" )  
    except IOError as xxx_todo_changeme1:
        ( strerror) = xxx_todo_changeme1
        print(strerror, file=sys.stderr)
        sys.exit()
        
    for i, currfeat in enumerate(trnaloci):
        allcoverages[currfeat.name] = readcoverage(trnaloci[i])
        readcounts[currfeat.name] = 0
        for currread in getbam(bamfile, trnaloci[i]):
            if maxmismatches is not None and currread.getmismatches() > maxmismatches:
                continue
            #print >>sys.stderr, "||**||"+str(currread.getmismatches())
            if trnaloci[i].coverage(currread) > 10:

                if minextend is not None and not (currread.start + minextend <= trnaloci[i].start or currread.end - minextend >= trnaloci[i].end):
                    continue

                readcounts[trnaloci[i].name] += 1
                allcoverages[trnaloci[i].name].addread(currread)
                
    return locicoverageinfo( readcounts, allcoverages)

def getsamplecoverage(currsample, sampledata, trnalist, trnaseqs,maxmismatches = None, minextend = None, removestart = True, uniqueonly = False): 
    
    currbam = sampledata.getbam(currsample)
    allcoverages = dict()
    multaminocoverages = dict()
    multaccoverages = dict()
    multtrnacoverages = dict()
    uniquecoverages = dict()
    uniquegenomecoverages = dict()
    multigenomecoverages = dict()
    #print >>sys.stderr, trnalist
    readmismatches = dict()
    
    adeninemismatches = dict()
    thyminemismatches = dict()
    cytosinemismatches = dict()
    guanosinemismatches = dict()
    readstarts = dict()
    readends = dict()
    readskips = dict()      
    
    trimreadcoverage =  dict()
    trimreadmismatches =  dict()
    
    readcounts = dict()
    
    try:
        #print >>sys.stderr, currbam
        if not os.path.isfile(currbam+".bai"):
            pysam.index(""+currbam)
        bamfile = pysam.Samfile(""+currbam, "rb" )  
    except IOError as xxx_todo_changeme2:
        ( strerror) = xxx_todo_changeme2
        print(strerror, file=sys.stderr)
        sys.exit()
        
    for i, currfeat in enumerate(trnalist):

        allcoverages[currfeat.name] = readcoverage(trnalist[i])
        multaminocoverages[currfeat.name] = readcoverage(trnalist[i])
        multaccoverages[currfeat.name] = readcoverage(trnalist[i])
        multtrnacoverages[currfeat.name] = readcoverage(trnalist[i])
        uniquecoverages[currfeat.name] = readcoverage(trnalist[i])
        uniquegenomecoverages[currfeat.name] = readcoverage(trnalist[i])
        multigenomecoverages[currfeat.name] = readcoverage(trnalist[i])
        
        readstarts[currfeat.name] = readcoverage(trnalist[i])
        readends[currfeat.name] = readcoverage(trnalist[i])
        readmismatches[currfeat.name] = readcoverage(trnalist[i])
        adeninemismatches[currfeat.name] =   readcoverage(trnalist[i])
        thyminemismatches[currfeat.name] =   readcoverage(trnalist[i])
        cytosinemismatches[currfeat.name] =  readcoverage(trnalist[i])
        guanosinemismatches[currfeat.name] = readcoverage(trnalist[i])
        readskips[currfeat.name] = readcoverage(trnalist[i])
        
        trimreadcoverage[currfeat.name] =  readcoverage(trnalist[i])
        trimreadmismatches[currfeat.name] =  readcoverage(trnalist[i])
        readcounts[currfeat.name] = 0

        #print >>sys.stderr, trnalist[i]
        for currread in getbam(bamfile, trnalist[i]):
            

            if maxmismatches is not None and currread.getmismatches() > maxmismatches:
                continue
            #print >>sys.stderr, "||**||"+str(currread.getmismatches())
            if trnalist[i].coverage(currread) > 10:

                if minextend is not None and not (currread.start + minextend <= trnalist[i].start or currread.end - minextend >= trnalist[i].end):
                    continue
                if uniqueonly and not currread.issinglemapped():
                    continue

                readstart = currread.getfirst(1)
                readend = currread.getlast(1)
                readcounts[trnalist[i].name] += 1
                allcoverages[trnalist[i].name].addread(currread)
                readstarts[trnalist[i].name].addread(readstart)
                readends[trnalist[i].name].addread(readend )

                if not currread.isuniqueaminomapping():
                    multaminocoverages[trnalist[i].name].addread(currread)
                elif not currread.isuniqueacmapping():
                    multaccoverages[trnalist[i].name].addread(currread)
                elif not currread.isuniquetrnamapping():
                    multtrnacoverages[trnalist[i].name].addread(currread)
                else:
                    uniquecoverages[trnalist[i].name].addread(currread)
                if currread.issinglemapped():
                    uniquegenomecoverages[trnalist[i].name].addread(currread)
                else:
                    multigenomecoverages[trnalist[i].name].addread(currread)
                
                currseq = currread.getseq()
                #allcoverages[currsample][trnaname].addread(currread)
                readcov = list(cigarrefcoverage(currread.getcigar()))
                trnaname = trnalist[i].name
                
                alignseq = "".join(currseq[sum(readcov[0:i])] if readcov[i] > 0 else "-" for i in range(cigarreflength(currread.getcigar())))
                refseq = trnaseqs[currfeat.name][currread.start:currread.start+ cigarreflength(currread.getcigar())]
                #if currread.name == "NB501427:156:H2F7MAFXY:3:21609:5222:18520":
                #if cigarreflength(currread.getcigar()) < cigarreadlength(currread.getcigar()):
                #    print >>sys.stderr, currread.name
                #    print >>sys.stderr, currread.getcigar()
                #    print >>sys.stderr, "".join(str(curr) for curr in readcov)
                #    print >>sys.stderr, alignseq 
                #    print >>sys.stderr,refseq
                #    print >>sys.stderr, currseq       
                for currpos in range(cigarreflength(currread.getcigar())): #30
                    currbase = alignseq[currpos]
                    refbase = refseq[currpos]
                    if refbase not in gapchars:
                        if currbase == "-":
                            readskips[trnaname].addbase(currread.start + currpos)
                        if refbase != currbase:
                            #if (currpos + currread.start) - readmismatches[trnaname].region.start < 0:
                            #    print >>sys.stderr, "before start: "+str(currpos)+"+"+str(currread.start) +"-"+str(readmismatches[trnaname].region.start)
                            #    #base - self.region.start
                            readmismatches[trnaname].addbase(currread.start + currpos)
                        if currpos > 3:
                            trimreadcoverage[trnaname].addbase(currread.start + currpos)
                            if refbase != currbase:
                                trimreadmismatches[trnaname].addbase(currread.start + currpos)
                        if currpos > 3 and removestart:
                            if currbase == "A":
                                adeninemismatches[trnaname].addbase(currread.start + currpos)
                            elif currbase == "T":
                                thyminemismatches[trnaname].addbase(currread.start + currpos)
                            elif currbase == "C":
                                cytosinemismatches[trnaname].addbase(currread.start + currpos)
                            elif currbase == "G":
                                guanosinemismatches[trnaname].addbase(currread.start + currpos)

    return coverageinfo( readcounts, allcoverages,readstarts, readends,multaminocoverages, multaccoverages, multtrnacoverages,uniquecoverages, uniquegenomecoverages,multigenomecoverages, readmismatches,adeninemismatches,thyminemismatches,cytosinemismatches, guanosinemismatches,readskips,trimmismatches = trimreadmismatches, trimcoverage = trimreadcoverage  )

def transcriptcoverage(samplecoverages, mismatchreport, trnalist,sampledata,sizefactor, mincoverage, trnastk, positionnums, skipgaps = True):

    #print >>mismatchreport, "\t".join(["Feature","Sample","position","coverage","readstarts","readends","uniquecoverage","multitrnacoverage","multianticodoncoverage","multiaminocoverage","tRNAreadstotal","actualbase","mismatchedbases","deletedbases","adenines","thymines","cytosines","guanines","deletions"])
    #print >>sys.stderr,mismatchreport
    #print >>sys.stderr,"||***"
    samples = sampledata.getsamples()
    for currfeat in trnalist:
        #print >>sys.stderr, samplecoverages[list(samples)[0]].allcoverages.keys()
        totalreads = sum(samplecoverages[currsample].allcoverages[currfeat.name].totalreads for currsample in samples)
        ambigreads = sum(samplecoverages[currsample].multaminocoverages[currfeat.name].totalreads for currsample in samples)
        if totalreads - ambigreads < mincoverage:
            continue
        reportpositions = set()  
        for currsample in samples:
            #readcounts = samplecoverages[currsample].readcounts
            
            #print >>sys.stderr, ",".join(list(str(curr) if curr is not None else 0 for curr in samplecoverages[currsample].readmismatches[currfeat.name].coveragealign(trnastk.aligns[currfeat.name])))
            #covcounts = list(curr/(1.*samplecoverages[currsample].readcounts[currfeat.name]) if curr is not None else 0 for curr in samplecoverages[currsample].readmismatches[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            
            covcounts =  list(curr/sizefactor[currsample] if curr is not None else 0 for curr in samplecoverages[currsample].allcoverages[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))

            mismatches =  list(curr/sizefactor[currsample] if curr is not None else 0 for curr in samplecoverages[currsample].readmismatches[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            deletions =  list(curr/sizefactor[currsample] if curr is not None else 0 for curr in samplecoverages[currsample].readskips[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))


            uniquecounts =  list(curr/sizefactor[currsample] if curr is not None else 0 for curr in samplecoverages[currsample].uniquecoverages[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            multitrna =  list(curr/sizefactor[currsample] if curr is not None else 0 for curr in samplecoverages[currsample].multtrnacoverages[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            multaccounts =  list(curr/sizefactor[currsample] if curr is not None else 0 for curr in samplecoverages[currsample].multaccoverages[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            multaminocounts =  list(curr/sizefactor[currsample] if curr is not None else 0 for curr in samplecoverages[currsample].multaminocoverages[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))



            allstarts  = list(curr/sizefactor[currsample] if curr is not None else 0 for curr in samplecoverages[currsample].readstarts[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            allends = list(curr/sizefactor[currsample] if curr is not None else 0 for curr in samplecoverages[currsample].readends[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            allcovcount  = list(curr/sizefactor[currsample] if curr is not None else 0 for curr in samplecoverages[currsample].allcoverages[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            adeninecount  = list(curr if curr is not None else 0 for curr in samplecoverages[currsample].adeninemismatches[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            thyminecount = list(curr if curr is not None else 0 for curr in samplecoverages[currsample].thyminemismatches[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            cytosinecount = list(curr if curr is not None else 0 for curr in samplecoverages[currsample].cytosinemismatches[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            guanosinecount  = list(curr if curr is not None else 0 for curr in samplecoverages[currsample].guanosinemismatches[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            readskipcount  = list(curr if curr is not None else 0 for curr in samplecoverages[currsample].readskips[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))

            for i, currcount in enumerate(allcovcount):
                #removing all non-canonical positions to save space/memory
                if skipgaps and "gap" in positionnums[i]:
                    continue
                mismatchthreshold = .1
                realbase = trnastk.aligns[currfeat.name][i].upper()
                if realbase in gapchars:
                    
                    realbase = "-"
                if realbase == "U":
                    realbase = "T"
                print("\t".join([currfeat.name,currsample,str(positionnums[i]),str(covcounts[i]),str(allstarts[i]),str(allends[i]),str(uniquecounts[i]),str(multitrna[i]),str(multaccounts[i]),str(multaminocounts[i]),str(1.*samplecoverages[currsample].readcounts[currfeat.name]/sizefactor[currsample]),realbase,str(mismatches[i]),str(deletions[i]),str(adeninecount[i]),str(thyminecount[i]),str(cytosinecount[i]),str(guanosinecount[i]), str(readskipcount[i])]), file=mismatchreport)
    #sys.exit(1)        

def locuscoverage(locicoverages, locicoveragetable, locilist,sampledata,sizefactor, mincoverage, locistk, locipositionnums, skipgaps = True):

    #print >>locicoveragetable, "\t".join(["tRNA_name","sample","position","coverage", "total"])
    samples = sampledata.getsamples()
    for currfeat in locilist:

        totalreads = sum(locicoverages[currsample].allcoverages[currfeat.name].totalreads for currsample in samples)
        if totalreads < mincoverage:
            continue
        reportpositions = set()  
        for currsample in samples:
            
            allcovcount  = list(curr/sizefactor[currsample] if curr is not None else 0 for curr in locicoverages[currsample].allcoverages[currfeat.name].coveragealign(locistk.aligns[currfeat.name]))

            for i, currcount in enumerate(allcovcount):
                mismatchthreshold = .1
                if skipgaps and "gap" in locipositionnums[i]:
                    continue

                print("\t".join([currfeat.name,currsample,str(locipositionnums[i]),str(currcount),str(locicoverages[currsample].readcounts[currfeat.name])]), file=locicoveragetable)
    #sys.exit(1)  


class positioninfo:
    def __init__(self, trnaname, positionname):
        self.trnaname = trnaname
        self.positionname = positionname
        self.sampletotal = dict()
        self.samplemism = dict()
        self.samplegaps = dict()
        self.samplebase = dict()
        self.sampleacounts = dict()
        self.sampletcounts = dict()
        self.sampleccounts = dict()
        self.samplegcounts = dict()   
        self.allsamples = set()
        
        self.trimreadcount = dict()
        self.trimmismatchcount = dict()   
        
    def addsample(self, currsample, totalcounts, mismatches, gaps, base, acounts, tcounts, ccounts,gcounts, trimreadcount, trimmismatchcount):
        self.allsamples.add(currsample)
        self.sampletotal[currsample]   = totalcounts
        self.samplemism[currsample]    = mismatches
        self.samplegaps[currsample]    = gaps
        self.samplebase[currsample]    = base
        self.sampleacounts[currsample] = acounts
        self.sampletcounts[currsample] = tcounts
        self.sampleccounts[currsample] = ccounts
        self.samplegcounts[currsample] = gcounts
        
        self.trimreadcount[currsample] = trimreadcount           
        self.trimmismatchcount[currsample] = trimmismatchcount
    '''
    get the median of the mismatchs while incorporating the 30 mismatch cutoff
    '''
    def getsigdiffs(self, repdict):
        
        for firsample, secsample in itertools.combinations(list(repdict.keys()),2):
            
            #firmismedian = (self.samplemism[firsample]/(1.*self.sampletotal[firsample] for currrep in repdict    
            if self.sampletotal[firsample] < 40 or self.sampletotal[secsample] < 40:
                continue
            firmism = self.samplemism[firsample]/(1.*self.sampletotal[firsample])
            secmism = self.samplemism[secsample]/(1.*self.sampletotal[secsample])
            
            if abs(firmism - secmism) > .2:
                yield tuple([self.trnaname,self.positionname,firsample,secsample,self.samplemism[firsample],self.sampletotal[firsample],self.samplemism[secsample],self.sampletotal[secsample]])
    def getsamplepair(self, firsample, secsample):
        #firmism = self.samplemism[firsample]/(1.*self.sampletotal[firsample])
        #secmism = self.samplemism[secsample]/(1.*self.sampletotal[secsample])
        if firsample in self.sampletotal and secsample in self.sampletotal:
            yield tuple([self.trnaname,self.positionname,firsample,secsample,self.samplemism[firsample],self.sampletotal[firsample],self.samplemism[secsample],self.sampletotal[secsample],        self.trimmismatchcount[firsample],self.trimreadcount[firsample],self.trimmismatchcount[secsample],self.trimreadcount[secsample]                        ])
            
        


def getdiffs(samplecoverages, mismatchreport, trnalist,sampledata, sizefactor,mincoverage, trnastk, positionnums):

    #print >>mismatchreport, "\t".join(["Feature","Sample","position","coverage","ends","uniquecoverage","multitrnacoverage","multianticodoncoverage","multiaminocoverage","tRNAreadstotal","actualbase","mismatchedbases","adenines","thymines","cytosines","guanines","deletions"])
    samples = sampledata.getsamples()
    posinfo = dict()
    for currfeat in trnalist:
        posinfo[currfeat.name] = dict()
        for currpos in positionnums:
            posinfo[currfeat.name][currpos] = positioninfo(currfeat.name, currpos) 
    
    for currfeat in trnalist:
        #print >>sys.stderr, samplecoverages[list(samples)[0]].allcoverages.keys()
        totalreads = sum(samplecoverages[currsample].allcoverages[currfeat.name].totalreads for currsample in samples)
        if totalreads < mincoverage:
            continue
        reportpositions = set()
        for currsample in samples:
            covcounts =  list(curr if curr is not None else 0 for curr in samplecoverages[currsample].allcoverages[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))

            mismatches =  list(curr if curr is not None else 0 for curr in samplecoverages[currsample].readmismatches[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))

            uniquecounts =  list(curr if curr is not None else 0 for curr in samplecoverages[currsample].uniquecoverages[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            multitrna =  list(curr if curr is not None else 0 for curr in samplecoverages[currsample].multtrnacoverages[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            multaccounts =  list(curr if curr is not None else 0 for curr in samplecoverages[currsample].multaccoverages[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            multaminocounts =  list(curr if curr is not None else 0 for curr in samplecoverages[currsample].multaminocoverages[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            
            allstarts  = list(curr if curr is not None else 0 for curr in samplecoverages[currsample].readstarts[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            allcovcount  = list(curr if curr is not None else 0 for curr in samplecoverages[currsample].allcoverages[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            adeninecount  = list(curr if curr is not None else 0 for curr in samplecoverages[currsample].adeninemismatches[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            thyminecount = list(curr if curr is not None else 0 for curr in samplecoverages[currsample].thyminemismatches[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            cytosinecount = list(curr if curr is not None else 0 for curr in samplecoverages[currsample].cytosinemismatches[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            guanosinecount  = list(curr if curr is not None else 0 for curr in samplecoverages[currsample].guanosinemismatches[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            readskipcount  = list(curr if curr is not None else 0 for curr in samplecoverages[currsample].readskips[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))

            trimcovcount  = list(curr if curr is not None else 0 for curr in samplecoverages[currsample].trimcoverage[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            trimmismatchcount  = list(curr if curr is not None else 0 for curr in samplecoverages[currsample].trimmismatches[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))

            for i, currcount in enumerate(allcovcount):
                mismatchthreshold = .1
                
                
                currposition = positionnums[i]
                actualbase = trnastk.aligns[currfeat.name][i].upper()
                if actualbase == "U":
                    actualbaset = "T"
                if actualbase in gapchars:
                    continue
                posinfo[currfeat.name][currposition].addsample(currsample, covcounts[i],mismatches[i],readskipcount[i],actualbase,adeninecount[i],thyminecount[i],cytosinecount[i],guanosinecount[i], trimcovcount[i], trimmismatchcount[i])
                #print >>mismatchreport, "\t".join([currfeat.name,currsample,str(positionnums[i]),str(covcounts[i]),str(allstarts[i]),str(uniquecounts[i]),str(multitrna[i]),str(multaccounts[i]),str(multaminocounts[i]),str(1.*samplecoverages[currsample].readcounts[currfeat.name]),trnastk.aligns[currfeat.name][i],str(mismatches[i]),str(adeninecount[i]),str(thyminecount[i]),str(cytosinecount[i]),str(guanosinecount[i]), str(readskipcount[i])])
    print("\t".join(["pos","firsample","secsample","firmismatches","firtotal","secmismatches","sectotal","firmismatchestrim","firtotaltrim","secmismatchestrim","sectotaltrim"]), file=mismatchreport)


                
                 
 
    
def genomeprint(samplecoverages, uniquegenome, trnalist,sampledata, mincoverage):
    covfiles = {uniquegenome + '-uniquegenomecoverages.txt':uniquegenomecoverages,uniquegenome + '-multgenomecoverages.txt':multigenomecoverages}
    
    for filename, currcoverage in covfiles.items():
        covfile = open(filename, "w")
        print("Feature"+"\t"+"Sample"+"\t"+"\t".join(positionnums), file=covfile)
            
        for currfeat in trnalist:
            totalreads = sum(allcoverages[currsample][currfeat.name].totalreads for currsample in samples)
            if totalreads < mincoverage:
                continue        
            replicates = sampledata.allreplicates()
            for currrep in replicates:
                print(currfeat.name+"\t"+currrep+"\t"+"\t".join(str(curr) for curr in sumsamples(currcoverage,sampledata,currrep,currfeat,sizefactor)), file=covfile)

        covfile.close()
        
def uniqcoverage(samplecoverages, uniquename, trnalist,sampledata, mincoverage):
    covfiles = {uniquename + '-uniquecoverages.txt':uniquecoverages,uniquename + '-multaccoverages.txt':multaccoverages,uniquename + '-multtrnacoverages.txt':multtrnacoverages,uniquename + '-multaminocoverages.txt':multaminocoverages}
    
    for filename, currcoverage in covfiles.items():
        covfile = open(filename, "w")
        print("Feature"+"\t"+"Sample"+"\t"+"\t".join(positionnums), file=covfile)
            
        for currfeat in trnalist:
            totalreads = sum(allcoverages[currsample][currfeat.name].totalreads for currsample in samples)
            if totalreads < mincoverage:
                continue
        
            replicates = sampledata.allreplicates()
            for currrep in replicates:
                print(currfeat.name+"\t"+currrep+"\t"+"\t".join(str(curr) for curr in sumsamples(currcoverage,sampledata,currrep,currfeat,sizefactor)), file=covfile)
                
                
        covfile.close()
def makecoveragepool(args):
    return getsamplecoverage(*args[0], **args[1])
def makelocicoveragepool(args):
    return getlocicoverage(*args[0], **args[1])
def compressargs( *args, **kwargs):
    return tuple([args, kwargs])
def testmain(**argdict):
    #print >>sys.stderr, argdict
    argdict = defaultdict(lambda: None, argdict)
    if "edgemargin" not in  argdict:                    
        edgemargin = 0
    else:
        edgemargin = int(argdict["edgemargin"])
    #currently crashes if set to zero
    if "mincoverage" not in  argdict or argdict["mincoverage"] is None:
        mincoverage = 30
    else:
        mincoverage = int(argdict["mincoverage"])  
    
    if "bamdir" not in argdict:
        bamdir = "./"
    bamdir = argdict["bamdir"]
    sampledata = samplefile(argdict["samplefile"], bamdir = bamdir)
    trnafasta = argdict["trnafasta"]
    
    trnaseqs = fastadict(trnafasta)
    for currname in list(trnaseqs.keys()):
        trnaseqs[currname] = ("N"*edgemargin)+trnaseqs[currname]+("N"*edgemargin)
        
        
    maxmismatches = argdict["maxmismatches"]
    #uniquename = argdict["uniquename"]
    cores = argdict["cores"]
    threadmode = True
    if cores == 1:
        threadmode = False
        
    #uniquegenome = argdict["uniquegenome"]
    trnastk = list(readrnastk(open(argdict["stkfile"], "r")))[0]
    bedfile = argdict["bedfile"]
    locibed = argdict["locibed"]
    
    numfile = argdict["numfile"]
    locinums = argdict["locinums"]
    
    
    orgtype = "euk"
    if "orgtype" in argdict:
        orgtype = argdict["orgtype"]
        
    locistk = list(readrnastk(open(argdict["locistk"], "r")))[0]
    sizefactor = defaultdict(lambda: 1)
    if argdict["sizefactors"]:
        sizefactor = getsizefactors(argdict["sizefactors"]) 
        for currsample in sampledata.getsamples():
            if currsample not in sizefactor:
                print("Size factor file "+argdict["sizefactors"]+" missing "+currsample, file=sys.stderr)
                sys.exit(1)
    combinereps = argdict["combinereps"]
    allcoveragefile = None
    if "allcoverage" not in argdict or argdict["allcoverage"] == "stdout":
        allcoveragefile = sys.stdout
    else:
        allcoveragefile = open(argdict["allcoverage"],"w")
    samples = sampledata.getsamples()
    minextend = None
    if argdict["minextend"]:
        minextend = int(argdict["minextend"])

    alltrnas = list()
    #gettnanums
    lociedgemargin = argdict["lociedgemargin"]
    if orgtype != "euk" and os.path.isfile(numfile):
        positionnums = list(readtrnanums(numfile, margin = edgemargin))
        locipositionnums = list(readtrnanums(locinums, margin = lociedgemargin))
    else:
        positionnums = gettnanums(trnastk, margin = edgemargin, orgtype = orgtype)
        locipositionnums = gettnanums(locistk, margin = lociedgemargin, orgtype = orgtype)

    #print(orgtype)
    #print >>sys.stderr, locipositionnums
    #print(trnastk.aligns["tRNA-Arg-TCG-1"])
    
    trnastk = trnastk.addmargin(edgemargin)
    locistk = locistk.addmargin(lociedgemargin)

    
    try:
        locitrnas = list()
        basetrnas = list()
        for currfile in bedfile:
            basetrnas.extend(list(currbed for currbed in readbed(currfile)))
        for currfile in locibed:
            locitrnas.extend(list(currbed for currbed in readbed(currfile)))
    except IOError as e:
        print(e, file=sys.stderr)
        sys.exit()
    trnalist = list(curr.addmargin(edgemargin) for curr in basetrnas)
    locilist = list(curr.addmargin(lociedgemargin) for curr in locitrnas)

    featcount = defaultdict(int)

    maxoffset = 10

    #threadmode = False

    
    
    #print >>sys.stderr, ",".join(curr.name for curr in locilist if "Ala" in curr.name)
    #sys.exit(1)
    coveragetable = open(argdict["allcoverage"], "w")
    print("\t".join(["Feature","Sample","position","coverage","readstarts","readends","uniquecoverage","multitrnacoverage","multianticodoncoverage","multiaminocoverage","tRNAreadstotal","actualbase","mismatchedbases","deletedbases","adenines","thymines","cytosines","guanines","deletions"]), file=coveragetable)

    locicoveragetable = open(argdict["locicoverage"], "w")
    print("\t".join(["tRNA_name","sample","position","coverage", "total"]), file=locicoveragetable)
    #uniqcoveragetable = open(argdict["uniqcoverage"], "w")
    mismatchcomparetable = open("mismatchcompare.txt", "w")

    #I'm chunking the tRNAs into 100-long pieces to avoid memory issues with large tRNA sets
    coveragepool = Pool(processes = cores)
    trnachunksize = 50
    for trnanum in range(0,len(trnalist),trnachunksize):
         
        trackargs = list()
        lociargs = list()
        trackuniqargs = list()
        locicoverages = dict()
        samplecoverages = dict()
        uniquecoverages = dict()
        endpoint = trnanum + trnachunksize
        if len(trnalist) < endpoint:
            endpoint = len(trnalist) - 1
        #print >>sys.stderr, "chunk: "+str(trnanum)+"-"+str(endpoint)
        trnasubset = trnalist[trnanum:endpoint]
        locisubset = locilist[trnanum:endpoint]
        if not threadmode:
            for currsample in samples:
                samplecoverages[currsample] = getsamplecoverage(currsample, sampledata, trnasubset,  trnaseqs,  maxmismatches = maxmismatches, minextend = minextend)
                #uniquecoverages[currsample] = getsamplecoverage(currsample, sampledata, trnalist,  trnaseqs,  maxmismatches = maxmismatches, minextend = minextend, uniqueonly = True)
                locicoverages[currsample] =   getlocicoverage(  currsample, sampledata, locisubset,  maxmismatches = maxmismatches, minextend = minextend)
        
        else:
            for currsample in samples:
                trackargs.append(compressargs(currsample, sampledata, trnasubset,  trnaseqs,  maxmismatches = maxmismatches, minextend = minextend))
                #trackuniqargs.append(compressargs(currsample, sampledata, trnasubset,  trnaseqs,  maxmismatches = maxmismatches, minextend = minextend, uniqueonly = True))
        
                lociargs.append(compressargs(  currsample, sampledata, locisubset,  maxmismatches = maxmismatches, minextend = minextend))
        
            results = coveragepool.map(makecoveragepool, trackargs)
            for i, currsample in enumerate(samples):
                samplecoverages[currsample] = results[i]
            lociresults = coveragepool.map(makelocicoveragepool, lociargs)
            for i, currsample in enumerate(samples):
                locicoverages[currsample] = lociresults[i]
            #uniqresults = coveragepool.map(uniquecoverages, trackuniqargs)
            #for i, currsample in enumerate(samples):
            #    uniquecoverages[currsample] = uniqresults[i]
        #print >>sys.stderr, samplecoverages.values()

        
        transcriptcoverage(samplecoverages, coveragetable, trnasubset,sampledata,sizefactor, mincoverage,trnastk, positionnums)
        
        #transcriptcoverage(uniquecoverages, uniqcoveragetable, trnalist,sampledata,sizefactor, mincoverage,trnastk, positionnums)
        
        
        if False:
            getdiffs(samplecoverages, mismatchcomparetable, trnasubset,sampledata, sizefactor,mincoverage,trnastk, positionnums)
        
        
        
        locuscoverage(locicoverages, locicoveragetable, locisubset,sampledata, sizefactor,mincoverage, locistk, locipositionnums)

    coveragetable.close()
    mismatchcomparetable.close()
    locicoveragetable.close()
        
def main(**argdict):
    #print >>sys.stderr, argdict
    argdict = defaultdict(lambda: None, argdict)
    if "edgemargin" not in  argdict:
        edgemargin = 0
    else:
        edgemargin = int(argdict["edgemargin"])
    #currently crashes if set to zero
    if "mincoverage" not in  argdict:
        mincoverage = 10
    else:
        mincoverage = int(argdict["mincoverage"])  
    
        
    sampledata = samplefile(argdict["samplefile"])

    maxmismatches = argdict["maxmismatches"]
    uniquename = argdict["uniquename"]
    uniquegenome = argdict["uniquegenome"]
    trnastk = list(readrnastk(open(argdict["stkfile"], "r")))[0]
    bedfile = argdict["bedfile"]
    sizefactor = defaultdict(lambda: 1)
    if argdict["sizefactors"]:
        sizefactor = getsizefactors(argdict["sizefactors"]) 
        for currsample in sampledata.getsamples():
            if currsample not in sizefactor:
                print("Size factor file "+argdict["sizefactors"]+" missing "+currsample, file=sys.stderr)
                sys.exit(1)
    combinereps = argdict["combinereps"]
    allcoveragefile = None
    if "allcoverage" not in argdict or argdict["allcoverage"] == "stdout":
        allcoveragefile = sys.stdout
    else:
        allcoveragefile = open(argdict["allcoverage"],"w")
    samples = sampledata.getsamples()
    minextend = None
    if argdict["minextend"]:
        minextend = int(argdict["minextend"])

    alltrnas = list()
    #gettnanums
    
    positionnums = gettnanums(trnastk, margin = edgemargin)
    trnastk = trnastk.addmargin(edgemargin)
    
    try:
        basetrnas = list()
        for currfile in bedfile:
            basetrnas.extend(list(currbed for currbed in readbed(currfile)))       
    except IOError as e:
        print(e, file=sys.stderr)
        sys.exit()
    
    trnalist = list(curr.addmargin(edgemargin) for curr in basetrnas)

    featcount = defaultdict(int)

    maxoffset = 10

    allcoverages = dict()
    multaminocoverages = dict()
    multaccoverages = dict()
    multtrnacoverages = dict()
    uniquecoverages = dict()    
    uniquegenomecoverages = dict()  
    multigenomecoverages = dict()
    for currsample in samples:
        currbam = sampledata.getbam(currsample)
        allcoverages[currsample] = dict()
        multaminocoverages[currsample] = dict()
        multaccoverages[currsample] = dict()
        multtrnacoverages[currsample] = dict()
        uniquecoverages[currsample] = dict()
        uniquegenomecoverages[currsample] = dict()
        multigenomecoverages[currsample] = dict()
        try:
            #print >>sys.stderr, currbam
            if not os.path.isfile(currbam+".bai"):
                pysam.index(""+currbam)
            bamfile = pysam.Samfile(""+currbam, "rb" )  
        except IOError as xxx_todo_changeme:
            ( strerror) = xxx_todo_changeme
            print(strerror, file=sys.stderr)
            sys.exit()
            
        for i, currfeat in enumerate(trnalist):
            allcoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
            multaminocoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
            multaccoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
            multtrnacoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
            uniquecoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
            uniquegenomecoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
            multigenomecoverages[currsample][currfeat.name] = readcoverage(trnalist[i])
            if trnalist[i].name == "tRNA-Leu-AAG-1-2":
                pass
                #print >>sys.stderr, "**|"
            for currread in getbamrange(bamfile, trnalist[i],maxmismatches = maxmismatches,allowindels = False):

                if  trnalist[i].coverage(currread) > 10:
                    

                    if minextend is not None and not (currread.start + minextend <= trnalist[i].start or currread.end - minextend >= trnalist[i].end):
                        continue

  
                    allcoverages[currsample][trnalist[i].name].addread(currread)
                    if not isuniqueaminomapping(currread):
                        multaminocoverages[currsample][trnalist[i].name].addread(currread)
                    elif not isuniqueacmapping(currread):
                        multaccoverages[currsample][trnalist[i].name].addread(currread)
                    elif not isuniquetrnamapping(currread):
                        multtrnacoverages[currsample][trnalist[i].name].addread(currread)
                    else:
                        uniquecoverages[currsample][trnalist[i].name].addread(currread)
                    if issinglemapped(currread):
                        uniquegenomecoverages[currsample][trnalist[i].name].addread(currread)
                    else:
                        multigenomecoverages[currsample][trnalist[i].name].addread(currread)
                    
                else:
                    pass
    covfiles = dict()
    print("Feature"+"\t"+"Sample"+"\t"+"\t".join(positionnums), file=allcoveragefile)
        
    for currfeat in trnalist:

        totalreads = sum(allcoverages[currsample][currfeat.name].totalreads for currsample in samples)
        if currfeat.name == "tRNA-Glu-TTC-5-1":
            pass
            #print >>sys.stderr, "**||"+str(totalreads)
        if totalreads < mincoverage:
            continue
        if combinereps:
    
            replicates = sampledata.allreplicates()
            for currrep in replicates:
                print(currfeat.name+"\t"+currrep+"\t"+"\t".join(str(curr) for curr in sumsamples(allcoverages,sampledata,currrep,currfeat,sizefactor)), file=allcoveragefile)
            
        else:
            for currsample in samples:
                print(currfeat.name+"\t"+currsample+"\t"+"\t".join(str(curr) for curr in allcoverages[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name],sizefactor = sizefactor[currsample])), file=allcoveragefile)
    if uniquegenome:
        covfiles = {uniquegenome + '-uniquegenomecoverages.txt':uniquegenomecoverages,uniquegenome + '-multgenomecoverages.txt':multigenomecoverages}

        for filename, currcoverage in covfiles.items():
            covfile = open(filename, "w")
            print("Feature"+"\t"+"Sample"+"\t"+"\t".join(positionnums), file=covfile)
                
            for currfeat in trnalist:
                totalreads = sum(allcoverages[currsample][currfeat.name].totalreads for currsample in samples)
                if totalreads < mincoverage:
                    continue
                if combinereps:
            
                    replicates = sampledata.allreplicates()
                    for currrep in replicates:
                        print(currfeat.name+"\t"+currrep+"\t"+"\t".join(str(curr) for curr in sumsamples(currcoverage,sampledata,currrep,currfeat,sizefactor)), file=covfile)
                    
                else:
                    for currsample in samples:
                        print(currfeat.name+"\t"+currsample+"\t"+"\t".join(str(curr) for curr in currcoverage[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name],sizefactor = sizefactor[currsample])), file=covfile)
            covfile.close()
    if uniquename:
        covfiles = {uniquename + '-uniquecoverages.txt':uniquecoverages,uniquename + '-multaccoverages.txt':multaccoverages,uniquename + '-multtrnacoverages.txt':multtrnacoverages,uniquename + '-multaminocoverages.txt':multaminocoverages}
        
        for filename, currcoverage in covfiles.items():
            covfile = open(filename, "w")
            print("Feature"+"\t"+"Sample"+"\t"+"\t".join(positionnums), file=covfile)
                
            for currfeat in trnalist:

                totalreads = sum(allcoverages[currsample][currfeat.name].totalreads for currsample in samples)

                if totalreads < mincoverage:
                    continue
                if combinereps:
            
                    replicates = sampledata.allreplicates()
                    for currrep in replicates:
                        print(currfeat.name+"\t"+currrep+"\t"+"\t".join(str(curr) for curr in sumsamples(currcoverage,sampledata,currrep,currfeat,sizefactor)), file=covfile)
                    
                else:
                    for currsample in samples:
                        print(currfeat.name+"\t"+currsample+"\t"+"\t".join(str(curr) for curr in currcoverage[currsample][currfeat.name].coveragealign(trnastk.aligns[currfeat.name],sizefactor = sizefactor[currsample])), file=covfile)
            covfile.close()
        
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
    parser.add_argument('--bedfile',  nargs='+', default=list(),
                       help='bed file with mature tRNA features')
    parser.add_argument('--samplefile',
                       help='Sample file in format')
    parser.add_argument('--stkfile',
                       help='Stockholm file')
    parser.add_argument('--sizefactors',
                       help='Optional file including size factors that will be used for normalization')
    parser.add_argument('--combinereps', action="store_true", default=False,
                       help='Sum samples that are replicates')
    parser.add_argument('--edgemargin', type=int, default=0,
                       help='margin to add to feature coordinates')
    
    parser.add_argument('--mincoverage', type=int, default=10,
                       help='Reads with less then this are filtered out (default 10)')
    parser.add_argument('--uniquename',
                       help='Name for files showing unique and non-unique tRNA reads')
    parser.add_argument('--uniquegenome',
                       help='Name for files showing unique and non-unique genome reads')
    parser.add_argument('--maxmismatches', default=None,
                       help='Set maximum number of allowable mismatches')
    '''
    parser.add_argument('--trnapositions', action="store_true", default=False,
                       help='Use tRNA positions')
    '''
    
    
    '''
    Perform check on sizefactor file to ensure it has all samples
    '''
    args = parser.parse_args()
    argdict = vars(args)
    main(**argdict)
