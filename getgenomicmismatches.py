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
gapchars = set("-._~")

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
                
def cigarreadcoverage(cigar):
    nextsum = 1
    for curr in cigar:
        if curr[0] == bam_cins:
            for i in range(curr[1]):
                yield 0
                
        elif curr[0] == bam_cdel:
            nextsum += curr[1]
        elif curr[0] == bam_match:
            for i in range(curr[1]):
                yield nextsum
                nextsum = 1
                
def cigarreadinserts(cigar):
    nextsum = 1
    for curr in cigar:
        if curr[0] == bam_cins:
            for i in range(curr[1]):
                yield False
                
        elif curr[0] == bam_cdel:
            pass
        elif curr[0] == bam_match:
            for i in range(curr[1]):
                yield True
    
def faitobed(fastafile):
    chromlist = list()
    for currline in open(fastafile):
        fields = currline.split("\t")
        currchrom = GenomeRange("dbname", fields[0], 0, fields[1])
        chromlist.append(currchrom)
        pass
    return chromlist
    

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
        

def transcriptcoverage(samplecoverages, mismatchreport, genelist,sampledata,geneseqs,sizefactor, mincoverage):

    print >>mismatchreport, "\t".join(["Feature","Sample","position","coverage","readstarts","readends","readstotal","expreadstotal","actualbase","mismatchedbases","deletedbases","adenines","thymines","cytosines","guanines","deletions"])
    #print >>sys.stderr,mismatchreport
    print >>sys.stderr,"||***"
    samples = sampledata.getsamples()
    for currfeat in genelist:
        #print >>sys.stderr, totalreads
        #print >>sys.stderr, samplecoverages[list(samples)[0]].allcoverages.keys()
        totalreads = sum(samplecoverages[currsample].allcoverages[currfeat.name].totalreads for currsample in samples)
        #print >>sys.stderr, totalreads
        if totalreads < mincoverage:
            continue
        reportpositions = set()  
        #print >>sys.st
        mismatchpos = dict()
        coveragepos = dict()
        trnalen = None
        
        for currsample in samples:
            
            
            covcounts =  list(curr if curr is not None else 0 for curr in samplecoverages[currsample].allcoverages[currfeat.name].coveragelist())
            mismatches =  list(curr if curr is not None else 0 for curr in samplecoverages[currsample].readmismatches[currfeat.name].coveragelist())
            trnalen = len(mismatches)
            mismatchpos[currsample] = list(mismatches[i] for i in range(trnalen))
            coveragepos[currsample] = list(covcounts[i] for i in range(trnalen))
            
        maxmismatch = dict()
        maxcov = dict()
        maxpercent = dict()
        for currpos in range(trnalen):
            maxmismatch[currpos] = max(mismatchpos[currsample][currpos] for currsample in samples)
            maxcov[currpos] = max(coveragepos[currsample][currpos] for currsample in samples)
            maxpercent[currpos] = max((0+1.*mismatchpos[currsample][currpos])/(10+coveragepos[currsample][currpos]) for currsample in samples)
            #print >>sys.stderr, currfeat.name+":"+str(currpos)+":"+str(maxmismatch[currpos])+"/"+str(maxcov[currpos])+":"+str(maxpercent[currpos])

        for currsample in samples:
            #readcounts = samplecoverages[currsample].readcounts
            
            #print >>sys.stderr, ",".join(list(str(curr) if curr is not None else 0 for curr in samplecoverages[currsample].readmismatches[currfeat.name].coveragealign(trnastk.aligns[currfeat.name])))
            #covcounts = list(curr/(1.*samplecoverages[currsample].readcounts[currfeat.name]) if curr is not None else 0 for curr in samplecoverages[currsample].readmismatches[currfeat.name].coveragealign(trnastk.aligns[currfeat.name]))
            
            covcounts =  list(curr/sizefactor[currsample] if curr is not None else 0 for curr in samplecoverages[currsample].allcoverages[currfeat.name].coveragelist())

            mismatches =  list(curr/sizefactor[currsample] if curr is not None else 0 for curr in samplecoverages[currsample].readmismatches[currfeat.name].coveragelist())
            deletions =  list(curr/sizefactor[currsample] if curr is not None else 0 for curr in samplecoverages[currsample].readskips[currfeat.name].coveragelist())


            #uniquecounts =  list(curr/sizefactor[currsample] if curr is not None else 0 for curr in samplecoverages[currsample].uniquecoverages[currfeat.name].coveragelist())
            #multitrna =  list(curr/sizefactor[currsample] if curr is not None else 0 for curr in samplecoverages[currsample].multtrnacoverages[currfeat.name].coveragelist())
            #multaccounts =  list(curr/sizefactor[currsample] if curr is not None else 0 for curr in samplecoverages[currsample].multaccoverages[currfeat.name].coveragelist())
            #multaminocounts =  list(curr/sizefactor[currsample] if curr is not None else 0 for curr in samplecoverages[currsample].multaminocoverages[currfeat.name].coveragelist())



            allstarts  = list(curr/sizefactor[currsample] if curr is not None else 0 for curr in samplecoverages[currsample].readstarts[currfeat.name].coveragelist())
            allends = list(curr/sizefactor[currsample] if curr is not None else 0 for curr in samplecoverages[currsample].readends[currfeat.name].coveragelist())
            allcovcount  = list(curr/sizefactor[currsample] if curr is not None else 0 for curr in samplecoverages[currsample].allcoverages[currfeat.name].coveragelist())
            adeninecount  = list(curr if curr is not None else 0 for curr in samplecoverages[currsample].adeninemismatches[currfeat.name].coveragelist())
            thyminecount = list(curr if curr is not None else 0 for curr in samplecoverages[currsample].thyminemismatches[currfeat.name].coveragelist())
            cytosinecount = list(curr if curr is not None else 0 for curr in samplecoverages[currsample].cytosinemismatches[currfeat.name].coveragelist())
            guanosinecount  = list(curr if curr is not None else 0 for curr in samplecoverages[currsample].guanosinemismatches[currfeat.name].coveragelist())
            readskipcount  = list(curr if curr is not None else 0 for curr in samplecoverages[currsample].readskips[currfeat.name].coveragelist())

            for i, currcount in enumerate(allcovcount):
                mismatchthreshold = .05
                realbase = geneseqs[currfeat.name][i]
                #print >>sys.stderr, len(geneseqs[currfeat.name])
                #print >>sys.stderr, 
                if maxpercent[i] < mismatchthreshold:
                    continue
                if realbase in gapchars:
                    
                    realbase = "-"
                if realbase == "U":
                    realbase = "T"
                print >>mismatchreport, "\t".join([currfeat.name,currsample,str(i),str(covcounts[i]),str(allstarts[i]),str(allends[i]),str(1.*samplecoverages[currsample].readcounts[currfeat.name]/sizefactor[currsample]),str(totalreads),realbase,str(mismatches[i]),str(deletions[i]),str(adeninecount[i]),str(thyminecount[i]),str(cytosinecount[i]),str(guanosinecount[i]), str(readskipcount[i])])
    #sys.exit(1)  
    
def getsamplecoverage(currsample, sampledata, genelist,geneseqs,maxmismatches = None, minextend = None): 
    
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
    
    skipped = 0
    total = 0
    
    try:
        #print >>sys.stderr, currbam
        if not os.path.isfile(currbam+".bai"):
            pysam.index(""+currbam)
        bamfile = pysam.Samfile(""+currbam, "rb" )  
    except IOError as ( strerror):
        print >>sys.stderr, strerror
        sys.exit()
        
    for i, currfeat in enumerate(genelist):

        allcoverages[currfeat.name] = readcoverage(genelist[i])
        uniquegenomecoverages[currfeat.name] = readcoverage(genelist[i])
        multigenomecoverages[currfeat.name] = readcoverage(genelist[i])
        
        readstarts[currfeat.name] = readcoverage(genelist[i])
        readends[currfeat.name] = readcoverage(genelist[i])
        readmismatches[currfeat.name] = readcoverage(genelist[i])
        adeninemismatches[currfeat.name] =   readcoverage(genelist[i])
        thyminemismatches[currfeat.name] =   readcoverage(genelist[i])
        cytosinemismatches[currfeat.name] =  readcoverage(genelist[i])
        guanosinemismatches[currfeat.name] = readcoverage(genelist[i])
        readskips[currfeat.name] = readcoverage(genelist[i])
        
        trimreadcoverage[currfeat.name] =  readcoverage(genelist[i])
        trimreadmismatches[currfeat.name] =  readcoverage(genelist[i])
        readcounts[currfeat.name] = 0

        #print >>sys.stderr, trnalist[i]
        for currread in getbam(bamfile, genelist[i]):
            

            if maxmismatches is not None and currread.getmismatches() > maxmismatches:
                continue
            #print >>sys.stderr, "||**||"+str(currread.getmismatches())
            if genelist[i].coverage(currread) > 10 and genelist[i].strand == currread.strand:
                
                if minextend is not None and not (currread.start + minextend <= trnalist[i].start or currread.end - minextend >= trnalist[i].end):
                    continue
                total += 1
                trnaname = genelist[i].name
                readstart = currread.getfirst(1)
                readend = currread.getlast(1)
                readcounts[trnaname] += 1
                allcoverages[genelist[i].name].addread(currread)
                readstarts[trnaname].addread(readstart)
                readends[trnaname].addread(readend )

                if currread.issinglemapped():
                    uniquegenomecoverages[genelist[i].name].addread(currread)
                else:
                    multigenomecoverages[genelist[i].name].addread(currread)
                
                currseq = currread.getseq()
                geneseq = geneseqs[genelist[i].name]
                
                
                #if currread.start < currfeat.start  or currread.end < currfeat.end  or genelist[i].strand == "+":  #and currread.end > currfeat.end
                    #pass
                
                if genelist[i].strand == "+":
                    genestart = max([0,currread.start - currfeat.start])
                    geneend = genestart+cigarreflength(currread.getcigar())  - max([0,currfeat.start - currread.start])
                    #geneend = min([genestart+cigarreflength(currread.getcigar()) + 1])  #minimum of either the 
                else:
                    genestart = max([0,currfeat.end - currread.end])
                    geneend = genestart+cigarreflength(currread.getcigar())  - max([0,currread.end - currfeat.end])
                    #geneend = min([genestart+len(currseq) + 1 ,genestart+ cigarreadlength(currread.getcigar()) - 1 ])
                refseq = geneseq[genestart:geneend]
                

                

                #print >>sys.stderr, genestart
                #print >>sys.stderr, geneend
                #print >>sys.stderr, currread.getcigar()
                #
                #print >>sys.stderr, cigarreadlength(currread.getcigar())
                #print >>sys.stderr, len(currseq)
                #print >>sys.stderr, refseq
                
                if genelist[i].strand == "+":
                    readstart = max([0,currfeat.start - currread.start])
                    
                    #readend = readstart + cigarreflength(currread.getcigar())
                    readend = min([readstart +len(refseq)  ,readstart + cigarreflength(currread.getcigar())])
                else:
                    readstart = max([0,currread.end - currfeat.end])
                    #readend = readstart + cigarreflength(currread.getcigar())
                    readend = min([readstart +len(refseq)  ,readstart + cigarreflength(currread.getcigar())])
                
                #print >>sys.stderr, readcov
                readcov = list(cigarrefcoverage(currread.getcigar()))
                if genelist[i].strand == "-":
                    readcov = list(cigarrefcoverage(reversed(currread.getcigar())))
                    #readcov = reversed(readcov)
                    pass
                
                alignseq = "".join(currseq[sum(readcov[0:i])] if readcov[i] > 0 else "-" for i in range(cigarreflength(currread.getcigar())))
                alignseq = alignseq[readstart:readend]
                #print >>sys.stderr, alignseq
                #
                #refseq = ""
                #if genestart < 0:
                #    refset = refseq + ("-"*(-genestart))
                #    genestart = 0
                #
                #if geneend >= len(geneseqs[currfeat.name]):
                #    
                #    geneend = len(geneseqs[currfeat.name]) - 1
                #refseq = refseq + geneseqs[currfeat.name][genestart:geneend]
                
                
                #if currread.name == "NB501427:156:H2F7MAFXY:3:21609:5222:18520":
                #if len(refseq) != 0:
                #    print >>sys.stderr, currread.getcigar()
                #    print >>sys.stderr, currfeat.start
                #    print >>sys.stderr, currread.start 
                #    print >>sys.stderr, genestart
                #    print >>sys.stderr, currfeat.name
                #    print >>sys.stderr, currfeat.length()
                #    print >>sys.stderr, currseq
                #    
                #if cigarreflength(currread.getcigar()) < cigarreadlength(currread.getcigar()):
                #    print >>sys.stderr, currread.name
                #    print >>sys.stderr, currread.getcigar()
                #    print >>sys.stderr, "".join(str(curr) for curr in readcov)
                #    print >>sys.stderr, alignseq 
                #    print >>sys.stderr,refseq
                #    print >>sys.stderr, currseq
                #if len(refseq) != len(currseq):
                #    print >>sys.stderr, currread.name
                #    print >>sys.stderr,refseq
                #    print >>sys.stderr, currseq
                skipends = True
                
                
                #need to check mismatches later
                #if  currread.name == "NB501427:404:HJ3WGAFX2:3:11602:17144:17931":
                
              
                #if refseq != alignseq and len(currread.getcigar()) > 1 and genelist[i].strand == "-":# or cigarreflength(currread.getcigar()) != len(refseq):
                #    #skipped += 1
                #    print >>sys.stderr, currread.name
                #    print >>sys.stderr, currread.getcigar() 
                #    print >>sys.stderr, genelist[i].strand
                #    print >>sys.stderr,currseq
                #    print >>sys.stderr, "gene: "+str(currfeat.start)+"-"+str(currfeat.end)
                #    print >>sys.stderr, "read:" +str(currread.start)+"-"+str(currread.end)
                #    print >>sys.stderr,  "gene: "+str(genestart)+"-"+str(geneend)
                #    print >>sys.stderr,  "read: " +str(readstart)+"-"+str(readend)
                #    print >>sys.stderr,refseq
                #    print >>sys.stderr, alignseq
                #    pass
                #    #continue
                #else:
                #    #continue
                #    pass
                for currpos in range(len(refseq)): #30
                    
                    if currpos < 0 or currpos >= len(refseq):
                        #print >>sys.stderr, currread.name
                        #print >>sys.stderr, currpos
                        #print >>sys.stderr, len(refseq)
                        #print >>sys.stderr, cigarreflength(currread.getcigar())
                        #
                        #print >>sys.stderr,refseq
                        #print >>sys.stderr, currseq
                        pass
                        
                        
                    currbase = alignseq[currpos]


                    refbase = refseq[currpos]
                    if refbase not in gapchars:

                        if refbase != currbase:
                            #if (currpos + currread.start) - readmismatches[trnaname].region.start < 0:
                            #    print >>sys.stderr, "before start: "+str(currpos)+"+"+str(currread.start) +"-"+str(readmismatches[trnaname].region.start)
                            #    #base - self.region.start
                            

                            readmismatches[trnaname].addbase(currread.start + currpos)
                        #allcoverages[trnaname].addbase(currread.start + currpos)
                        if skipends:
                            if currpos < 3 or currpos > cigarreflength(currread.getcigar()) - 3 :
                                continue
                        #allcoverages[genelist[i].name].addbase(currread.start + currpos)
                        #if currpos > 3:
                        #    trimreadcoverage[trnaname].addbase(currread.start + currpos)
                        #    if refbase != currbase:
                        #        trimreadmismatches[trnaname].addbase(currread.start + currpos)
                        if currbase == "-":
                            readskips[trnaname].addbase(currread.start + currpos)
                        if currbase == "A":
                            adeninemismatches[trnaname].addbase(currread.start + currpos)
                        elif currbase == "T":
                            thyminemismatches[trnaname].addbase(currread.start + currpos)
                        elif currbase == "C":
                            cytosinemismatches[trnaname].addbase(currread.start + currpos)
                        elif currbase == "G":
                            guanosinemismatches[trnaname].addbase(currread.start + currpos)
    #print >>sys.stderr, currsample+":" +str(skipped)+"/"+str(total)+":"+str(((1.*skipped)/total))
    return coverageinfo( readcounts, allcoverages,readstarts, readends,multaminocoverages, multaccoverages, multtrnacoverages,uniquecoverages, uniquegenomecoverages,multigenomecoverages, readmismatches,adeninemismatches,thyminemismatches,cytosinemismatches, guanosinemismatches,readskips,trimmismatches = trimreadmismatches, trimcoverage = trimreadcoverage  )


def makecoveragepool(args):
    return getsamplecoverage(*args[0], **args[1])
def makelocicoveragepool(args):
    return getlocicoverage(*args[0], **args[1])
def compressargs( *args, **kwargs):
    return tuple([args, kwargs])
def main(**argdict):
    #print >>sys.stderr, argdict
    argdict = defaultdict(lambda: None, argdict)
    if "edgemargin" not in  argdict:                    
        edgemargin = 0
    else:
        edgemargin = int(argdict["edgemargin"])
    #currently crashes if set to zero
    if "mincoverage" not in  argdict:
        mincoverage = 100
    else:
        mincoverage = int(argdict["mincoverage"])  
    
    if "bamdir" not in argdict:
        bamdir = "./"
    bamdir = argdict["bamdir"]
    sampledata = samplefile(argdict["samplefile"], bamdir = bamdir)
    samples = sampledata.getsamples()
    genomefasta = argdict["genomefasta"]
    chroms = faitobed(genomefasta+".fai")
    ensemblgtf = argdict["ensemblgtf"]
    bedfiles = argdict["bedfile"]
    if bedfiles is None:
        bedfiles = list()
    coveragefile = argdict["covfile"]
    if coveragefile is None:
        print >>sys.stderr, "no output file (--covfile)"
        sys.exit(1)
    embllist = list()
    ensemblgtf = os.path.expanduser(ensemblgtf)
    bedlist = list()
    if ensemblgtf is not None:    
        embllist = list(readgtf(ensemblgtf, filtertypes = set(), seqfile= genomefasta))
    for currbed in bedfiles:
        bedlist.extend(readbed(currbed, seqfile= genomefasta))
    genelist = embllist + bedlist
    
    chromnames = set(curr.chrom for curr in chroms)
    genelist = list(curr for curr in genelist if curr.chrom in chromnames)
    geneseqs = getseqdict(genelist, faifiles = {"genome":genomefasta+".fai"})
    maxmismatches = argdict["maxmismatches"]
    cores = 8
    if argdict["cores"] is not None:
        cores = int(argdict["cores"])
    
    threadmode = True
    if cores == 1:
        threadmode = False
        
    #uniquegenome = argdict["uniquegenome"]

    #print >>sys.stderr, geneseqs
    
    #sys.exit()
    sizefactor = defaultdict(lambda: 1)
    if argdict["sizefactors"]:
        sizefactor = getsizefactors(argdict["sizefactors"]) 
        for currsample in sampledata.getsamples():
            if currsample not in sizefactor:
                print >>sys.stderr, "Size factor file "+argdict["sizefactors"]+" missing "+currsample
                sys.exit(1)
    combinereps = True


    #gettnanums

    #print(orgtype)
    #print(positionnums)
    #print(trnastk.aligns["tRNA-Arg-TCG-1"])

    featcount = defaultdict(int)

    maxoffset = 10
    samplecoverages = dict()
    #threadmode = False
    trackargs = list()
    lociargs = list()
    
    locicoverages = dict()
    print >>sys.stderr, threadmode
    #sys.exit(1)
    if not threadmode:
        for currsample in samples:
            samplecoverages[currsample] = getsamplecoverage(currsample, sampledata, genelist,geneseqs,maxmismatches = maxmismatches)
            
    else:
        coveragepool = Pool(processes = cores)
        for currsample in samples:
            trackargs.append(compressargs(currsample, sampledata, genelist,geneseqs,maxmismatches = maxmismatches))
 
        results = coveragepool.map(makecoveragepool, trackargs)
        for i, currsample in enumerate(samples):
            samplecoverages[currsample] = results[i]
    coveragetable = open(coveragefile, "w")
    transcriptcoverage(samplecoverages, coveragetable, genelist,sampledata,geneseqs,sizefactor, mincoverage)
 
            
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')

    parser.add_argument('--samplefile',
                       help='Sample file in format')
    parser.add_argument('--genomefasta',
                       help='Genome file')
    parser.add_argument('--ensemblgtf',
                       help='ensemblgtf')
    parser.add_argument('--covfile',
                       help='output coverage file')
    parser.add_argument('--cores',
                       help='cores')
    parser.add_argument('--sizefactors',
                       help='Optional file including size factors that will be used for normalization')
    parser.add_argument('--bedfile',  nargs='*',
                   help='Additional bed files for feature list')

    parser.add_argument('--mincoverage', type=int, default=10,
                       help='Reads with less then this are filtered out (default 10)')
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