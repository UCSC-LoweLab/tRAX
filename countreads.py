#!/usr/bin/env python

import pysam
import sys
import argparse
import os.path
from collections import defaultdict
from trnasequtils import *
import itertools
import threading
import time
from multiprocessing import Process, Queue, Pool



def getdupes(namelist):
    allset = set()
    for currname in namelist:
        if currname in allset:
            yield currname
        else:
            allset.add(currname)

def enddict():
    return defaultdict(int)
class featurecount:
    def __init__(self, samplename, bamfile, trnas = list(), trnaloci = list(), emblgenes = list(), otherfeats = list()):
        self.samplename = samplename
        self.bamfile = bamfile
        self.trnas = trnas
        self.trnaloci = trnaloci
        self.emblgenes = emblgenes
        self.otherfeats = otherfeats
        
        self.counts = defaultdict(int)
        self.trnacounts = defaultdict(int)
        self.antitrnacount = defaultdict(int)
        self.trnawholecounts = defaultdict(int)
        self.trnafivecounts = defaultdict(int)
        self.trnathreecounts = defaultdict(int)
        self.trnalocuscounts = defaultdict(int)
        self.trnalocustrailercounts = defaultdict(int)
        self.partialtrnalocuscounts = defaultdict(int)
        self.fulltrnalocuscounts  = defaultdict(int)
        self.trnauniquecounts = defaultdict(int)
        self.aminocounts  = defaultdict(int)
        self.anticodoncounts =  defaultdict(int) 
        self.trnaendtypecounts = defaultdict(enddict)
        self.lengthsum = defaultdict(int) 
        self.lengthtotal = defaultdict(int) 
        
        self.gcpercent = defaultdict(int) 
        self.gctotal = defaultdict(int) 
        
        self.genetypes = dict()
        
    def setgenetype(self, genename, genetype):
        self.genetypes[genename] = genetype
    def addcount(self, genename):
       self.counts[genename] += 1
    def addantitrnacount(self, genename):
       self.antitrnacount[genename] += 1
       
    def addlocuscount(self, genename):
       self.trnalocuscounts[genename] += 1
    def addpartiallocuscount(self, genename):
       self.partialtrnalocuscounts[genename] += 1
    def addfulllocuscount(self, genename):
       self.fulltrnalocuscounts[genename] += 1
    def addlocustrailercount(self, genename):
       self.trnalocustrailercounts[genename] += 1
    def addtrnacount(self, genename):
       self.trnacounts[genename] += 1
    def adduniquecount(self, genename):
       self.trnauniquecounts[genename] += 1
    def addaminocount(self, amino):
       self.aminocounts[amino] += 1
    def addanticodoncount(self, anticodon):
       self.anticodoncounts[anticodon] += 1
    def addfragcount(self, featname, fragtype):    
        if fragtype == "Whole":
            self.trnawholecounts[featname] += 1
        elif fragtype == "Fiveprime":
            self.trnafivecounts[featname] += 1
        elif fragtype == "Threeprime":
            self.trnathreecounts[featname] += 1
    def addendcount(self, featname, endtype):    
        if endtype is not None:
            self.trnaendtypecounts[featname][endtype] += 1
            
    def addreadlength(self, genename, length):
       self.lengthsum[genename] += length
       self.lengthtotal[genename] += 1
       
    def addgc(self, genename, gc, length):
       self.gcpercent[genename] += gc
       self.gctotal[genename] += length
       
       
    def getgenecount(self, genename):
       return self.counts[genename]
    def getantitrnacount(self, genename):
       return self.antitrnacount[genename]
    def getlocuscount(self, genename):
       return self.trnalocuscounts[genename]
    def getpartiallocuscount(self, genename):
       return self.partialtrnalocuscounts[genename]
    def getfulllocuscount(self, genename):
       return self.fulltrnalocuscounts[genename] 
    def getlocustrailercount(self, genename):
       return self.trnalocustrailercounts[genename]
    def gettrnacount(self, genename):
       return self.trnacounts[genename]
    def getuniquecount(self, genename):
       return self.trnauniquecounts[genename]
    def getaminocount(self, amino):
       return self.aminocounts[amino]
    def getanticodoncount(self, anticodon):
       return self.anticodoncounts[anticodon]

    def getfivecount(self, genename):
       return self.trnafivecounts[genename]
    def getthreecount(self, genename):
       return self.trnathreecounts[genename]
    def getwholecount(self, genename):
       return self.trnawholecounts[genename]
    def getendtypecount(self, genename):
       return self.trnaendtypecounts[genename]
    def getreadavglength(self, genename):
        if self.lengthtotal[genename] == 0:
            return 0
        else:
            return self.lengthsum[genename] / self.lengthtotal[genename]
    def getreadavggc(self, genename):
        if self.gctotal[genename] == 0:
            return 0
        else:
            return self.gcpercent[genename] / self.gctotal[genename]

def getbamcounts(bamfile, samplename,trnainfo, trnaloci, trnalist,featurelist = dict(),otherseqdict = dict(), embllist = list(), bedfiles = list(),nomultimap = False, allowindels = True, maxmismatches = None):
    samplecounts = featurecount(samplename, bamfile, trnas = trnalist, trnaloci = trnaloci, emblgenes = embllist, otherfeats = featurelist)
    fullpretrnathreshold = 2
    minpretrnaextend = 5
    #minimum mapq
    #nomultimap = False
    minmapq = 0
    if nomultimap:
        minmapq = 2
    #minimum number of reads for a feature to be reported
    minreads = 5
    #print >>sys.stderr, embllist
    
    genetypes = dict()
    currbam = bamfile
    
    for currfile in bedfiles:
        bedfeatures = list(readfeatures(currfile, removepseudo = False))
        for curr in bedfeatures:
            genetypes[curr.name] = os.path.basename(currfile)
            
        featurelist[currfile] = bedfeatures
    
    #print >>sys.stderr, currsample
    #doing this thing here why I only index the bamfile if the if the index file isn't there or is older than the map file
    try:
        if not os.path.isfile(currbam+".bai") or os.path.getmtime(currbam+".bai") < os.path.getmtime(currbam):
            pysam.index(""+currbam)
        bamfile = pysam.Samfile(""+currbam, "rb" )  
    except IOError as ( strerror):
        print >>sys.stderr, strerror
        sys.exit(1)
    except pysam.utils.SamtoolsError:
        print >>sys.stderr, "Can not index "+currbam
        print >>sys.stderr, "Exiting..."
        sys.exit(1)
        
    
    for currfile in featurelist.iterkeys():
        for currfeat in featurelist[currfile]:
            #try catch is to account for weird chromosomes and the like that aren't in the genome
            #means that if I can't find a feature, I record no counts for it rather than bailing
            try:
                for currread in getbamrange(bamfile, currfeat, singleonly = nomultimap, maxmismatches = maxmismatches,allowindels = allowindels):
                    if currfeat.coverage(currread) > 10:
                        samplecounts.addcount(currfeat.name)
                        samplecounts.addreadlength(currfeat.name, currread.length())
                        #samplecounts.addgc(currfeat.name, currread.getgc(), currread.length())
                        samplecounts.setgenetype(currfeat.name,os.path.basename(currfile))
            except ValueError:
                pass

    #extra sequences built during database creation (experimental)
    for currtype in otherseqdict.iterkeys():
        for currfeat in otherseqdict[currtype]:
            for currread in getbamrange(bamfile, currfeat, singleonly = nomultimap, maxmismatches = maxmismatches,allowindels = allowindels):
                #print >>sys.stderr, currfeat.name
                samplecounts.addcount(currfeat.name)
                samplecounts.addreadlength(currfeat.name, currread.length())
                #samplecounts.addgc(currfeat.name, currread.getgc(), currread.length())
                samplecounts.setgenetype(currfeat.name,currtype)
    for genename, featset in itertools.groupby(embllist,lambda x: x.data["genename"]):
        #print >>sys.stderr, "**"
        #pass 
        try:
            allreads = set()
            for currfeat in list(featset):
                
                for currread in getbamrangeshort(bamfile, currfeat, singleonly = nomultimap, maxmismatches = maxmismatches,allowindels = allowindels, skiptags = True):
                    #print >>sys.stderr, "**"+currread.name 
                    #continue
                    
                    if currfeat.coverage(currread) > 10:
                        
                        samplecounts.addcount(genename)
                        samplecounts.addreadlength(currfeat.name, currread.length())
                        #samplecounts.addgc(currfeat.name, currread.getgc(), currread.length())
                        #print >>sys.stderr, "**"+currread.name
                        samplecounts.setgenetype(genename,currfeat.data["biotype"])
                        #print >>sys.stderr, currfeat.bedstring()
        except ValueError:
            pass
    for currfeat in trnaloci:
        #print >>sys.stderr,  currfeat.bedstring()
        #print >>sys.stderr,  currfeat.getdownstream(30).bedstring()
        for currread in getbamrangeshort(bamfile, currfeat.addmargin(30), singleonly = nomultimap, maxmismatches = maxmismatches,allowindels = allowindels, skiptags = True):
            #gotta be more than 5 bases off one end to be a true pre-tRNA
            #might want to shove these to the real tRNA at some point, but they are for now just ignored
            
            if currfeat.coverage(currread) > 10 and (currread.start + minpretrnaextend <= currfeat.start or currread.end - minpretrnaextend >= currfeat.end):
                samplecounts.addlocuscount(currfeat.name)
                samplecounts.addreadlength(currfeat.name, currread.length())
                #samplecounts.addgc(currfeat.name, currread.getgc(), currread.length())
                if currread.start + fullpretrnathreshold <  currfeat.start and currread.end - fullpretrnathreshold + 3 >  currfeat.end:
                    samplecounts.addfulllocuscount(currfeat.name)
                else:
                    #partialtrnalocuscounts[currsample][currfeat.name] += 1
                    samplecounts.addpartiallocuscount(currfeat.name)
            elif currfeat.getdownstream(30).coverage(currread) > 10:  #need the elif otherwise fragments that include trailer get in there
                samplecounts.addlocuscount(currfeat.name)
                samplecounts.addlocustrailercount(currfeat.name)
            else:
                #print >>sys.stderr,  currfeat.getdownstream(30).coverage(currread)
                pass
    #print >>sys.stderr, samplename+" threadA "+str(time.time())        
    for currfeat in trnalist:
        #print >>sys.stderr, samplename+":"+currfeat.name
        
        featreads = 0
        for currread in getbam(bamfile, currfeat, singleonly = nomultimap, allowindels = allowindels):
            #samplecounts.addgc(currfeat.name, currread.getgc(), currread.length())
            if maxmismatches is not None and currread.getmismatches() > maxmismatches:
                continue
            samplecounts.addreadlength(currfeat.name, currread.length())

            featreads += 1
            if not currfeat.strand == currread.strand:
                samplecounts.addantitrnacount(currfeat.name)
                continue
            if not currfeat.coverage(currread) > 10:
                continue

            curramino = trnainfo.getamino(currfeat.name)
            curranticodon = trnainfo.getanticodon(currfeat.name)
            #samplecounts.addfragcount(currfeat.name, fragtype)
            samplecounts.addtrnacount(currfeat.name)
                
            fragtype = getfragtype(currfeat, currread)
            samplecounts.addfragcount(currfeat.name, fragtype)
            endtype = getendtype(currfeat, currread)
            #print >>sys.stderr, endtype
            samplecounts.addendcount(currfeat.name, endtype)
            if currread.isuniquetrnamapping():
                samplecounts.adduniquecount(currfeat.name)
            if currread.isuniqueaminomapping():
                pass
            if not currread.isuniqueaminomapping():
                pass
            elif currread.isuniqueacmapping():
                samplecounts.addaminocount(curramino)
            else:
                samplecounts.addaminocount(curramino)
                samplecounts.addanticodoncount(curramino)
        #print >>sys.stderr, str(featreads)+"/"+str(samplecounts.gettrnacount(currfeat.name))
    #print >>sys.stderr, samplename+" thread "+str(time.time())
            
    return samplecounts

def printcountfile(countfile, samples,  samplecounts, trnalist, trnaloci, featurelist, embllist, otherseqdict = dict(),minreads = 5, includebase = False):
    print >>countfile, "\t".join(samples)
    trnanames = set()
    for currfeat in trnalist:
        #print >>sys.stderr, samplecounts
        if max(itertools.chain((samplecounts[currsample].gettrnacount(currfeat.name) for currsample in samples), [0])) < minreads:
            continue
        if includebase:
            print >>countfile, currfeat.name+"\t"+"\t".join(str(samplecounts[currsample].gettrnacount(currfeat.name)) for currsample in samples)
            print >>countfile, currfeat.name+"_antisense\t"+"\t".join(str(samplecounts[currsample].getantitrnacount(currfeat.name)) for currsample in samples)

        else:
            print >>countfile, currfeat.name+"_wholecounts\t"+"\t".join(str(samplecounts[currsample].getwholecount(currfeat.name)) for currsample in samples)
            print >>countfile, currfeat.name+"_fiveprime\t"+"\t".join(str(samplecounts[currsample].getfivecount(currfeat.name) ) for currsample in samples)
            print >>countfile, currfeat.name+"_threeprime\t"+"\t".join(str(samplecounts[currsample].getthreecount(currfeat.name)) for currsample in samples)
            print >>countfile, currfeat.name+"_other\t"+"\t".join(str(samplecounts[currsample].gettrnacount(currfeat.name) - (samplecounts[currsample].getwholecount(currfeat.name) + samplecounts[currsample].getfivecount(currfeat.name) + samplecounts[currsample].getthreecount(currfeat.name))) for currsample in samples)
            

            print >>countfile, currfeat.name+"_antisense\t"+"\t".join(str(samplecounts[currsample].getantitrnacount(currfeat.name)) for currsample in samples)

    for currfeat in trnaloci:
        if max(itertools.chain((samplecounts[currsample].getlocuscount(currfeat.name) for currsample in samples),[0])) < minreads:
            continue
        if includebase:
            print >>countfile, currfeat.name+"\t"+"\t".join(str(samplecounts[currsample].getlocuscount(currfeat.name)) for currsample in samples)
        else:
            print >>countfile, currfeat.name+"_wholeprecounts\t"+"\t".join(str(samplecounts[currsample].getfulllocuscount(currfeat.name) ) for currsample in samples)
            print >>countfile, currfeat.name+"_partialprecounts\t"+"\t".join(str(samplecounts[currsample].getpartiallocuscount(currfeat.name) ) for currsample in samples)
            print >>countfile, currfeat.name+"_trailercounts\t"+"\t".join(str(samplecounts[currsample].getlocustrailercount(currfeat.name)) for currsample in samples)        
    for currbed in featurelist.iterkeys():
        for currfeat in featurelist[currbed]:
            if currfeat.name in trnanames:
                continue
            trnanames.add(currfeat.name)
            if max(samplecounts[currsample].getgenecount(currfeat.name) for currsample in samples) > minreads:
                print >>countfile, currfeat.name+"\t"+"\t".join(str(samplecounts[currsample].getgenecount(currfeat.name)) for currsample in samples)
    for currtype in otherseqdict.iterkeys():        
        for currfeat in otherseqdict[currtype] :
        
            trnanames.add(currfeat.name)
            if max(samplecounts[currsample].getgenecount(currfeat.name) for currsample in samples) > minreads:
                print >>countfile, currfeat.name+"\t"+"\t".join(str(samplecounts[currsample].getgenecount(currfeat.name)) for currsample in samples)
    for currfeat in embllist:
        
        genename = currfeat.data['genename']
        if genename in trnanames:
            continue
        trnanames.add(genename)
        
        if genename is None:
            print >>sys.stderr, currfeat.name
            sys.exit(1)
        #print >>sys.stderr, list(samplecounts[currsample].getgenecount(currfeat.name) for currsample in samples)
        if max(samplecounts[currsample].getgenecount(genename) for currsample in samples) > minreads:
            print >>countfile, genename+"\t"+"\t".join(str(samplecounts[currsample].getgenecount(genename)) for currsample in samples)


def averagesamples(allcounts, genename,samples):
    
    return str(sum(allcounts[currsample].lengthsum[genename] for currsample in samples)/(.01+1.*sum(allcounts[currsample].lengthtotal[genename] for currsample in samples)))
    
def gcsamples(allcounts, genename,samples):
        return str(sum(allcounts[currsample].gcpercent[genename] for currsample in samples)/(.01+1.*sum(allcounts[currsample].gctotal[genename] for currsample in samples)))

def printtypefile(genetypeout,samples, allcounts,trnalist, trnaloci, featurelist, embllist , otherseqdict = dict(),minreads = 5):
    trnanames = set()
    genetypes = dict()
    genelengths = dict()
    for currsample in samples:
        genetypes.update(allcounts[currsample].genetypes)
    for currbed in featurelist.iterkeys():
        for currfeat in featurelist[currbed] :
            if currfeat.name in trnanames:
                continue
            trnanames.add(currfeat.name)
            if max(allcounts[currsample].counts[currfeat.name] for currsample in samples) > minreads:
                print >>genetypeout, currfeat.name+"\t"+genetypes[currfeat.name]   +"\t"+currfeat.chrom+"\t"+averagesamples(allcounts, currfeat.name, samples)
    
        
    for currfeat in trnaloci:
        print >>genetypeout, currfeat.name+"_wholeprecounts"+"\t"+"trna_wholeprecounts" +"\t"+currfeat.chrom +"\t"+averagesamples(allcounts, currfeat.name, samples)
        print >>genetypeout, currfeat.name+"_partialprecounts"+"\t"+"trna_partialprecounts"+"\t"+currfeat.chrom+"\t"+averagesamples(allcounts, currfeat.name, samples)
        print >>genetypeout, currfeat.name+"_trailercounts"+"\t"+"trna_trailercounts"+"\t"+currfeat.chrom+"\t"+averagesamples(allcounts, currfeat.name, samples)
        print >>genetypeout, currfeat.name+""+"\t"+"tRNA_locus"+"\t"+currfeat.chrom+"\t"+averagesamples(allcounts, currfeat.name, samples)
    for currfeat in trnalist:
        print >>genetypeout, currfeat.name+"_wholecounts"+"\t"+"trna_wholecounts"+"\t"+"tRNA"+"\t"+averagesamples(allcounts, currfeat.name, samples)
        print >>genetypeout, currfeat.name+"_fiveprime"+"\t"+"trna_fiveprime"+"\t"+"tRNA"+"\t"+averagesamples(allcounts, currfeat.name, samples)
        print >>genetypeout, currfeat.name+"_threeprime"+"\t"+"trna_threeprime"+"\t"+"tRNA"+"\t"+averagesamples(allcounts, currfeat.name, samples)
        print >>genetypeout, currfeat.name+"_other"+"\t"+"trna_other"+"\t"+"tRNA"+"\t"+averagesamples(allcounts, currfeat.name, samples)
        print >>genetypeout, currfeat.name+"_antisense"+"\t"+"trna_antisense"+"\t"+"tRNA"+"\t"+averagesamples(allcounts, currfeat.name, samples)
        print >>genetypeout, currfeat.name+""+"\t"+"tRNA"+"\t"+"tRNA"+"\t"+averagesamples(allcounts, currfeat.name, samples)
    
    for currfeat in embllist:
        genename = currfeat.data['genename']
        if genename in trnanames:
            continue
        trnanames.add(genename)
        if genename is None:
            #print >>sys.stderr, currfeat.name
            continue
        if max(allcounts[currsample].counts[genename] for currsample in samples) > minreads:
            print >>genetypeout, genename+"\t"+genetypes[genename]        +"\t"+currfeat.chrom+"\t"+averagesamples(allcounts, currfeat.name, samples)
    for currtype in otherseqdict.iterkeys():
        for currfeat in otherseqdict[currtype]:
            genename = currfeat.name
            if genename in trnanames:
                continue
            trnanames.add(genename)
            if genename is None:
                #print >>sys.stderr, currfeat.name
                continue
            if max(allcounts[currsample].counts[genename] for currsample in samples) > minreads:
                print >>genetypeout, genename+"\t"+genetypes[genename]    +"\t"+currfeat.chrom +"\t"+averagesamples(allcounts, currfeat.name, samples)         
              
     

def printtrnauniquecountcountfile(trnauniquefile,samples,  samplecounts, trnalist, trnaloci , minreads = 5):
    trnauniquefile = open(trnauniquefile, "w")
    print >>trnauniquefile, "\t".join(currsample for currsample in samples)
    for currfeat in trnalist:
        if max(samplecounts[currsample].getuniquecount(currfeat.name) for currsample in samples) < minreads:
            continue
        print  >>trnauniquefile, currfeat.name+"\t"+"\t".join(str(samplecounts[currsample].getuniquecount(currfeat.name)) for currsample in samples)
    trnauniquefile.close() 
def printtrnacountfile(trnacountfilename,samples,  samplecounts, trnalist, trnaloci , minreads = 5):
    trnacountfile = open(trnacountfilename, "w")
    print >>trnacountfile, "\t".join(currsample for currsample in samples)
    for currfeat in trnaloci:
        if max(samplecounts[currsample].getlocuscount(currfeat.name) for currsample in samples) < minreads:
            continue
        print >>trnacountfile, currfeat.name+"\t"+"\t".join(str(samplecounts[currsample].getlocuscount(currfeat.name)) for currsample in samples)
    for currfeat in trnalist:
        if max(samplecounts[currsample].gettrnacount(currfeat.name) for currsample in samples) < minreads:
            continue
        print  >>trnacountfile, currfeat.name+"\t"+"\t".join(str(samplecounts[currsample].gettrnacount(currfeat.name)) for currsample in samples)
    trnacountfile.close()
trnaends = list(["CCA","CC","C",""])    
def printtrnaendfile(trnaendfilename,samples,  samplecounts, trnalist, trnaloci , minreads = 5):
    trnaendfile = open(trnaendfilename, "w")
    print >>trnaendfile, "end\t"+"\t".join(currsample for currsample in samples)

    for currfeat in trnalist:
        if max(samplecounts[currsample].gettrnacount(currfeat.name) for currsample in samples) < minreads:
            continue
        for currend in trnaends:
            endstring = currend
            if currend == "":
                endstring = "Trimmed"
            print  >>trnaendfile, currfeat.name+"\t"+endstring+"\t"+"\t".join(str(samplecounts[currsample].getendtypecount(currfeat.name)[currend]) for currsample in samples)
    trnaendfile.close()   
    
def getbamcountsthr(results,currsample, *args, **kwargs):
    results[currsample] = getbamcounts(*args, **kwargs)
def getbamcountsqueue(countqueue,currsample, *args, **kwargs):
    countqueue.put([currsample,getbamcounts(*args, **kwargs)])
    
    
def countreadspool(args):
    return getbamcounts(*args[0], **args[1])
    
def compressargs( *args, **kwargs):
    return tuple([args, kwargs])

def testmain(**argdict):
    trnauniquefilename = None
    argdict = defaultdict(lambda: None, argdict)
    includebase = argdict["nofrag"]
    fullpretrnasonly = argdict["onlyfullpretrnas"]
    trnatable = argdict["trnatable"]
    removepseudo = argdict["removepseudo"]
    ensemblgtf = argdict["ensemblgtf"]
    nomultimap = argdict["nomultimap"]
    if argdict["maxmismatches"] is not None:
        maxmismatches = int(argdict["maxmismatches"])
    else:
        maxmismatches = None
    cores = argdict["cores"]
    trnaendfilename = argdict["trnaends"]
    threadmode = True
    if cores == 1:
        threadmode = False
    otherseqs = extraseqfile(argdict["otherseqs"])
    typefile = None
    
    if "bamdir" not in argdict:
        bamdir = "./"
    bamdir = argdict["bamdir"]
    sampledata = samplefile(argdict["samplefile"], bamdir = bamdir)
    
    bedfiles = list() 
    if "trnauniquecounts" in argdict:
        trnauniquefilename = argdict["trnauniquecounts"]
    if "bedfile"  in argdict:
        bedfiles = argdict["bedfile"]
    trnalocifiles = list()
    if "trnaloci"  in argdict:
        trnalocifiles = argdict["trnaloci"]
    maturetrnas = list()
    if "maturetrnas" in argdict:
        maturetrnas = argdict["maturetrnas"]
        
    #trnalocifiles = argdict["trnaloci"]
    #maturetrnas=argdict["maturetrnas"]
    genetypefile = argdict["genetypefile"]
    trnacountfilename = argdict["trnacounts"]

    
    trnainfo = transcriptfile(trnatable)

    #print >>sys.stderr, bedfiles
    
    alltrnas = list()
    
    
    samplefiles = dict()
    
    
    samples = sampledata.getsamples()
    genetypes = dict()
    fullpretrnathreshold = 2
    otherseqdict = dict()
    #Grabbing all the features to count
    #print >>sys.stderr, otherseqs
    try:
        featurelist = dict()
        trnaloci = list()
        for currfile in bedfiles:
            bedfeatures = list(readfeatures(currfile, removepseudo = removepseudo))
            for curr in bedfeatures:
                genetypes[curr.name] = os.path.basename(currfile)
                
            featurelist[currfile] = bedfeatures
        trnalist = list()
        for currfile in trnalocifiles:
            trnaloci.extend(list(readbed(currfile)))
        for currfile in maturetrnas:
            trnalist.extend(list(readbed(currfile)))
        if ensemblgtf is not None:    
            embllist = list(readgtf(ensemblgtf, filterpsuedo = removepseudo))
        else:
            embllist = list()
        for name, currfile in otherseqs.getseqbeds().iteritems():
            otherseqdict[name] = list(readbed(currfile))

    except IOError as e:
        print >>sys.stderr, e
        sys.exit()
    featcount = defaultdict(int)
    allfeats = trnaloci+trnalist
    if len(set(curr.name for curr in allfeats)) < len(list(curr.name for curr in allfeats )):
        print >>sys.stderr, "Duplicate names in feature list"
    
    
    #featurelist = list(curr for curr in featurelist if curr.name == 'unknown20') 
    #alltrnas = list(curr.name for curr in featurelist)
    #print >>sys.stderr, "***"
    #setting up all the feature count dictionaries
                            
    allcounts = dict()
    threads = dict()
    #threadmode = False
    starttime = time.time()
    #print  list(curr.name for curr in trnalist)
    print >>sys.stderr, "**||"
    print >>sys.stderr, maxmismatches
    
    #sys.exit()
    if threadmode:

        countpool = Pool(processes=cores)
        arglist = list()
        for currsample in samples:
            currbam = sampledata.getbam(currsample)
            arglist.append(compressargs(currbam, currsample,trnainfo, trnaloci, trnalist, otherseqdict = otherseqdict, embllist = embllist, featurelist = featurelist, bedfiles = bedfiles, maxmismatches = maxmismatches))
        #arglist = list((tuple([currsample, sampledata.getbam(currsample)]) for currsample in samples))
        results = countpool.map(countreadspool, arglist)
        for i, curr in enumerate(samples):
            allcounts[curr] = results[i]

    else:

        for currsample in samples:
            
            currbam = sampledata.getbam(currsample)
            allcounts[currsample] = getbamcounts(currbam, currsample,trnainfo, trnaloci, trnalist, otherseqdict = otherseqdict,embllist = embllist, featurelist = featurelist, bedfiles = bedfiles, maxmismatches = maxmismatches)
            #getbamcountsthr(allcounts, allcounts)
            #threads[currsample] = threading.Thread(target=getbamcountsthr, args=(allcounts,currsample,currbam, currsample,trnainfo, trnaloci, trnalist), kwargs = {'embllist' : embllist, 'featurelist' : featurelist, 'maxmismatches' : maxmismatches})
            #threads[currsample].start()
    endtime = time.time()
    #print >>sys.stderr, "time:" +str(endtime-starttime)
    if "countfile" not in argdict or argdict["countfile"] == "stdout":
        countfile = sys.stdout
    else:
        countfile = open(argdict["countfile"], "w")
    printcountfile(countfile, samples, allcounts,trnalist, trnaloci, featurelist, embllist, otherseqdict = otherseqdict,includebase = includebase)
    if genetypefile is not None:
        genetypeout = open(genetypefile, "w")
        printtypefile(genetypeout,samples, allcounts,trnalist, trnaloci, featurelist, embllist,otherseqdict = otherseqdict )
    #it's currently not used, but here is where I could count by amino acid or anticodon
    if typefile:
        trnacountfile = open(trnacountfilename, "w")
        for curramino in trnainfo.allaminos():
                print >>typefile, "AminoTotal_"+curramino+"\t"+"\t".join(str(aminocounts[currsample][curramino]) for currsample in samples)
        for currac in trnainfo.allanticodons():
                print >>typefile, "AnticodonTotal_"+currac+"\t"+"\t".join(str(anticodoncounts[currsample][currac]) for currsample in samples)


            
    if trnacountfilename is not None:
        #trnauniquefile = open(trnauniquefilename, "w")
        #printtrnacountfile()
        printtrnacountfile(trnacountfilename,samples,  allcounts, trnalist, trnaloci )
        printtrnaendfile(trnaendfilename,samples,  allcounts, trnalist, trnaloci )
       
        
    if trnauniquefilename is not None:
        
        #trnauniquefile = open(trnauniquefilename, "w")
        printtrnauniquecountcountfile(trnauniquefilename,samples,  allcounts, trnalist, trnaloci )
        #print >>trnauniquefile, "\t".join(currsample for currsample in samples)
        #for currfeat in trnalist:
        #    if max(trnacounts[currsample][currfeat.name] for currsample in samples) < minreads:
        #        continue
        #    print  >>trnauniquefile, currfeat.name+"\t"+"\t".join(str(trnacounts[currsample][currfeat.name]) for currsample in samples)
        #trnauniquefile.close()          
        pass


def oldmain(**argdict):
    trnauniquefilename = None
    argdict = defaultdict(lambda: None, argdict)
    includebase = argdict["nofrag"]
    fullpretrnasonly = argdict["onlyfullpretrnas"]
    trnatable = argdict["trnatable"]
    removepseudo = argdict["removepseudo"]
    ensemblgtf = argdict["ensemblgtf"]
    nomultimap = argdict["nomultimap"]
    maxmismatches = int(argdict["maxmismatches"])
    typefile = None
    sampledata = samplefile(argdict["samplefile"])
    bedfiles = list()
    if "trnauniquecounts" in argdict:
        trnauniquefilename = argdict["trnauniquecounts"]
    if "bedfile"  in argdict:
        bedfiles = argdict["bedfile"]
    trnalocifiles = list()
    if "trnaloci"  in argdict:
        trnalocifiles = argdict["trnaloci"]
    maturetrnas = list()
    if "maturetrnas" in argdict:
        maturetrnas = argdict["maturetrnas"]
        
    #trnalocifiles = argdict["trnaloci"]
    #maturetrnas=argdict["maturetrnas"]
    genetypefile = argdict["genetypefile"]
    trnacountfilename = argdict["trnacounts"]
    trnaendfilename = argdict["trnaends"]
    if "countfile" not in argdict or argdict["countfile"] == "stdout":
        countfile = sys.stdout
    else:
        countfile = open(argdict["countfile"], "w")
    
    trnacountfilename = argdict["trnacounts"]
    trnainfo = transcriptfile(trnatable)
    
    wholetrnas = dict()
    fivefrags = dict()
    threefrags = dict()
    trailerfrags = dict()
    otherfrags = dict()
    allfrags = dict()
    allowindels = False
    
    
    alltrnas = list()
    
    
    samplefiles = dict()
    
    
    samples = sampledata.getsamples()
    genetypes = dict()
    fullpretrnathreshold = 2
    #Grabbing all the features to count
    try:
        featurelist = list()
        trnaloci = list()
        for currfile in bedfiles:
            
            bedfeatures = list(readfeatures(currfile, removepseudo = removepseudo))
            for curr in bedfeatures:
                genetypes[curr.name] = os.path.basename(currfile)
                
            featurelist.extend(bedfeatures)
        trnalist = list()
        for currfile in trnalocifiles:
            trnaloci.extend(list(readbed(currfile)))
        for currfile in maturetrnas:
            trnalist.extend(list(readbed(currfile)))
        if ensemblgtf is not None:    
            embllist = list(readgtf(ensemblgtf, filterpsuedo = removepseudo))
        else:
            embllist = list()
    except IOError as e:
        print >>sys.stderr, e
        sys.exit()
    featcount = defaultdict(int)
    allfeats = featurelist+trnaloci+trnalist
    if len(set(curr.name for curr in allfeats)) < len(list(curr.name for curr in allfeats )):
        print >>sys.stderr, "Duplicate names in feature list"
    
    
    #featurelist = list(curr for curr in featurelist if curr.name == 'unknown20')
    alltrnas = list(curr.name for curr in featurelist)
    #print >>sys.stderr, "***"
    #setting up all the feature count dictionaries
    counts = defaultdict(lambda: defaultdict(int))
    trnacounts = defaultdict(lambda: defaultdict(int))
    trnawholecounts = defaultdict(lambda: defaultdict(int))
    trnafivecounts = defaultdict(lambda: defaultdict(int))
    trnathreecounts = defaultdict(lambda: defaultdict(int))
    trnalocuscounts = defaultdict(lambda: defaultdict(int))
    trnalocustrailercounts = defaultdict(lambda: defaultdict(int))
    partialtrnalocuscounts = defaultdict(lambda: defaultdict(int))
    fulltrnalocuscounts  = defaultdict(lambda: defaultdict(int))
    trnauniquecounts = defaultdict(lambda: defaultdict(int))
    aminocounts  = defaultdict(lambda: defaultdict(int))
    anticodoncounts =  defaultdict(lambda: defaultdict(int))                                         
    
    #how much a pre-tRNA must extend off the end
    minpretrnaextend = 5
    #minimum mapq
    minmapq = 0
    if nomultimap:
        minmapq = 2
    #minimum number of reads for a feature to be reported
    minreads = 5
    for currsample in samples:
        
        currbam = sampledata.getbam(currsample)
        #print >>sys.stderr, currsample
        #doing this thing here why I only index the bamfile if the if the index file isn't there or is older than the map file
        try:
            if not os.path.isfile(currbam+".bai") or os.path.getmtime(currbam+".bai") < os.path.getmtime(currbam):
                pysam.index(""+currbam)
            bamfile = pysam.Samfile(""+currbam, "rb" )  
        except IOError as ( strerror):
            print >>sys.stderr, strerror
            sys.exit(1)
        except pysam.utils.SamtoolsError:
            print >>sys.stderr, "Can not index "+currbam
            print >>sys.stderr, "Exiting..."
            sys.exit(1)
            
        
        for currfeat in featurelist:
            #try catch is to account for weird chromosomes and the like that aren't in the genome
            #means that if I can't find a feature, I record no counts for it rather than bailing
            try:
                for currread in getbamrange(bamfile, currfeat, singleonly = nomultimap, maxmismatches = maxmismatches,allowindels = allowindels):
                    if currfeat.coverage(currread) > 10:
                        counts[currsample][currfeat.name] += 1
            except ValueError:
                pass
         
        for genename, featset in itertools.groupby(embllist,lambda x: x.data["genename"]):
            
            #pass 
            try:
                allreads =set()
                for currfeat in list(featset):
                    for currread in getbamrangeshort(bamfile, currfeat, singleonly = nomultimap, maxmismatches = maxmismatches,allowindels = allowindels, skiptags = True):
                        #print >>sys.stderr, "**"+currread.name
                        #continue
                        if currfeat.coverage(currread) > 10:
                            counts[currsample][genename] += 1 
                            genetypes[genename] = currfeat.data["biotype"]
                            #print >>sys.stderr, currfeat.bedstring()
            except ValueError:
                pass
        for currfeat in trnaloci:
            #print >>sys.stderr,  currfeat.bedstring()
            #print >>sys.stderr,  currfeat.getdownstream(30).bedstring()
            for currread in getbamrangeshort(bamfile, currfeat.addmargin(30), singleonly = nomultimap, maxmismatches = maxmismatches,allowindels = allowindels, skiptags = True):
                #gotta be more than 5 bases off one end to be a true pre-tRNA
                #might want to shove these to the real tRNA at some point, but they are for now just ignored

                if currfeat.coverage(currread) > 10 and (currread.start + minpretrnaextend <= currfeat.start or currread.end - minpretrnaextend >= currfeat.end):
                    trnalocuscounts[currsample][currfeat.name] += 1
                    if currread.start + fullpretrnathreshold <  currfeat.start and currread.end - fullpretrnathreshold + 3 >  currfeat.end:
                        fulltrnalocuscounts[currsample][currfeat.name] += 1
                    else:
                        partialtrnalocuscounts[currsample][currfeat.name] += 1
                elif currfeat.getdownstream(30).coverage(currread) > 10:  #need the elif otherwise fragments that include trailer get in there
                    trnalocuscounts[currsample][currfeat.name] += 1
                    trnalocustrailercounts[currsample][currfeat.name] += 1
                else:
                    #print >>sys.stderr,  currfeat.getdownstream(30).coverage(currread)
                    pass
        
        for currfeat in trnalist:
            for currread in getbamrange(bamfile, currfeat, singleonly = nomultimap, maxmismatches = maxmismatches,allowindels = allowindels):

                if not currfeat.strand == currread.strand:
                    continue
                if not currfeat.coverage(currread) > 10:
                    continue
                curramino = trnainfo.getamino(currfeat.name)
                curranticodon = trnainfo.getanticodon(currfeat.name)
                trnacounts[currsample][currfeat.name] += 1
                    
                fragtype = getfragtype(currfeat, currread)

                if fragtype == "Whole":
                    trnawholecounts[currsample][currfeat.name] += 1
                elif fragtype == "Fiveprime":
                    trnafivecounts[currsample][currfeat.name] += 1
                elif fragtype == "Threeprime":
                    trnathreecounts[currsample][currfeat.name] += 1
                if isuniqueaminomapping(currread):
                    trnauniquecounts[currsample][currfeat.name] += 1
                if not isuniqueaminomapping(currread):
                    pass
                elif not isuniqueacmapping(currread):
                    aminocounts[currsample][curramino] += 1
                else:
                    aminocounts[currsample][curramino] += 1
                    anticodoncounts[currsample][curranticodon] += 1
                  
                            
    
    print >>countfile, "\t".join(samples)
    
    trnanames = set()
    for currfeat in trnalist:
        if max(itertools.chain((trnacounts[currsample][currfeat.name] for currsample in samples), [0])) < minreads:
            continue
        if includebase:
            print >>countfile, currfeat.name+"\t"+"\t".join(str(trnacounts[currsample][currfeat.name]) for currsample in samples)
        else:
            print >>countfile, currfeat.name+"_wholecounts\t"+"\t".join(str(trnawholecounts[currsample][currfeat.name]) for currsample in samples)
            print >>countfile, currfeat.name+"_fiveprime\t"+"\t".join(str(trnafivecounts[currsample][currfeat.name]) for currsample in samples)
            print >>countfile, currfeat.name+"_threeprime\t"+"\t".join(str(trnathreecounts[currsample][currfeat.name]) for currsample in samples)
            print >>countfile, currfeat.name+"_other\t"+"\t".join(str(trnacounts[currsample][currfeat.name] - (trnathreecounts[currsample][currfeat.name] + trnafivecounts[currsample][currfeat.name] + trnawholecounts[currsample][currfeat.name])) for currsample in samples)
        
        
    
    for currfeat in trnaloci:
        if max(itertools.chain((trnalocuscounts[currsample][currfeat.name] for currsample in samples),[0])) < minreads:
            continue
        if includebase:
            print >>countfile, currfeat.name+"\t"+"\t".join(str(trnalocuscounts[currsample][currfeat.name]) for currsample in samples)
        else:
            print >>countfile, currfeat.name+"_wholeprecounts\t"+"\t".join(str(fulltrnalocuscounts[currsample][currfeat.name]) for currsample in samples)
            print >>countfile, currfeat.name+"_partialprecounts\t"+"\t".join(str(partialtrnalocuscounts[currsample][currfeat.name]) for currsample in samples)
            print >>countfile, currfeat.name+"_trailercounts\t"+"\t".join(str(trnalocustrailercounts[currsample][currfeat.name]) for currsample in samples)
    
    #it's currently not used, but here is where I could count by amino acid or anticodon
    if typefile:
        for curramino in trnainfo.allaminos():
                print >>typefile, "AminoTotal_"+curramino+"\t"+"\t".join(str(aminocounts[currsample][curramino]) for currsample in samples)
        for currac in trnainfo.allanticodons():
                print >>typefile, "AnticodonTotal_"+currac+"\t"+"\t".join(str(anticodoncounts[currsample][currac]) for currsample in samples)
    if genetypefile is not None:
        genetypeout = open(genetypefile, "w")
    for currfeat in featurelist :
        if currfeat.name in trnanames:
            continue
        trnanames.add(currfeat.name)
        if max(counts[currsample][currfeat.name] for currsample in samples) > minreads:
            print >>countfile, currfeat.name+"\t"+"\t".join(str(counts[currsample][currfeat.name]) for currsample in samples)
            if genetypefile is not None:
                print >>genetypeout, currfeat.name+"\t"+genetypes[currfeat.name]   
    
    if genetypefile is not None:
        
        for currfeat in trnaloci:
            print >>genetypeout, currfeat.name+"_wholeprecounts"+"\t"+"trna_wholeprecounts"
            print >>genetypeout, currfeat.name+"_partialprecounts"+"\t"+"trna_partialprecounts"
            print >>genetypeout, currfeat.name+"_trailercounts"+"\t"+"trna_trailercounts"
            print >>genetypeout, currfeat.name+""+"\t"+"tRNA_locus"
        for currfeat in trnalist:
            print >>genetypeout, currfeat.name+"_wholecounts"+"\t"+"trna_wholecounts"
            print >>genetypeout, currfeat.name+"_fiveprime"+"\t"+"trna_fiveprime"
            print >>genetypeout, currfeat.name+"_threeprime"+"\t"+"trna_threeprime"
            print >>genetypeout, currfeat.name+"_other"+"\t"+"trna_other"
            print >>genetypeout, currfeat.name+""+"\t"+"tRNA"
    
    for currfeat in embllist:
        genename = currfeat.data['genename']
        if genename in trnanames:
            continue
        trnanames.add(genename)
        if genename is None:
            print >>sys.stderr, currfeat.name
            continue
        if max(counts[currsample][genename] for currsample in samples) > minreads:
            print >>countfile, genename+"\t"+"\t".join(str(counts[currsample][genename]) for currsample in samples)
            if genetypefile is not None:
                print >>genetypeout, genename+"\t"+genetypes[genename]          
            


                    
    if genetypefile is not None:
        genetypeout.close()    
            
            
            
    if trnacountfilename is not None:
        trnacountfile = open(trnacountfilename, "w")
        print >>trnacountfile, "\t".join(currsample for currsample in samples)
        for currfeat in trnaloci:
            if max(trnalocuscounts[currsample][currfeat.name] for currsample in samples) < minreads:
                continue
            print >>trnacountfile, currfeat.name+"\t"+"\t".join(str(trnalocuscounts[currsample][currfeat.name]) for currsample in samples)
        for currfeat in trnalist:
            if max(trnacounts[currsample][currfeat.name] for currsample in samples) < minreads:
                continue
            print  >>trnacountfile, currfeat.name+"\t"+"\t".join(str(trnacounts[currsample][currfeat.name]) for currsample in samples)
        trnacountfile.close()           
        
        
    if trnauniquefilename is not None:
        trnauniquefile = open(trnauniquefilename, "w")
        print >>trnauniquefile, "\t".join(currsample for currsample in samples)
        for currfeat in trnalist:
            if max(trnacounts[currsample][currfeat.name] for currsample in samples) < minreads:
                continue
            print  >>trnauniquefile, currfeat.name+"\t"+"\t".join(str(trnacounts[currsample][currfeat.name]) for currsample in samples)
        trnauniquefile.close()          

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
    parser.add_argument('--samplefile',
                       help='Sample file in format')
    parser.add_argument('--bedfile',  nargs='+', default=list(),
                       help='bed file with non-tRNA features')
    parser.add_argument('--gtffile',  nargs='+', default=list(),
                       help='gtf file with non-tRNA features')
    parser.add_argument('--ensemblgtf',
                       help='ensembl gtf file with tRNA features')
    parser.add_argument('--trnaloci',  nargs='+', default=list(),
                       help='bed file with tRNA features')
    parser.add_argument('--maturetrnas',  nargs='+', default=list(),
                       help='bed file with mature tRNA features')
    parser.add_argument('--onlyfullpretrnas', action="store_true", default=False,
                       help='only include full pre-trnas')
    parser.add_argument('--trnatable',
                       help='table of tRNA features')
    parser.add_argument('--removepseudo', action="store_true", default=False,
                       help='remove psuedogenes from ensembl GTFs')
    parser.add_argument('--genetypefile',
                       help='Output file with gene types')
    parser.add_argument('--trnacounts',
                       help='Output file with just trna gene counts')
    parser.add_argument('--nofrag', action="store_true", default=False,
                       help='disable fragment determination')
    parser.add_argument('--nomultimap', action="store_true", default=False,
                       help='do not count multiply mapped reads')
    parser.add_argument('--maxmismatches', default=None,
                       help='Set maximum number of allowable mismatches')
    
    
    args = parser.parse_args()
        
    #main(samplefile=args.samplefile, bedfile=args.bedfile, gtffile=args.bedfile, ensemblgtf=args.ensemblgtf, trnaloci=args.trnaloci, onlyfullpretrnas=args.onlyfullpretrnas,removepseudo=args.removepseudo,genetypefile=args.genetypefile,trnacounts=args.trnacounts,maturetrnas=args.maturetrnas,nofrag=args.nofrag,nomultimap=args.nomultimap,maxmismatches=args.maxmismatches)
    argvars = vars(args)
    #argvars["countfile"] = "stdout"
    testmain(**argvars)
        