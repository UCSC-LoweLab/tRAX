#!/usr/bin/env python

import pysam
import sys
import argparse
import os.path
from collections import defaultdict
from trnasequtils import *
import itertools
from multiprocessing import Process, Queue, Pool
import time


def getchromdict(features):
    chromdict = defaultdict(list)
    for curr in features:
        chromdict[curr.chrom].append(curr)
    return chromdict
class counttypes:
    def __init__(self, samplename, bamfile, trnas = list(), trnaloci = list(), emblgenes = list(), otherfeats = list()):
        self.samplename = samplename
        self.bamfile = bamfile
        self.trnas = trnas
        self.trnaloci = trnaloci
        self.emblgenes = emblgenes
        self.otherfeats = otherfeats
        self.emblbiotypes = set()
        self.aminos = set()
        self.bedtypes = set()
        self.extraseqtypes = set()
        self.anticodons = set()
        
        self.embltypecounts =defaultdict(int)
        self.bedtypecounts =defaultdict(int)
        self.trnafragtypes =defaultdict(int)
        self.trnafragtypes =defaultdict(int)
        self.totalreads = 0
        self.trnareads = 0
        self.otherreads = 0
        self.readlengths = defaultdict(int)
        self.trnareadlengths = defaultdict(int)
        self.trnacounts = defaultdict(int)
        self.aminocounts = defaultdict(int)
        self.anticodoncounts = defaultdict(int)
        self.indelreads = defaultdict(int)
        self.trnaanticounts = defaultdict(int)
        self.trnaemblcounts = defaultdict(int)
        self.pretrnareadlengths = defaultdict(int)
        self.trnalocuscounts = defaultdict(int)
        self.partiallocuscounts = defaultdict(int)
        self.fulllocuscounts = defaultdict(int)
        self.extraseqcounts = defaultdict(int)
        self.trnaantilocuscounts  = defaultdict(int)

    def addsamplecounts(self):
        self.totalreads += 1
    def addreadlengths(self, length):
        self.readlengths[length] += 1
    def addtrnareadlengths(self, length):
        self.trnareadlengths[length] += 1
    def addpretrnareadlengths(self, length):
        self.pretrnareadlengths[length] += 1
    def addpartiallocuscounts(self, currbed):
        self.partiallocuscounts[currbed] += 1     
    def addtrnasamplecounts(self):
        self.trnareads += 1
    def addtrnacounts(self, currbed):
        self.trnacounts[currbed] += 1
    def addtrnaantilocuscounts(self, currbed):
        self.trnaantilocuscounts[currbed] += 1
    def addtrnalocuscounts(self, currbed):
        self.trnalocuscounts[currbed] += 1
    def addfulllocuscounts(self, currbed):
        self.fulllocuscounts[currbed] += 1
    def addaminocounts(self, curramino):
        self.aminos.add(curramino)
        self.aminocounts[curramino] += 1
    def addanticodoncounts(self, curranticodon):
        self.anticodons.add(curranticodon)
        self.anticodoncounts[curranticodon] += 1
    def addindelreads(self, curramino):
        self.indelreads[curramino] += 1
    def addotherreads(self):
        self.otherreads += 1
    def addtrnaantisense(self, currbed):
        self.trnaanticounts[currbed] += 1
    def addemblcounts(self, currtype):
        self.emblbiotypes.add(currtype)
        self.embltypecounts[currtype] += 1     
    def addbedcounts(self, genetype):
        self.bedtypes.add(genetype)
        self.bedtypecounts[genetype] += 1
    def addextracounts(self, genetype):
        self.extraseqtypes.add(genetype)
        self.extraseqcounts[genetype] += 1


def counttypereads(bamfile, samplename,trnainfo, trnaloci, trnalist,maturenames,featurelist = dict(), otherseqlist = list(), embllist = list(),nomultimap = False, allowindels = True, maxmismatches = None, bamnofeature = False, countfrags = False):
    
    bedlist = featurelist.keys()
    readtypecounts = counttypes(samplename, bamfile, trnas = trnalist, trnaloci = trnaloci, emblgenes = embllist, otherfeats = bedlist)
    mitochrom = None
    fullpretrnathreshold = 2
    minpretrnaextend = 5
    ncrnaorder = defaultdict(int)
    currbam = bamfile
    dumpotherreads = True

    for i, curr in enumerate(reversed(list(["snoRNA","miRNA", "rRNA","snRNA","misc_RNA","lincRNA", "protein_coding"]))):
        ncrnaorder[curr] = i + 1

        
    try:
        #print >>sys.stderr, currbam
        
        
        if not os.path.isfile(currbam+".bai") or  os.path.getmtime(currbam+".bai") < os.path.getmtime(currbam):
            pysam.index(""+currbam)
        bamfile = pysam.Samfile(""+currbam, "rb" )
        if bamnofeature:
            outname = os.path.splitext(currbam)[0]+"_nofeat.bam"
            outbamnofeature =  pysam.Samfile( outname, "wb", template =  bamfile)
    except IOError as ( strerror):
        print >>sys.stderr, strerror
        sys.exit()
    #continue #point0
    #print >>sys.stderr, "**||"+currbam

    for i, currread in enumerate(getbam(bamfile, primaryonly = True)):

        isindel = False
        readlength = currread.getlength()
        gotread = False
        #continue #point1
        readtypecounts.addsamplecounts()
        if currread.hasindel():
            readtypecounts.addindelreads(readlength)
            isindel = True
            #continue

        else:
            pass
            #continue
        readtypecounts.addreadlengths(readlength)
        #readlengths[currsample][readlength] += 1
        #continue #point2
        #if currread.name == "NB501427:473:H3YJ2BGXG:3:11512:1635:4586":
        #    print >>sys.stderr, "**foundread"
        for currbed in trnaloci:
            for currfeat in trnaloci[currbed].getbin(currread):
                expandfeat = currfeat.addmargin(30)
                #if currread.name == "NB501427:473:H3YJ2BGXG:3:11512:1635:4586"  and currfeat.name == "tRNA-Val-AAC-5-1":
                #    print >>sys.stderr, "**trnaread:" +str(currfeat.coverage(currread))
                #    print >>sys.stderr, "**trnaread:" +str(currfeat.coverage(currread))
                if currfeat.coverage(currread) > 10: # and (currread.start + minpretrnaextend <= currfeat.start or currread.end - minpretrnaextend >= currfeat.end):
                    #if currread.name == "NB501427:473:H3YJ2BGXG:3:11512:1635:4586":
                    #    print >>sys.stderr, "**foundread:" +str(currfeat.name)
                    if currfeat.strand != currread.strand:
                        readtypecounts.addantilocuscounts(currbed)
                        break
                    if (currread.start + minpretrnaextend <= currfeat.start or currread.end - minpretrnaextend >= currfeat.end):
                        pass
                    readtypecounts.addpretrnareadlengths(readlength)
                    readtypecounts.addtrnalocuscounts(currbed)
                    if currread.start + fullpretrnathreshold <  currfeat.start and currread.end - fullpretrnathreshold + 3 >  currfeat.end:
                        #fulltrnalocuscounts[currsample][currbed] += 1
                        readtypecounts.addfulllocuscounts(currbed)
                    else:# currread.start + fullpretrnathreshold <  currfeat.start or currread.end - fullpretrnathreshold +3 >  currfeat.end:
                        #partialtrnalocuscounts[currsample][currbed] += 1
                        readtypecounts.addpartiallocuscounts(currbed)
                        #print >>sys.stderr, "***"
                    gotread = True
                    break
                if currfeat.getdownstream(30).coverage(currread) > 10:
                    readtypecounts.addtrnaantisense(currbed)
                    #readtypecounts.addpretrnareadlengths(readlength)
                    #print >>sys.stderr, currfeat.bedstring()
                    gotread = True
                    break
                elif expandfeat.antisense().coverage(currread) > 5:
                    #trnaantisense[currsample][currbed] += 1
                    readtypecounts.addtrnaantisense(currbed)
                    gotread = True
                    break
        #if currread.name == "NB501427:473:H3YJ2BGXG:3:11512:1635:4586":
        #    print >>sys.stderr, "**foundread2"            
        if gotread: 
            continue
        #continue #point3

        for currbed in trnalist:
            if currread.chrom in maturenames[currbed]:
                currfeat = maturenames[currbed][currread.chrom]
                if currread.strand == "+":
                    readtypecounts.addtrnareadlengths(readlength)
                    readtypecounts.addtrnasamplecounts()
                    readtypecounts.addtrnacounts(currbed)
                    readtypecounts.addaminocounts(trnainfo.getamino(currfeat.name))
                    readtypecounts.addanticodoncounts(trnainfo.getanticodon(currfeat.name))
                    

                    #aminos.add(trnainfo.getamino(currfeat.name))
                    if countfrags:
                        fragtype = getfragtype(currfeat, currread)
                        
                        if fragtype == "Whole":
                            
                            trnawholecounts[currsample][currbed] += 1
                        elif fragtype == "Fiveprime":
                            trnafivecounts[currsample][currbed] += 1
                        elif fragtype == "Threeprime":
                            trnathreecounts[currsample][currbed] += 1
                        elif fragtype == "Trailer":
                            trnatrailercounts[currsample][currbed] += 1
                    gotread = True
                    break
                        #print >>sys.stderr, str(currread.start - currfeat.start)+"-"+str(currread.end - currfeat.start)  
                        #print >>sys.stderr, str(currfeat.start - currfeat.start)+"-"+str(currfeat.end - currfeat.start)
                        #print >>sys.stderr, "****"
                elif currfeat.antisense().coverage(currread) > 10:
                    readtypecounts.addtrnaantisense(currbed)
                    #trnaantisense[currsample][currbed] += 1
                    gotread = True
                    break
        if gotread: 
            continue
        #continue #point4
        if embllist is not None:
            currtype = None
            for currfeat in embllist.getbin(currread):
                if currfeat.coverage(currread) > 10: 
                    if currfeat.data["biotype"] == "processed_transcript":
                        #print >>sys.stderr, currfeat.bedstring()
                        
                        pass

                    if currtype is None or ncrnaorder[currfeat.data["biotype"]] > ncrnaorder[currtype]:
                        currtype= currfeat.data["biotype"]
                        #if mitochrom == currread.chrom:
                            #currtype = "mito"+currtype
                    
                    
                    
            if currtype is not None:
                readtypecounts.addemblcounts(currtype)
                #emblcounts[currsample][currtype] += 1
                #emblbiotypes.add(currtype)
                gotread = True
                    #print >>sys.stderr, currbam +":"+ currbed
        if gotread: 
            continue
        #continue #point5
        #print >>sys.stderr, "**||"
        
        for currbed in bedlist:

            #if currread.name == "SRR10038183.1660151":
            #    print >>sys.stderr, "||"+currbed
            #    
            #    print >>sys.stderr, list(featurelist[currbed].getbin(currread))
            #    print >>sys.stderr, list(featurelist[currbed].getfeatbin("12-qE-23911.2"))
            #    print >>sys.stderr, list(featurelist[currbed].getbinnums(currread))
            
            for currfeat in featurelist[currbed].getbin(currread):
                
                if currfeat.coverage(currread) > 10:
                    #print >>sys.stderr, currbam +":"+ currbed
                    readtypecounts.addbedcounts(currbed)
                    #counts[currsample][currbed] += 1
                    gotread = True
                    break
                    #print >>sys.stderr, currbam +":"+ currbed
        if gotread: 
            continue
        for currbed in otherseqlist:
            #print >>sys.stderr, len(list(otherseqlist[currbed].getbin(currread)))
            #print >>sys.stderr, "**"+ otherseqlist[currbed].keys()[0]
            #print >>sys.stderr, "||"+ currfeat.chrom
            for currfeat in otherseqlist[currbed][currread.chrom]:
                #print >>sys.stderr, "**"+currbed
                if currfeat.coverage(currread) > 10:
                    
                    readtypecounts.addextracounts(currbed)
                    gotread = True
                    break 
                    #print >>sys.stderr, currbam +":"+ currbed
        #print >>sys.stderr, "**||"
        if gotread: 
            continue
        readtypecounts.addotherreads()
        if not gotread and embllist is not None and mitochrom == currread.chrom:
            currtype = "Mitochondrial_other"
            readtypecounts.addemblcounts(currtype)
            #emblbiotypes.add(currtype)
        if not gotread and bamnofeature:
            outbamnofeature.write(currread.bamline)
    return readtypecounts



def printtypefile(countfile,samples, sampledata,allcounts,trnalist, trnaloci, bedtypes, emblbiotypes, sizefactor,extraseqtypes = set(),countfrags = False, combinereps = True):

    def sumsamples(countdict,sampledata, repname, currfeat = None, sizefactors = defaultdict(lambda: 1)):
        if currfeat is None: #To account for the "other" counts, which don't have a feature
            return sum(countdict[currsample]/sizefactors[currsample] for currsample in sampledata.getrepsamples(repname))
        else:
            return sum(countdict[currsample][currfeat]/sizefactors[currsample] for currsample in sampledata.getrepsamples(repname))
  
     
    
    if combinereps:
        replicates = list(sampledata.allreplicates())
        print >>countfile, "\t".join(replicates)
        #print >>sys.stderr, allcounts[sampledata.getrepsamples(replicates[0])[0]].embltypecounts 

        #print  >>countfile, "other"+"\t"+"\t".join(str(sumsamples(othercounts,sampledata,currrep, sizefactors = sizefactor)) for currrep in replicates)
        print  >>countfile, "other"+"\t"+"\t".join(str(sum(allcounts[currsample].otherreads/sizefactor[currsample] for currsample in sampledata.getrepsamples(currrep))) for currrep in replicates)
        for currbed in bedtypes:
            print  >>countfile, os.path.basename(currbed).split(".")[0]+"\t"+"\t".join(str(sum(allcounts[currsample].bedtypecounts[currbed]/sizefactor[currsample] for currsample in sampledata.getrepsamples(currrep))) for currrep in replicates)
        #print  >>countfile, "other"+"\t"+"\t".join(str(sumsamples(othercounts,sampledata,currrep, sizefactors = sizefactor)) for currrep in replicates)
        #sys.exit()
        
        
        for currname in extraseqtypes:
            #print >>sys.stderr, currname
            print >>countfile, currname+"_seq\t"+"\t".join(str(sum(allcounts[currsample].extraseqcounts[currname]/sizefactor[currsample] for currsample in sampledata.getrepsamples(currrep))) for currrep in replicates)
 
        biotypefirst = ['snoRNA','snRNA','scaRNA','sRNA','miRNA']         
        biotypelast = ['Mt_rRNA','Mt_tRNA','rRNA']
        otherbiotypes = list(set(emblbiotypes) - (set(biotypefirst) | set(biotypelast)))
        biotypeorder = biotypefirst + otherbiotypes + biotypelast
        #print >>sys.stderr, biotypeorder

        for currbiotype in biotypeorder:
            #print  >>countfile, currbiotype+"\t"+"\t".join(str(sumsamples(emblcounts,sampledata,currrep, currbiotype, sizefactors = sizefactor)) for currrep in replicates)
            print  >>countfile, currbiotype+"\t"+"\t".join(str(sum(allcounts[currsample].embltypecounts[currbiotype]/sizefactor[currsample] for currsample in sampledata.getrepsamples(currrep))) for currrep in replicates)
               

        for currbed in trnaloci:
            if countfrags:
                pass
                #print  >>countfile, "pretRNA_full\t"+"\t".join(str(sumsamples(fulltrnalocuscounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
                #print  >>countfile, "pretRNA_partial\t"+"\t".join(str(sumsamples(partialtrnalocuscounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
                #print  >>countfile, "pretRNA_trailer\t"+"\t".join(str(sumsamples(trnalocustrailercounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
            else:
                print  >>countfile, "pretRNA_amtosemse\t"+"\t".join(str(sum(allcounts[currsample].trnaantilocuscounts[currbed]/sizefactor[currsample] for currsample in sampledata.getrepsamples(currrep))) for currrep in replicates)
                print  >>countfile, "pretRNA\t"+"\t".join(str(sum(allcounts[currsample].trnalocuscounts[currbed]/sizefactor[currsample] for currsample in sampledata.getrepsamples(currrep))) for currrep in replicates)

        for currbed in trnalist:
            
            if countfrags:
                pass
                #print  >>countfile, "tRNA_wholecounts\t"+"\t".join(str(sumsamples(trnawholecounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
                #print  >>countfile, "tRNA_fiveprime\t"+"\t".join(str(sumsamples(trnafivecounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
                #print  >>countfile, "tRNA_threeprime\t"+"\t".join(str(sumsamples(trnathreecounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
                #print  >>countfile, "tRNA_other\t"+"\t".join(str(sumsamples(trnacounts,sampledata,currrep, currbed, sizefactors = sizefactor) - (sumsamples(trnafivecounts,sampledata,currrep, currbed, sizefactors = sizefactor) + sumsamples(trnathreecounts,sampledata,currrep, currbed, sizefactors = sizefactor) + sumsamples(trnawholecounts,sampledata,currrep, currbed, sizefactors = sizefactor))) for currrep in replicates)
                #print  >>countfile, "tRNA_antisense\t"+"\t".join(str(sumsamples(trnaantisense,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
            else:
                print  >>countfile, "tRNA_antisense\t"+"\t".join(str(sum(allcounts[currsample].trnaanticounts[currbed]/sizefactor[currsample] for currsample in sampledata.getrepsamples(currrep))) for currrep in replicates)
                print  >>countfile, "tRNA\t"+"\t".join(str(sum(allcounts[currsample].trnacounts[currbed]/sizefactor[currsample] for currsample in sampledata.getrepsamples(currrep))) for currrep in replicates)
                

                #print  >>countfile, "tRNA\t"+"\t".join(str(sumsamples(allcounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
            
        
    else:
        pass
        '''
        print  >>countfile, "\t".join(samples)
        




        print  >>countfile, "other"+"\t"+"\t".join(str(othercounts[currsample]/sizefactor[currsample]) for currsample in samples)
        
        for currbed in bedlist:
            print  >>countfile, os.path.basename(currbed).split(".")[0]+"\t"+"\t".join(str(counts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
        biotypefirst = ['snoRNA','snRNA','scaRNA','sRNA','miRNA']         
        biotypelast = ['Mt_rRNA','Mt_tRNA','rRNA']
        otherbiotypes = list(set(sampledata.emblbiotypes) - (set(biotypefirst) | set(biotypelast)))
        biotypeorder = biotypefirst + otherbiotypes + biotypelast
        #print >>sys.stderr, biotypeorder
        #print >>sys.stderr, "****"
        for currbiotype in biotypeorder:
            print  >>countfile, currbiotype+"\t"+"\t".join(str(emblcounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
        for currbed in locilist:
            if countfrags:
                print  >>countfile, "pretRNA_full\t"+"\t".join(str(fulltrnalocuscounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
                print  >>countfile, "pretRNA_partial\t"+"\t".join(str(partialtrnalocuscounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
                print  >>countfile, "pretRNA_trailer\t"+"\t".join(str(trnalocustrailercounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
            else:
                print  >>countfile, "pretRNA\t"+"\t".join(str(trnalocuscounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
        for currbed in trnalist:
            
            if countfrags:
                print  >>countfile, "tRNA_wholecounts\t"+"\t".join(str(trnawholecounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
                print  >>countfile, "tRNA_fiveprime\t"+"\t".join(str(trnafivecounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
                print  >>countfile, "tRNA_threeprime\t"+"\t".join(str(trnathreecounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
                print  >>countfile, "tRNA_other\t"+"\t".join(str((trnacounts[currsample][currbed] - (trnathreecounts[currsample][currbed] + trnafivecounts[currsample][currbed] + trnawholecounts[currsample][currbed]))/sizefactor[currsample]) for currsample in samples)
                print  >>countfile, "tRNA_antisense\t"+"\t".join(str(trnaantisense[currsample][currbed]/sizefactor[currsample]) for currsample in samples)

                
                
            else:
                print  >>countfile, "tRNA\t"+"\t".join(str(trnacounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
                
        '''        
def printrealcounts(countfile,samples, sampledata,allcounts,trnalist, trnaloci, bedtypes, emblbiotypes,extraseqtypes = set()):

    biotypefirst = ['snoRNA','snRNA','scaRNA','sRNA','miRNA']         
    biotypelast = ['Mt_rRNA','Mt_tRNA','rRNA']
    otherbiotypes = list(set(emblbiotypes) - (set(biotypefirst) | set(biotypelast)))
    biotypeorder = biotypefirst + otherbiotypes + biotypelast
        
    replicates = list(sampledata.allreplicates())
    print >>countfile, "\t".join(samples)


    #print >>sys.stderr, allcounts[sampledata.getrepsamples(replicates[0])[0]].embltypecounts 


    #print  >>countfile, "other"+"\t"+"\t".join(str(sumsamples(othercounts,sampledata,currrep, sizefactors = sizefactor)) for currrep in replicates)

    
    print  >>countfile, "other"+"\t"+"\t".join(str(allcounts[currsample].otherreads) for currsample in samples)
    for currbed in bedtypes:
        print  >>countfile, os.path.basename(currbed)+"\t"+"\t".join(str(allcounts[currsample].bedtypecounts[currbed]) for currsample in samples)
    for currname in extraseqtypes:
         print >>countfile, currname+"_seq\t"+"\t".join(str(allcounts[currsample].extraseqcounts[currname]) for currsample in samples)

    for currbiotype in reversed(biotypeorder):
        print  >>countfile, currbiotype+"\t"+"\t".join(str(allcounts[currsample].embltypecounts[currbiotype]) for currsample in samples)
        
    for currbed in trnaloci:
 
        print  >>countfile, "pretRNA\t"+"\t".join(str(allcounts[currsample].trnalocuscounts[currbed]) for currsample in samples)
    for currbed in trnalist:     
        print  >>countfile, "tRNA_antisense\t"+"\t".join(str(allcounts[currsample].trnaanticounts[currbed]) for currsample in samples)
        print  >>countfile, "tRNA\t"+"\t".join(str(allcounts[currsample].trnacounts[currbed]) for currsample in samples)
        

def printaminocounts(trnaaminofilename, sampledata,allcounts, sizefactor):
    #print >>sys.stderr, trnaaminocounts
    replicates = list(sampledata.allreplicates())
    trnaaminofile = open(trnaaminofilename, "w")
    aminos = allaminos
    #print >>sys.stderr, aminos
    print  >>trnaaminofile, "\t".join(replicates)
    for curramino in aminos:
        #print >>sys.stderr, curramino
        print >>trnaaminofile, curramino+"\t"+"\t".join(str(sum(allcounts[currsample].aminocounts[curramino]/sizefactor[currsample] for currsample in sampledata.getrepsamples(currrep))) for currrep in replicates)

def printanticodoncounts(trnaanticodonfilename, sampledata,allcounts, sizefactor):
    replicates = list(sampledata.allreplicates())
    trnaanticodonfile = open(trnaanticodonfilename, "w")
    anticodons = set(itertools.chain.from_iterable(allcounts[currsample].anticodons for currsample in sampledata.getsamples()))
    print  >>trnaanticodonfile, "\t".join(replicates)
    for curranticodon in anticodons:
        print >>trnaanticodonfile, curranticodon+"\t"+"\t".join(str(sum(allcounts[currsample].anticodoncounts[curranticodon]/sizefactor[currsample] for currsample in sampledata.getrepsamples(currrep))) for currrep in replicates)


def printtrnanormfile(samples, allcounts):         
    #samples trnasamplecounts.keys()
    trnanormfile = open(trnanormfile, "w")
    mean = 1.*sum(trnasamplecounts.values())/len(trnasamplecounts.values())
    print >>trnanormfile, "\t".join(samples)
    print >>trnanormfile, "\t".join(str(trnasamplecounts[currsample]/mean) for currsample in samples)
def printallreadsnormfile(samples, allcounts):
    allreadsnormfile = open(allreadsnormfile, "w")
    mean = 1.*sum(totalsamplecounts.values())/len(totalsamplecounts.values())
    print >>allreadsnormfile,"\t".join(samples)
    print >>allreadsnormfile,"\t".join(str(totalsamplecounts[currsample]/mean) for currsample in samples)
def printlengthfile(readlengthfile, samples,allcounts):
    readlengthfile = open(readlengthfile, "w")
    print >>readlengthfile, "Length\tSample\tother\ttrnas\tpretrnas"
    for currsample in samples:
        for curr in range(0,max(allcounts[currsample].readlengths.keys())+1):
            othercount = allcounts[currsample].trnareadlengths[curr] + allcounts[currsample].pretrnareadlengths[curr]
            print >>readlengthfile, str(curr)+"\t"+currsample+"\t"+str(allcounts[currsample].readlengths[curr] - othercount)+"\t"+str(allcounts[currsample].trnareadlengths[curr]) +"\t"+str(allcounts[currsample].pretrnareadlengths[curr])
        
def counttypereadsqueue(countqueue,currsample, *args, **kwargs):
    countqueue.put([currsample,counttypereads(*args, **kwargs)])

def counttypereadspool(args):
    return counttypereads(*args[0], **args[1])
    
def compressargs( *args, **kwargs):
    return tuple([args, kwargs])
    
def testmain(**argdict):
    argdict = defaultdict(lambda: None, argdict)
    countfrags = argdict["countfrags"]
    combinereps = argdict["combinereps"]
    ensemblgtf = argdict["ensemblgtf"]
    bamnofeature = argdict["bamnofeature"]
    trnatable = argdict["trnatable"]
    trnaaminofile = argdict["trnaaminofile"]
    trnaanticodonfile = argdict["trnaanticodonfile"]
    otherseqs = extraseqfile(argdict["otherseqs"])
    #print >>sys.stderr, argdict["otherseqs"]
    if "bamdir" not in argdict:
        bamdir = "./"
    bamdir = argdict["bamdir"]
    sampledata = samplefile(argdict["samplefile"], bamdir = bamdir)
    cores = argdict["cores"]
    threadmode = True
    if cores == 1:
        threadmode = False
    minpretrnaextend = 5
    mitochrom = None
    if argdict["mitochrom"]:
        mitochrom = argdict["mitochrom"]
    sizefactor = defaultdict(lambda: 1)
    if argdict["sizefactors"]:
        sizefactor = getsizefactors(argdict["sizefactors"]) 
        for currsample in sampledata.getsamples():
            if currsample not in sizefactor:
                print >>sys.stderr, "Size factor file "+argdict["sizefactors"]+" missing "+currsample
                sys.exit(1)
        
    bedfiles = list()
    
    if argdict["bedfile"]  is not None:
        bedfiles = argdict["bedfile"]
    trnalocifiles = list()
    if argdict["trnaloci"] is not None:
        trnalocifiles = argdict["trnaloci"]
    maturetrnas = list()
    if argdict["maturetrnas"] is not None:
        maturetrnas = argdict["bedfile"]
        
    #trnalocifiles = argdict["trnaloci"]
    #maturetrnas=argdict["maturetrnas"]
    
    genetypefile = argdict["genetypefile"]
    locifiles = argdict["trnaloci"]
    maturetrnafiles = argdict["maturetrnas"]     
    trnaaminofilename = argdict["trnaaminofile"]
    trnaanticodonfilename = argdict["trnaanticodonfile"]
    trnanormfile = argdict["trnanormfile"]
    allreadsnormfile = argdict["allreadsnormfile"]
    readlengthfile = argdict["readlengthfile"]    
    
    if argdict["realcountfile"] == "stdout":
        realcountfile = sys.stdout
    else:
        realcountfile = open(argdict["realcountfile"],"w")
    
    
    if argdict["countfile"] == "stdout":
        countfile = sys.stdout
    else:
        countfile = open(argdict["countfile"],"w")
    

    
    
    wholetrnas = dict()
    fivefrags = dict()
    threefrags = dict()
    trailerfrags = dict()
    otherfrags = dict()
    allfrags = dict()
    alltrnas = list()
    maturenames = dict()
    otherseqlist = dict()
    


    
    trnainfo = transcriptfile(trnatable)
    samplefiles = dict()
    
    samples = list(sampledata.getsamples())
    
        

    #print >>sys.stderr, "**READCOUNT**"
    try:
        featurelist = dict()
        trnaloci = dict()
        trnalist = dict()
        ensembllist = dict()
        otherseqlist = dict()
        for currfile in bedfiles:
            featurelist[currfile] = RangeBin(readfeatures(currfile))
        
        for currfile in locifiles:
            trnaloci[currfile] = RangeBin(readbed(currfile), binfactor = 10000)
        for currfile in maturetrnafiles:
            matlist = list(readbed(currfile))
            trnalist[currfile] = list(matlist)
            maturenames[currfile] = {curr.name:curr for curr in matlist}
        if ensemblgtf is not None:    
            embllist = RangeBin(readgtf(ensemblgtf, filtertypes = set()))
        else:
            embllist = None
        
        for currname, currfile in otherseqs.seqbed.iteritems():
            #print >>sys.stderr, currfile
            #print >>sys.stderr, len(list(readbed(currfile)))
            otherseqlist[currname] = getchromdict(readbed(currfile))
    except IOError as e:
        print >>sys.stderr, e
        sys.exit()
    
    #print >>sys.stderr, otherseqlist.keys()
    #print >>sys.stderr, "**"
    #sys.exit()
    featcount = defaultdict(int)
    bedlist = list(featurelist.iterkeys())
    maxmismatches = None
    allcounts = dict()
    poolmode = True
    starttime = time.time()
    #threadmode = False
    if threadmode:
        countqueue = Queue()
        threads = dict()
        if poolmode:
            countpool = Pool(processes=cores)
            arglist = list()
            for currsample in samples:
                currbam = sampledata.getbam(currsample)
                arglist.append(compressargs(currbam,currsample, trnainfo, trnaloci, trnalist,maturenames, otherseqlist = otherseqlist, embllist = embllist, featurelist = featurelist, maxmismatches = maxmismatches, bamnofeature = bamnofeature))
            #arglist = list((tuple([currsample, sampledata.getbam(currsample)]) for currsample in samples))
            results = countpool.map(counttypereadspool, arglist)
            for i, curr in enumerate(samples):
                allcounts[curr] = results[i]
        else:
            for currsample in samples:
                currbam = sampledata.getbam(currsample)
            
                threads[currsample] = Process(target=counttypereadsqueue,args = (countqueue,currsample,currbam, currsample,trnainfo, trnaloci, trnalist,maturenames), kwargs = { "embllist" : embllist, "featurelist" : featurelist, "maxmismatches" : maxmismatches, "bamnofeature" : bamnofeature})
                
                #threads[currsample] = Process(target=testqueue,args =  (countqueue,currsample,currbam, currsample,trnainfo, trnaloci, trnalist,maturenames), kwargs = { "embllist" : embllist, "featurelist" : featurelist, "maxmismatches" : maxmismatches, "bamnofeature" : bamnofeature})
                threads[currsample].start()
            for sample in threads.iterkeys():
                
                currsample, counts = countqueue.get()
                allcounts[currsample] = counts
                print >>sys.stderr, currsample+":" +str(time.time()-starttime)
            
            
            pass
    else:
        for currsample in samples:
            currbam = sampledata.getbam(currsample)
            allcounts[currsample] = counttypereads(currbam, currsample,trnainfo, trnaloci, trnalist,maturenames, otherseqlist = otherseqlist, embllist = embllist, featurelist = featurelist, maxmismatches = maxmismatches, bamnofeature = bamnofeature)
        
        
    emblbiotypes  = set(itertools.chain.from_iterable(curr.emblbiotypes for curr in allcounts.values()))        
    bedtypes  = set(itertools.chain.from_iterable(curr.bedtypes for curr in allcounts.values())) 
    extraseqtypes  = set(itertools.chain.from_iterable(curr.extraseqtypes for curr in allcounts.values()))   
    print >>sys.stderr, bedtypes
    printtypefile(countfile, samples, sampledata,allcounts,trnalist, trnaloci, bedtypes, emblbiotypes,sizefactor, countfrags = countfrags , extraseqtypes = extraseqtypes)
    printrealcounts(realcountfile, samples, sampledata,allcounts,trnalist, trnaloci, bedtypes, emblbiotypes , extraseqtypes = extraseqtypes)

    #printrealcounts()
    if readlengthfile is not None:
        printlengthfile(readlengthfile, samples, allcounts)

    if trnaaminofilename is not None:
        printaminocounts(trnaaminofilename, sampledata, allcounts, sizefactor)
    if trnaanticodonfile is not None:
        printanticodoncounts(trnaanticodonfilename, sampledata, allcounts, sizefactor)
def main(**argdict):
    argdict = defaultdict(lambda: None, argdict)
    countfrags = argdict["countfrags"]
    combinereps = argdict["combinereps"]
    ensemblgtf = argdict["ensemblgtf"]
    bamnofeature = argdict["bamnofeature"]
    trnatable = argdict["trnatable"]
    trnaaminofile = argdict["trnaaminofile"]
    sampledata = samplefile(argdict["samplefile"])
    minpretrnaextend = 5
    mitochrom = None
    if argdict["mitochrom"]:
        mitochrom = argdict["mitochrom"]
    sizefactor = defaultdict(lambda: 1)
    if argdict["sizefactors"]:
        sizefactor = getsizefactors(argdict["sizefactors"]) 
        for currsample in sampledata.getsamples():
            if currsample not in sizefactor:
                print >>sys.stderr, "Size factor file "+argdict["sizefactors"]+" missing "+currsample
                sys.exit(1)
        
    bedfiles = list()
    
    if argdict["bedfile"]  is not None:
        bedfiles = argdict["bedfile"]
    trnalocifiles = list()
    if argdict["trnaloci"] is not None:
        trnalocifiles = argdict["trnaloci"]
    maturetrnas = list()
    if argdict["maturetrnas"] is not None:
        maturetrnas = argdict["bedfile"]
        
    #trnalocifiles = argdict["trnaloci"]
    #maturetrnas=argdict["maturetrnas"]
    
    genetypefile = argdict["genetypefile"]
    locifiles = argdict["trnaloci"]
    maturetrnafiles = argdict["maturetrnas"]    
    trnaaminofilename = argdict["trnaaminofile"]
    trnanormfile = argdict["trnanormfile"]
    allreadsnormfile = argdict["allreadsnormfile"]
    readlengthfile = argdict["readlengthfile"]    
    
    if argdict["realcountfile"] == "stdout":
        realcountfile = sys.stdout
    else:
        realcountfile = open(argdict["realcountfile"],"w")
    
    
    if argdict["countfile"] == "stdout":
        countfile = sys.stdout
    else:
        countfile = open(argdict["countfile"],"w")
    

    
    
    wholetrnas = dict()
    fivefrags = dict()
    threefrags = dict()
    trailerfrags = dict()
    otherfrags = dict()
    allfrags = dict()
    alltrnas = list()
    maturenames = dict()
    
    
    ncrnaorder = defaultdict(int)
    #This decides the priority for ensembl genes, the rest I don't really care about yet
    for i, curr in enumerate(reversed(list(["snoRNA","miRNA", "rRNA","snRNA","misc_RNA","lincRNA", "protein_coding"]))):
        ncrnaorder[curr] = i + 1
    fullpretrnathreshold = 2
    
    trnainfo = transcriptfile(trnatable)
    samplefiles = dict()
    
    samples = list(sampledata.getsamples())
    
        

    #print >>sys.stderr, "**READCOUNT**"
    try:
        featurelist = dict()
        trnaloci = dict()
        trnalist = dict()
        ensembllist = dict()
        for currfile in bedfiles:
            featurelist[currfile] = RangeBin(readfeatures(currfile))
        
        for currfile in locifiles:
            trnaloci[currfile] = RangeBin(readbed(currfile), binfactor = 10000)
        for currfile in maturetrnafiles:
            matlist = list(readbed(currfile))
            trnalist[currfile] = list(matlist)
            maturenames[currfile] = {curr.name:curr for curr in matlist}
        if ensemblgtf is not None:    
            embllist = RangeBin(readgtf(ensemblgtf, filtertypes = set()))
        else:
            embllist = None
    
    except IOError as e:
        print >>sys.stderr, e
        sys.exit()
    
    
    #featurelist = list(readbed(sys.argv[1]))
    #featurelist = list(readbed(sys.argv[1]))
    
    #./countcomplete.py hg19-nontrnas.bed hg19-tRNAs.bed hg19-complete-tRNApad.fa
    '''
    ./countcomplete.py ../combinedb/sacCer3-fatRNAs.bed sacCer3-agingtranscripts.bed >sacCer3-agingcount.txt
    '''
    
    featcount = defaultdict(int)
    
    
    
    
    #featurelist = list(curr for curr in featurelist if curr.name == 'unknown20')
    #print >>sys.stderr, "***"
    
    #lengths = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    counts = defaultdict(lambda: defaultdict(int))
    emblcounts = defaultdict(lambda: defaultdict(int))
    trnacounts = defaultdict(lambda: defaultdict(int))
    trnawholecounts = defaultdict(lambda: defaultdict(int))
    trnafivecounts = defaultdict(lambda: defaultdict(int))
    trnathreecounts = defaultdict(lambda: defaultdict(int))
    trnatrailercounts  = defaultdict(lambda: defaultdict(int))
    trnaantisense  = defaultdict(lambda: defaultdict(int))
    
    trnalocuscounts = defaultdict(lambda: defaultdict(int))
    trnatrailercounts = defaultdict(lambda: defaultdict(int))
    
    fulltrnalocuscounts  = defaultdict(lambda: defaultdict(int))
    partialtrnalocuscounts  = defaultdict(lambda: defaultdict(int))
    trnalocustrailercounts = defaultdict(lambda: defaultdict(int))
    trnaaminocounts = defaultdict(lambda: defaultdict(int))
    
    readlengths = defaultdict(lambda: defaultdict(int))
    trnareadlengths = defaultdict(lambda: defaultdict(int))
    pretrnareadlengths = defaultdict(lambda: defaultdict(int))
    
    aminos = set()
    
    othercounts = defaultdict(int)
    
    maxoffset = 10
    
    bedlist = list(featurelist.iterkeys())
    locilist = list(trnaloci.iterkeys())
    #print >>sys.stderr, bedlist
    #sys.exit()
    '''
    I use readbins rather than searching the bamfiles for regions here because I want to ensure that each read is only mapped to a single feature
    
    '''
    
    trnasamplecounts  = defaultdict(int)
    totalsamplecounts = defaultdict(int)
    emblbiotypes = set()
    
    indelaligns = defaultdict(int)
    
    for currsample in samples:
        currbam = sampledata.getbam(currsample)
        
        try:
            #print >>sys.stderr, currbam
            
            
            if not os.path.isfile(currbam+".bai") or  os.path.getmtime(currbam+".bai") < os.path.getmtime(currbam):
                pysam.index(""+currbam)
            bamfile = pysam.Samfile(""+currbam, "rb" )
            if bamnofeature:
                outname = os.path.splitext(currbam)[0]+"_nofeat.bam"
                outbamnofeature =  pysam.Samfile( outname, "wb", template =  bamfile)
        except IOError as ( strerror):
            print >>sys.stderr, strerror
            sys.exit()
        #continue #point0
        for i, currread in enumerate(getbamrange(bamfile, primaryonly = True)):

            isindel = False
    
            gotread = False
            #continue #point1
            totalsamplecounts[currsample] += 1
            if len(currread.data["CIGAR"]) > 1:
                
                #print >>sys.stderr, currread.bedstring()
                #print >>sys.stderr, currread.data["CIGAR"]
                #deletions
                if sum(curr[1] for curr in currread.data["CIGAR"] if curr[0] == 2):
                    indelaligns[currsample] += 1
                    pass
                #insertions
                if sum(curr[1] for curr in currread.data["CIGAR"] if curr[0] == 2):
                    #indelaligns[currsample] += 1
                    pass
                isindel = True
                #continue

            else:
                pass
                #continue
            readlength = len(currread.data['seq'])
            readlengths[currsample][readlength] += 1
            #continue #point2
            for currbed in locilist:
                for currfeat in trnaloci[currbed].getbin(currread):
                    expandfeat = currfeat.addmargin(30)
                    if currfeat.coverage(currread) > 10 and (currread.start + minpretrnaextend <= currfeat.start or currread.end - minpretrnaextend >= currfeat.end):
                        pretrnareadlengths[currsample][readlength] += 1
                        trnalocuscounts[currsample][currbed] += 1
                        if currread.start + fullpretrnathreshold <  currfeat.start and currread.end - fullpretrnathreshold + 3 >  currfeat.end:
                            fulltrnalocuscounts[currsample][currbed] += 1
                        else:# currread.start + fullpretrnathreshold <  currfeat.start or currread.end - fullpretrnathreshold +3 >  currfeat.end:
                            partialtrnalocuscounts[currsample][currbed] += 1
                            #print >>sys.stderr, "***"
                        gotread = True
                        break
                    if currfeat.getdownstream(30).coverage(currread) > 10:
                        pretrnareadlengths[currsample][readlength] += 1
                        trnalocustrailercounts[currsample][currbed] += 1 
                        #print >>sys.stderr, currfeat.bedstring()
                        gotread = True
                        break
                    elif expandfeat.antisense().coverage(currread) > 5:
                        trnaantisense[currsample][currbed] += 1
                        gotread = True
                        break
            if gotread: 
                continue
            #continue #point3
                
            for currbed in trnalist:
                if currread.chrom in maturenames[currbed]:
                    currfeat = maturenames[currbed][currread.chrom]
                    if currread.strand == "+":
                        trnareadlengths[currsample][readlength] += 1
                        trnasamplecounts[currsample] += 1
                        trnacounts[currsample][currbed] += 1
                        
                        trnaaminocounts[currsample][trnainfo.getamino(currfeat.name)] += 1
                        aminos.add(trnainfo.getamino(currfeat.name))
                        if countfrags:
                            fragtype = getfragtype(currfeat, currread)
                            if fragtype == "Whole":
                                
                                trnawholecounts[currsample][currbed] += 1
                            elif fragtype == "Fiveprime":
                                trnafivecounts[currsample][currbed] += 1
                            elif fragtype == "Threeprime":
                                trnathreecounts[currsample][currbed] += 1
                            elif fragtype == "Trailer":
                                trnatrailercounts[currsample][currbed] += 1
                        gotread = True
                        break
                            #print >>sys.stderr, str(currread.start - currfeat.start)+"-"+str(currread.end - currfeat.start)  
                            #print >>sys.stderr, str(currfeat.start - currfeat.start)+"-"+str(currfeat.end - currfeat.start)
                            #print >>sys.stderr, "****"
                    elif currfeat.antisense().coverage(currread) > 10:
                        trnaantisense[currsample][currbed] += 1
                        gotread = True
                        break
            if gotread: 
                continue
            #continue #point4
            if embllist is not None:
                currtype = None
                for currfeat in embllist.getbin(currread):
                    if currfeat.coverage(currread) > 10: 
                        if currfeat.data["biotype"] == "processed_transcript":
                            #print >>sys.stderr, currfeat.bedstring()
                            
                            pass
    
                        if currtype is None or ncrnaorder[currfeat.data["biotype"]] > ncrnaorder[currtype]:
                            currtype= currfeat.data["biotype"]
                            if mitochrom == currread.chrom:
                                currtype = "mito"+currtype
                        
                        
                        
                if currtype is not None:
                    emblcounts[currsample][currtype] += 1
                    emblbiotypes.add(currtype)
                    gotread = True
                        #print >>sys.stderr, currbam +":"+ currbed
            if gotread: 
                continue
            #continue #point5

            for currbed in bedlist:
                
                for currfeat in featurelist[currbed].getbin(currread):
                    if currfeat.coverage(currread) > 10:
                        counts[currsample][currbed] += 1
                        gotread = True
                        break
                        #print >>sys.stderr, currbam +":"+ currbed
    
            othercounts[currsample] += 1
            if not gotread and embllist is not None and mitochrom == currread.chrom:
                currtype = "Mitochondrial_other"
                emblcounts[currsample][currtype] += 1
                emblbiotypes.add(currtype)
            if not gotread and bamnofeature:
                outbamnofeature.write(currread.data["bamline"])
    #currently not counting indels, but might later.
    for currsample in samples:
        pass
        #print >>sys.stderr, currsample+" indels: "+ str(indelaligns[currsample])+"/"+str(totalsamplecounts[currsample])+":"+str(1.*indelaligns[currsample]/totalsamplecounts[currsample])
    def sumsamples(countdict,sampledata, repname, currfeat = None, sizefactors = defaultdict(lambda: 1)):
        if currfeat is None: #To account for the "other" counts, which don't have a feature
            return sum(countdict[currsample]/sizefactors[currsample] for currsample in sampledata.getrepsamples(repname))
        else:
            sum(sizefactors[currsample] for currsample in sampledata.getrepsamples(repname))
            return sum(countdict[currsample][currfeat]/sizefactors[currsample] for currsample in sampledata.getrepsamples(repname))
    #tRNA-Ser-AGA-1-1    
        
    if combinereps:
        replicates = list(sampledata.allreplicates())
        print >>countfile, "\t".join(replicates)
        for currbed in trnalist:
            
            if countfrags:
                #sumsamples(trnafivecounts,sampledata,currrep, currfeat)
                
                print  >>countfile, "tRNA_wholecounts\t"+"\t".join(str(sumsamples(trnawholecounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
                print  >>countfile, "tRNA_fiveprime\t"+"\t".join(str(sumsamples(trnafivecounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
                print  >>countfile, "tRNA_threeprime\t"+"\t".join(str(sumsamples(trnathreecounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
                print  >>countfile, "tRNA_other\t"+"\t".join(str(sumsamples(trnacounts,sampledata,currrep, currbed, sizefactors = sizefactor) - (sumsamples(trnafivecounts,sampledata,currrep, currbed, sizefactors = sizefactor) + sumsamples(trnathreecounts,sampledata,currrep, currbed, sizefactors = sizefactor) + sumsamples(trnawholecounts,sampledata,currrep, currbed, sizefactors = sizefactor))) for currrep in replicates)
                print  >>countfile, "tRNA_antisense\t"+"\t".join(str(sumsamples(trnaantisense,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
            else:
                
                print  >>countfile, "tRNA\t"+"\t".join(str(sumsamples(trnacounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)        
            
        
        for currbed in locilist:
            #print >>sys.stderr, currbed 
            if countfrags:
                print  >>countfile, "pretRNA_full\t"+"\t".join(str(sumsamples(fulltrnalocuscounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
                print  >>countfile, "pretRNA_partial\t"+"\t".join(str(sumsamples(partialtrnalocuscounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
                print  >>countfile, "pretRNA_trailer\t"+"\t".join(str(sumsamples(trnalocustrailercounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
            else:
                print  >>countfile, "pretRNA\t"+"\t".join(str(sumsamples(trnalocuscounts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
        for currbiotype in emblbiotypes:
            print  >>countfile, currbiotype+"\t"+"\t".join(str(sumsamples(emblcounts,sampledata,currrep, currbiotype, sizefactors = sizefactor)) for currrep in replicates)
        for currbed in bedlist:  
            print  >>countfile, os.path.basename(currbed)+"\t"+"\t".join(str(sumsamples(counts,sampledata,currrep, currbed, sizefactors = sizefactor)) for currrep in replicates)
        print  >>countfile, "other"+"\t"+"\t".join(str(sumsamples(othercounts,sampledata,currrep, sizefactors = sizefactor)) for currrep in replicates)
    else:
        print  >>countfile, "\t".join(samples)
        
        for currbed in trnalist:
            
            if countfrags:
                print  >>countfile, "tRNA_wholecounts\t"+"\t".join(str(trnawholecounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
                print  >>countfile, "tRNA_fiveprime\t"+"\t".join(str(trnafivecounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
                print  >>countfile, "tRNA_threeprime\t"+"\t".join(str(trnathreecounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
                print  >>countfile, "tRNA_other\t"+"\t".join(str((trnacounts[currsample][currbed] - (trnathreecounts[currsample][currbed] + trnafivecounts[currsample][currbed] + trnawholecounts[currsample][currbed]))/sizefactor[currsample]) for currsample in samples)
                print  >>countfile, "tRNA_antisense\t"+"\t".join(str(trnaantisense[currsample][currbed]/sizefactor[currsample]) for currsample in samples)

                
                
            else:
                print  >>countfile, "tRNA\t"+"\t".join(str(trnacounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
        for currbed in locilist:
            if countfrags:
                print  >>countfile, "pretRNA_full\t"+"\t".join(str(fulltrnalocuscounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
                print  >>countfile, "pretRNA_partial\t"+"\t".join(str(partialtrnalocuscounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
                print  >>countfile, "pretRNA_trailer\t"+"\t".join(str(trnalocustrailercounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
            else:
                print  >>countfile, "pretRNA\t"+"\t".join(str(trnalocuscounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
        for currbiotype in emblbiotypes:
            print  >>countfile, currbiotype+"\t"+"\t".join(str(emblcounts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
        for currbed in bedlist:
            print  >>countfile, os.path.basename(currbed)+"\t"+"\t".join(str(counts[currsample][currbed]/sizefactor[currsample]) for currsample in samples)
        print  >>countfile, "other"+"\t"+"\t".join(str(othercounts[currsample]/sizefactor[currsample]) for currsample in samples)
        
    if realcountfile:      
        print  >>realcountfile, "\t".join(samples)
        
        for currbed in trnalist:
            
            if countfrags:
                print  >>realcountfile, "tRNA_wholecounts\t"+"\t".join(str(trnawholecounts[currsample][currbed]) for currsample in samples)
                print  >>realcountfile, "tRNA_fiveprime\t"+"\t".join(str(trnafivecounts[currsample][currbed]) for currsample in samples)
                print  >>realcountfile, "tRNA_threeprime\t"+"\t".join(str(trnathreecounts[currsample][currbed]) for currsample in samples)
                print  >>realcountfile, "tRNA_other\t"+"\t".join(str((trnacounts[currsample][currbed] - (trnathreecounts[currsample][currbed] + trnafivecounts[currsample][currbed] + trnawholecounts[currsample][currbed]))) for currsample in samples)
                print  >>realcountfile, "tRNA_antisense\t"+"\t".join(str(trnaantisense[currsample][currbed]) for currsample in samples)

                
                
            else:
                print  >>realcountfile, "tRNA\t"+"\t".join(str(trnacounts[currsample][currbed]) for currsample in samples)

        for currbed in locilist:
            if countfrags:
                print  >>realcountfile, "pretRNA_full\t"+"\t".join(str(fulltrnalocuscounts[currsample][currbed]) for currsample in samples)
                print  >>realcountfile, "pretRNA_partial\t"+"\t".join(str(partialtrnalocuscounts[currsample][currbed]) for currsample in samples)
                print  >>realcountfile, "pretRNA_trailer\t"+"\t".join(str(trnalocustrailercounts[currsample][currbed]) for currsample in samples)
            else:
                print  >>realcountfile, "pretRNA\t"+"\t".join(str(trnalocuscounts[currsample][currbed]) for currsample in samples)
        for currbiotype in emblbiotypes:
            #print >>sys.stderr, emblcounts[currsample]
            print  >>realcountfile, currbiotype+"\t"+"\t".join(str(emblcounts[currsample][currbiotype]) for currsample in samples)
        for currbed in bedlist:
            print  >>realcountfile, os.path.basename(currbed)+"\t"+"\t".join(str(counts[currsample][currbed]) for currsample in samples)
        print  >>realcountfile, "other"+"\t"+"\t".join(str(othercounts[currsample]) for currsample in samples)
        
        
    realcountfile.close()
    
    #print >>sys.stderr, trnaaminocounts
    if trnaaminofilename is not None:
        trnaaminofile = open(trnaaminofilename, "w")
        if combinereps:
            print  >>trnaaminofile, "\t".join(replicates)
            for curramino in aminos:
                #print >>sys.stderr, curramino
                print >>trnaaminofile, curramino+"\t"+"\t".join(str(sumsamples(trnaaminocounts,sampledata,currrep, curramino, sizefactors = sizefactor)) for currrep in replicates)
        else:
            print  >>trnaaminofile, "\t".join(samples)
            for curramino in aminos:
                print >>trnaaminofile, curramino+"\t"+"\t".join(str(trnaaminocounts[currsample][curramino]/sizefactor[currsample]) for currsample in samples)
            
            
    if trnanormfile is not None:
        #samples trnasamplecounts.keys()
        trnanormfile = open(trnanormfile, "w")
        mean = 1.*sum(trnasamplecounts.values())/len(trnasamplecounts.values())
        print >>trnanormfile, "\t".join(samples)
        print >>trnanormfile, "\t".join(str(trnasamplecounts[currsample]/mean) for currsample in samples)
        
    if allreadsnormfile is not None:    
        allreadsnormfile = open(allreadsnormfile, "w")
        mean = 1.*sum(totalsamplecounts.values())/len(totalsamplecounts.values())
        print >>allreadsnormfile,"\t".join(samples)
        print >>allreadsnormfile,"\t".join(str(totalsamplecounts[currsample]/mean) for currsample in samples)
    if readlengthfile is not None:    
        readlengthfile = open(readlengthfile, "w")
        print >>readlengthfile, "Length\tSample\tother\ttrnas\tpretrnas"
        for currsample in readlengths.keys():
            for curr in range(0,max(readlengths[currsample].keys())+1):
                othercount = trnareadlengths[currsample][curr] + pretrnareadlengths[currsample][curr]
                print >>readlengthfile, str(curr)+"\t"+currsample+"\t"+str(readlengths[currsample][curr] - othercount)+"\t"+str(trnareadlengths[currsample][curr]) +"\t"+str(pretrnareadlengths[currsample][curr])
        
        
        
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
    parser.add_argument('--samplefile',
                       help='Sample file in format')
    parser.add_argument('--sizefactors',
                       help='Optional file including size factors that will be used for normalization')
    parser.add_argument('--bedfile',  nargs='*', default=list(),
                       help='bed file with non-tRNA features')
    parser.add_argument('--gtffile',  nargs='+', default=list(),
                       help='gtf file with non-tRNA features')
    parser.add_argument('--ensemblgtf',
                       help='ensembl gtf file with tRNA features')
    parser.add_argument('--gtftrnas',  nargs='+', default=list(),
                       help='gtf file with tRNA features')
    parser.add_argument('--trnaloci',  nargs='+', default=list(),
                       help='bed file with tRNA features')
    parser.add_argument('--maturetrnas',  nargs='+', default=list(),
                       help='bed file with mature tRNA features')
    parser.add_argument('--countfrags', action="store_true", default=False,
                       help='Seperate tRNA fragment types')
    parser.add_argument('--combinereps', action="store_true", default=False,
                       help='Sum samples that are replicates')
    parser.add_argument('--trnatable',
                       help='table of tRNA features')
    parser.add_argument('--trnaaminofile',
                       help='table of tRNAs by amino acid')
    parser.add_argument('--trnanormfile',
                       help='Create normalization file to use to normalize to total tRNA reads')
    parser.add_argument('--allreadsnormfile',
                       help='Create normalization file to use to normalize to total reads')
    parser.add_argument('--readlengthfile',
                       help='optional read lengths table')
    
    
    parser.add_argument('--bamnofeature', action="store_true", default=False,
                       help='Create bam file output for reads without a feature')
    parser.add_argument('--mitochrom',
                       help='Optional name of mitochondrial chromosome in database (Used to specially label mitchondrial features)')
    
    args = parser.parse_args()
    argvars = vars(args)
    main(**argvars)
    #main(samplefile = args.samplefile, sizefactors = args.sizefactors, bedfile = args.bedfile, gtffile=args.gtffile,ensemblgtf=args.ensemblgtf,gtfrnas=args.gtfrnas,trnaloci=args.trnaloci)        
