#!/usr/bin/env python3

import re
import random
import math
import os
import tempfile
import sys
import getopt
import itertools
import codecs
import collections
import time
import string
import gzip
import os.path
from collections import defaultdict
from trnasequtils import *




class tRNAlocus:
    def __init__(self, loc, seq, score, amino, anticodon, intronseq, intron, rawseq = None):
        self.name = loc.name
        self.loc = loc
        self.seq = seq
        self.score = score
        self.amino = amino
        self.anticodon = anticodon
        self.intronseq = intronseq
        self.intron = intron
        self.rawseq = rawseq

class tRNAtranscript:
    def __init__(self, seq, score, amino, anticodon, loci, intronseq, name = None, rawseq = None, artificialtrna = False):
        self.seq = seq
        self.score = score
        self.amino = amino
        self.anticodon = anticodon
        self.loci = tuple(loci)
        self.intronseqs = intronseq
        self.name = name
        self.artificialtrna = artificialtrna
        self.rawseq = rawseq
    def getmatureseq(self, addcca = True):
        prefix = ""
        #print >>sys.stderr, self.amino
        if self.amino == "His":
            prefix = "G"
        end = ""
        if addcca and not self.artificialtrna:
            end = "CCA"
        return prefix + self.seq + end 
    
            

def getuniquetRNAs(trnalist):
    sequencedict = defaultdict(list)
    scoredict = defaultdict(set)
    for curr in trnalist:
        sequencedict[curr.seq].append(curr)
    for currtrans in sequencedict.keys():
        scores = set(curr.score for curr in sequencedict[currtrans])
        anticodon = set(curr.anticodon for curr in sequencedict[currtrans])
        amino = set(curr.amino for curr in sequencedict[currtrans])
        introns = set(curr.intronseq for curr in sequencedict[currtrans])
        #remove psuedo if there's a real one somewheres
        if len(anticodon) > 1:
            anticodon.discard('Xxx')
            amino.discard('X')
        if introns == set([""]):
            introns = set()
        loci = list(curr for curr in sequencedict[currtrans])
        if len(scores) > 1:
            #print >>sys.stderr, "Multiple scores"
            pass
        if len(anticodon) > 1:
            print("tRNA file contains identical tRNAs with seperate anticodons, cannot continue", file=sys.stderr)
            sys.exit()
        yield tRNAtranscript(currtrans, scores,list(amino)[0],list(anticodon)[0],loci, introns)
        
        

def readrnacentral(scanfile,chromnames, mode = 'locus'):
    reftochrom = dict()
    convertfile = open(chromnames)
    
    for line in convertfile:
        if line.startswith("#"):
            continue
        fields = line.split()
        reftochrom[fields[1]] = fields[0]
    orgname = "genome"
   
    #print reftochrom
    trnafile = open(scanfile)
    transcriptinfo = defaultdict(list)
    for line in trnafile:
        fields = line.split(",")
        #print len(fields)
        #print "**"
        if fields[0] == 'Entry number':
            continue
        #print fields[1]
        if True:
            name = fields[2]
            amino = fields[11]
            anticodon = fields[12]
            sequence = fields[21]
            score = float(fields[18])
            #print >>sys.stderr, score
            #print name+":"+amino+":"+anticodon
            #print sequence
            gbchrom = re.sub(r'\.\d+$', '', fields[7])
            if gbchrom not in reftochrom:
                continue
            chrom = reftochrom[gbchrom]
            #print >>sys.stderr, chrom
            start = int(fields[8].split('-')[0]) - 1
            end = fields[8].split('-')[1]
            transcriptname = re.sub(r'-\d+$', '', fields[2])
            if fields[9] == "yes":
                strand = "-"
                
            elif fields[9] == "no":
                strand = "+"
            currtRNA = GenomeRange(orgname, chrom,start,end, name = name,strand = strand,orderstrand = True)
            #print >>sys.stderr, chrom
            intronnums = set()
            intronseqs = ""
            for currintron in fields[14:16]:
                intronmatch  = re.search(r'(\d+)\.\.(\d+)',currintron)
                if intronmatch:
                    intronstart = int(intronmatch.group(1))-1
                    intronend = int(intronmatch.group(2))
                    intronnums |= set(range(intronstart, intronend))
                    intronseqs = sequence[intronstart:intronend]
            
            newseq = ''

            
            for i in range(len(sequence)):
                if i in set(intronnums):
                    newseq += '-'
                else:
                    newseq += sequence[i]
                    
            rawseq = newseq.replace('-','')
            currlocus = tRNAlocus(currtRNA,rawseq, score,amino,anticodon,intronseqs, None)
            transcriptinfo[transcriptname].append(currlocus)
            if mode == 'locus':
                yield currlocus
        allseqs = dict()
    print(len(list(transcriptinfo.keys())), file=sys.stderr)
    for currtrans in transcriptinfo.keys():
        if len(set(curr.seq for curr in transcriptinfo[currtrans])) > 1:
            print("multiple", file=sys.stderr)
        #print transcriptinfo[currtrans][0].seq
        if transcriptinfo[currtrans][0].seq in allseqs:
            print("duplicate:" + currtrans + ":"+allseqs[transcriptinfo[currtrans][0].seq], file=sys.stderr)
        allseqs[transcriptinfo[currtrans][0].seq] = currtrans
        if mode == 'transcript':
            #print >>sys.stderr,currtrans 
            yield tRNAtranscript( transcriptinfo[currtrans][0].seq, set(curr.score for curr in transcriptinfo[currtrans]), transcriptinfo[currtrans][0].amino, transcriptinfo[currtrans][0].anticodon, set(transcriptinfo[currtrans]), transcriptinfo[currtrans][0].intronseq, name = currtrans)
        

def readtRNAscan(scanfile, genomefile, mode = None):
    #mode = 'gtRNAdb'
    trnalist = list()
    orgname = "genome"
    if  hasattr(scanfile ,'read'):
        trnascan = scanfile
    else:
        trnascan = open(scanfile, "r")
    trnascore = dict()
    trnaanticodon = dict()
    trnaamino = dict()
    tRNAintron = dict()
    trnas = dict()
    
    for currline in trnascan:
        
        if  currline.startswith("Sequence") or  currline.startswith("Name") or  currline.startswith("------"):
            continue
        fields = currline.split()
        
        if mode == "gtRNAdb":
            #print >>sys.stderr, fields[6:8]
            del fields[6:8]
        curramino = fields[4]
        currac = fields[5]
        
        if currac == "???":
            currac = 'Xxx'
        if fields[2] > fields[3]:
            end = int(fields[3]) - 1
            start = int(fields[2])
            
        else:
            end = int(fields[3])
            start = int(fields[2]) - 1
        currchrom = fields[0]
        trnanum = fields[1]
        currtRNA = GenomeRange(orgname, currchrom,start,end, name = currchrom+"."+"tRNA"+trnanum+"-"+curramino+currac,strand = "+",orderstrand = True)
        currtrans = currtRNA

        trnaamino[currtrans.name] = curramino
        trnaanticodon[currtrans.name] = currac
        #print >>sys.stderr, "**".join(fields)#currline
        trnascore[currtrans.name] =  float(fields[8])
        trnas[currtrans.name] =  currtrans
    
    
        
        currtRNA.fastafile = genomefile

        trnalist.append(currtRNA)
        #print >>sys.stderr,trnalist
        if int(fields[6]) != 0:
            if currtRNA.strand ==  "-":
                intronstart = int(fields[2]) - int(fields[6]) 
                intronend = int(fields[2]) - int(fields[7]) +1
            else:
                intronstart = int(fields[6]) - int(fields[2]) - 1
                intronend = int(fields[7]) - int(fields[2])
            tRNAintron[currtRNA.name] = tuple([intronstart, intronend])
    #should add check to make sure all trnas are grabbed
    trnaseqs = getseqdict(trnalist, faifiles = {orgname:genomefile+".fai"})
    intronseq = defaultdict(str)
    #print >>sys.stderr,trnalist
    for curr in trnaseqs.keys():
        currintron = None
        if curr in tRNAintron:
            start = tRNAintron[curr][0]
            end = tRNAintron[curr][1]
            intronseq[curr] = trnaseqs[curr][start:end]
            trnaseqs[curr] = trnaseqs[curr][:start] + trnaseqs[curr][end:]
            currintron = tRNAintron[curr]
            
        yield tRNAlocus(trnas[curr], trnaseqs[curr], trnascore[curr],trnaamino[curr],trnaanticodon[curr],intronseq[curr], currintron)
def striplocus(trnaname):
    return re.sub(r"\-\d+$", "",trnaname)
#Sequence		tRNA     	Bounds   	tRNA	Anti	Intron Bounds	Inf	HMM	2'Str	Hit	      	Isotype	Isotype	Type
#Name    	tRNA #	Begin    	End      	Type	Codon	Begin	End	Score	Score	Score	Origin	Pseudo	CM	Score	
#--------
def readtRNAdb(scanfile, genomefile, trnamap):

    
    #mode = 'gtRNAdb'
    trnalist = list()
    orgname = "genome"
    if  hasattr(scanfile ,'read'):
        trnascan = scanfile
    else:
        trnascan = open(scanfile, "r")
    trnascore = dict()
    trnaanticodon = dict()
    trnaamino = dict()
    tRNAintron = dict()
    trnas = dict()
    for linenum, currline in enumerate(trnascan):
        if  currline.startswith("Sequence") or  currline.startswith("Name") or  currline.startswith("------"):
            continue
        if len(currline) < 5:
            print("cannot read line: " +str(linenum) +" of "+scanfile, file=sys.stderr)
            continue
        fields = currline.split()
        #print (fields[0], file = sys.stderr)
        curramino = fields[4]
        currac = fields[5]
        
        if currac == "???":
            pass
            #currac = 'Xxx'
        if fields[2] > fields[3]:
            end = int(fields[3]) - 1
            start = int(fields[2])
            
        else:
            end = int(fields[3])
            start = int(fields[2]) - 1
        currchrom = fields[0]
        trnanum = fields[1]
        
        trnascanname = currchrom+"."+"trna"+trnanum+"-"+curramino+currac
        #print >>sys.stderr, trnamap.keys()
        #print >>sys.stderr, trnascanname
        shorttrnascanname = currchrom+"."+"trna"+trnanum
        if trnascanname in trnamap:
            currtRNA = GenomeRange(orgname, currchrom,start,end, name = trnamap[trnascanname],strand = "+",orderstrand = True)
        elif shorttrnascanname in trnamap:
            currtRNA = GenomeRange(orgname, currchrom,start,end, name = trnamap[shorttrnascanname],strand = "+",orderstrand = True)
        else:
            print("Skipping "+trnascanname+", has no transcript name", file=sys.stderr)
            continue
        currtrans = currtRNA

        trnaamino[currtrans.name] = curramino
        trnaanticodon[currtrans.name] = currac
        #print >>sys.stderr, "**".join(fields)#currline
        trnascore[currtrans.name] =  float(fields[8])
        trnas[currtrans.name] =  currtrans
    
        #print (currtrans.name, file = sys.stderr)
        currtRNA.fastafile = genomefile
        #print >>sys.stderr, genomefile
        trnalist.append(currtRNA)
        if int(fields[6]) != 0:
            if currtRNA.strand ==  "-":
                intronstart = int(fields[2]) - int(fields[6])
                intronend = int(fields[2]) - int(fields[7])  + 1
            else:
                intronstart = int(fields[6]) - int(fields[2]) 
                intronend = int(fields[7]) - int(fields[2]) + 1
            tRNAintron[currtRNA.name] = tuple([intronstart, intronend])
    trnaseqs = getseqdict(trnalist, faifiles = {orgname:genomefile+".fai"})
    intronseq = defaultdict(str)
    trnaloci = list()
    #print >>sys.stderr, "tRNAs: "+str(len(trnalist))
    #print >>sys.stderr, "tRNA seqs: "+str(len(trnaseqs.keys()))
    for curr in trnaseqs.keys():
        currintron = None
        if curr in tRNAintron:
            
            start = tRNAintron[curr][0]
            end = tRNAintron[curr][1]

            intronseq[curr] = trnaseqs[curr][start:end]
            trnaseqs[curr] = trnaseqs[curr][:start] + trnaseqs[curr][end:]
            currintron = tRNAintron[curr]

        trnaloci.append( tRNAlocus(trnas[curr], trnaseqs[curr], trnascore[curr],trnaamino[curr],trnaanticodon[curr],intronseq[curr],currintron))
    '''
    for curr in trnaseqs.iterkeys():
        allintronseqs[curr] = ''
        if curr in tRNAintron:
            start = tRNAintron[curr][0]
            end = tRNAintron[curr][1]
            intronseq[curr] = trnaseqs[curr][start:end]
            allintronseqs[curr] =  trnaseqs[curr][start:end]
            intronsplices[curr] = [start,end]
            trnaseqs[curr] = trnaseqs[curr][:start] + trnaseqs[curr][end:]
    '''        
            
        
    trnaloci.sort(key = lambda x: striplocus(x.name))
    for transname, currloci in itertools.groupby(trnaloci, lambda x: striplocus(x.name)):
        currlocus = list(currloci)
        yield tRNAtranscript( currlocus[0].seq, set(curr.score for curr in currlocus), currlocus[0].amino,currlocus[0].anticodon, set(currlocus), currlocus[0].intronseq, name = transname)

        
    

