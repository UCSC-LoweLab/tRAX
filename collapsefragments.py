#!/usr/bin/env python3

import re
import sys
import os.path
import itertools
import subprocess
from trnasequtils import *
from distutils.spawn import find_executable
import time


from collections import defaultdict

def readmultifastq(fqfile, fullname = False):
    #print chrom+":"+ chromstart+"-"+ chromend
    if fqfile == "stdin":
        fqfile = sys.stdin
    elif fqfile.endswith(".gz"):
        fqfile = gzip.open(fqfile, "rb")
    else:
        fqfile = open(fqfile, "r")
    currloc = 0
    currseq = None
    sequence = ""
    quality = ""
    if fullname:
        reheader = re.compile(r"\@(.+)")
    else:
        reheader = re.compile(r"\@([^\s\,]+)")
    qualheader = re.compile(r"\+([^\s\,]*)")
    readqual = False
    for line in fqfile:
        #print line
        line = line.rstrip("\n")
        seqheader = reheader.match(line)
        qheader = qualheader.match(line)
        if readqual:
            quality += line
            if len(quality) == len(sequence):
                yield currseq, sequence, quality
                readqual = False
        elif seqheader:
            currseq = seqheader.groups(1)[0]
            sequence = ""
            quality = ""
        elif qheader and readqual == False:
            readqual = True
            pass
        else:
            sequence += line
            
class prunedict:
    def __init__(self, maxkeys = 10000000):
        self.counts = defaultdict(int) 
        self.maxkeys = maxkeys 
        self.trimcutoff = max([10,self.maxkeys/100000])
        self.totalkeys = 0
        self.trimmed = 0
    def trim(self):
        #print >>sys.stderr, "**"
        newdict = defaultdict(int) 
        trimmed = 0
        for curr in self.counts.keys():
            if self.counts[curr] > self.trimcutoff:
                newdict[curr] = self.counts[curr]
            else:
                trimmed += 1
        self.trimmed += trimmed
        #print >>sys.stderr, str(trimmed)+"/"+str(self.totalkeys)+" at "+str(self.trimcutoff)
        self.totalkeys = len(list(self.counts.keys()))
        self.counts = newdict
        
    def __getitem__(self, key):
        if key not in self.counts:
            self.totalkeys += 1
        return self.counts[key]
    def __setitem__(self, key, count):
        if key not in self.counts:
            self.totalkeys += 1
        self.counts[key] = count
        if self.totalkeys > self.maxkeys:
            #self.trim()

            #print >>sys.stderr, str(len(self.counts.keys())) +"/"+ str(len(self.newdict.keys())) +":"+str(1.*len(self.counts.keys())/len(self.newdict.keys()))
            if self.totalkeys > self.maxkeys:
                self.trimcutoff *= 1.1
    def iterkeys(self):
        return list(self.counts.keys())
            

sampledata = samplefile(sys.argv[1])
samples = sampledata.getsamples()

seqcount = dict()
allmode = False

for currsample in samples: 
    maxseqs = 10000000 
    seqcount[currsample] = prunedict()            
    total = 0
    for name, seq, qual in readmultifastq(sampledata.getfastq(currsample)):
        seqcount[currsample][seq] += 1
        total += 1
        #if total % 100000 == 0:
            #print >>sys.stderr, str(total)
    if not allmode:  
        pass
        #seqcount[currsample].trim()

seqfile = open(sys.argv[2], "w")

allseqs = defaultdict(int)
totalreads = defaultdict(int)

for currsample in samples:
    for currseq in seqcount[currsample].keys():
        allseqs[currseq] += 1
        totalreads[currseq] += seqcount[currsample][currseq]

maxmissing = 0
print("\t".join(samples))    
for i, currseq in enumerate(allseqs.keys()):
    
    if not allmode and (allseqs[currseq] < 2 or totalreads[currseq] < 40):   #allseqs[currseq] < len(samples)
        
        currmax = max(seqcount[currsample][currseq] for currsample in samples) 
        if currmax > maxmissing:
            maxmissing = currmax
        continue
    seqname = "frag"+str(i+1)+"_"+str(len(currseq))
    print(seqname+"\t"+ "\t".join(str(seqcount[currsample][currseq]) for currsample in samples))

    print(">"+ seqname, file=seqfile)
    print(currseq, file=seqfile)
    
