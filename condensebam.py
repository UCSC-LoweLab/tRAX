#!/usr/bin/env python3

import pysam
import sys
import argparse
import os.path
from collections import defaultdict
from trnasequtils import *
import itertools



class prunedict:
    def __init__(self, maxkeys = 1000000, autotrim = True):
        self.counts = defaultdict(int) 
        self.maxkeys = maxkeys 
        self.trimcutoff = max([10,self.maxkeys/100000])
        self.totalkeys = 0
        self.trimmed = 0
        self.autotrim = autotrim
    def trim(self):
        #print >>sys.stderr, "**"
        newdict = defaultdict(int) 
        trimmed = 0
        #currcutoff = max([trimcutoff, minreads])
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
        if self.totalkeys > self.maxkeys and self.autotrim:
            self.trim()

            #print >>sys.stderr, str(len(self.counts.keys())) +"/"+ str(len(self.newdict.keys())) +":"+str(1.*len(self.counts.keys())/len(self.newdict.keys()))
            if self.totalkeys > self.maxkeys:
                self.trimcutoff *= 1.1                        
    def checktrim(self):
        return self.totalkeys > self.maxkeys
    def iterkeys(self):
        return iter(self.counts.keys())
    def toseqnums(self):
        numdict = dict()
        return {'seq'+str(i)+":"+str(self.counts[currseq]):currseq  for i, currseq in self.keys()}
    def resetmax(self, newmax):
        if len(list(self.counts.values())) > newmax:
            
            allcounts = sorted(self.counts.values())
            self.trimcutoff = allcounts[newmax]
            self.trim()
        

            
class nameseq:
    def __init__(self, maxkeys = 1000):
        self.names = dict()
        self.counts = defaultdict(lambda: prunedict(maxkeys = maxkeys, autotrim = False))
        self.maxkeys = maxkeys
    def addread(self,featname, seqname, seq): #maybe I only need to store the first seqname?
        self.counts[featname][seq] += 1
        if seq not in self.names:
            self.names[seq] = seqname  #only store one name for each sequence
        
        if self.counts[featname].checktrim():
            self.trim(featname) 
    def trim(self,featname, newmax = None):
        if newmax is not None:
            self.counts[featname].resetmax(newmax)
        self.counts[featname].trim()
        for curr in (set(self.names.keys()) - set(itertools.chain.from_iterable(iter(self.counts[currfeat].keys()) for currfeat in self.counts.keys()))):
            pass
            del self.names[curr]
        #print >>sys.stderr, str(len(self.counts.keys())) +"/"+ str(len(self.newdict.keys())) +":"+str(1.*len(self.counts.keys())/len(self.newdict.keys()))
                
        if self.counts[featname].totalkeys > self.maxkeys:
            self.counts[featname].trimcutoff *= 1.1
    def setmax(self, newmax):
        for featname in self.counts.keys():
            self.trim(featname, newmax  = newmax)

    def getseqnames(self):
        seqset = dict()
        seqnum = 0
        for currfeat in self.counts.keys():
            for currseq in self.counts[currfeat].keys():
                if currseq not in seqset:
                    seqset[self.names[currseq]] = "frag"+str(seqnum+1)+":"+str(self.counts[currfeat][currseq])
                    seqnum += 1
        return seqset
        #{list(self.names[seq])[0]:"frag"+str(i+1)+":"+str(self.counts[seq]) for i, curr in enumerate(self.counts.iterkeys())}
    def gettotal(self):
        count = 0
        for featname in self.counts.keys():
            count += len(list(self.counts[featname].counts.keys()))
        return count
    

           
   
            
def main(**argdict):

    trnatable = argdict["trnatable"]

    
    sampledata = samplefile(argdict["samplefile"])
    if "trnaloci"  in argdict:
        trnalocifiles = argdict["trnaloci"]
    maturetrnas = list()
    if "maturetrnas" in argdict:
        maturetrnas = argdict["maturetrnas"]
        
    
    samples = sampledata.getsamples()
    trnainfo = transcriptfile(trnatable)
    trnaseqcounts = nameseq()
    trnaloci = list()
    trnalist = list()
    
    try:
        for currfile in trnalocifiles:
            trnaloci.extend(list(readbed(currfile)))
        for currfile in maturetrnas:
            trnalist.extend(list(readbed(currfile)))

    except IOError as e:
        print(e, file=sys.stderr)
        sys.exit()
    nomultimap = False
    maxmismatches = None
    allowindels = False
    minpretrnaextend = 5

    for currsample in samples:
        
        currbam = sampledata.getbam(currsample)
        print(currsample, file=sys.stderr)
        #doing this thing here why I only index the bamfile if the if the index file isn't there or is older than the map file
        try:
            if not os.path.isfile(currbam+".bai") or os.path.getmtime(currbam+".bai") < os.path.getmtime(currbam):
                pysam.index(""+currbam)
            bamfile = pysam.Samfile(""+currbam, "rb" )  
        except IOError as xxx_todo_changeme:
            ( strerror) = xxx_todo_changeme
            print(strerror, file=sys.stderr)
            sys.exit(1)
        except pysam.utils.SamtoolsError:
            print("Can not index "+currbam, file=sys.stderr)
            print("Exiting...", file=sys.stderr)
            sys.exit(1)
                
        for currfeat in trnaloci:
            for currread in getbamrange(bamfile, currfeat.addmargin(30), singleonly = nomultimap, maxmismatches = maxmismatches,allowindels = allowindels):
                #gotta be more than 5 bases off one end to be a true pre-tRNA
                #might want to shove these to the real tRNA at some point, but they are for now just ignored


                if currfeat.coverage(currread) > 10 and (currread.start + minpretrnaextend <= currfeat.start or currread.end - minpretrnaextend >= currfeat.end):
                    trnaseqcounts.addread(currfeat.name, currread.name , currread.data["seq"])
                elif currfeat.getdownstream(30).coverage(currread) > 10:
                    trnaseqcounts.addread(currfeat.name, currread.name , currread.data["seq"])

        
        for currfeat in trnalist:
            for currread in getbamrange(bamfile, currfeat, singleonly = nomultimap, maxmismatches = maxmismatches,allowindels = allowindels):
                
                if not currfeat.strand == currread.strand:
                    continue
                if not currfeat.coverage(currread) > 10:
                    continue
                trnaseqcounts.addread(currfeat.name, currread.name , currread.data["seq"])
                pass
        break
    print("Got sequences", file=sys.stderr)
    trnaseqcounts.setmax(20)
    print("Total: "+str(trnaseqcounts.gettotal()), file=sys.stderr)
    seqnames = trnaseqcounts.getseqnames()
    headbamfile =   pysam.Samfile(""+sampledata.getbam(samples[0]), "rb" )
    headerout = False
    outfile = pysam.Samfile( "-", "wb", template = headbamfile )
    for currsample in samples:
        
        currbam = sampledata.getbam(currsample)
        #print >>sys.stderr, currsample
        #doing this thing here why I only index the bamfile if the if the index file isn't there or is older than the map file
        try:
            if not os.path.isfile(currbam+".bai") or os.path.getmtime(currbam+".bai") < os.path.getmtime(currbam):
                pysam.index(""+currbam)
            bamfile = pysam.Samfile(""+currbam, "rb" )
        except IOError as xxx_todo_changeme1:
            ( strerror) = xxx_todo_changeme1
            print(strerror, file=sys.stderr)
            sys.exit(1)
        except pysam.utils.SamtoolsError:
            print("Can not index "+currbam, file=sys.stderr)
            print("Exiting...", file=sys.stderr)
            sys.exit(1)
        #indexing takes longer than just stepping through the files            
        #nameindex = pysam.IndexedReads(bamfile)
        #nameindex.build()
        #nameindex.find()
        '''
        for i, currname in enumerate(seqnames.iterkeys()):
            if i % 1000 == 0:
                #print >>sys.stderr, i
                pass
            currseqname = seqnames[currname]
            
            #bamline = nameindex.find(currname)
            try:
                for bamline in nameindex.find(currname):
                    bamline.qname = currseqname
                    outfile.write(bamline)
            except KeyError:
                #print >>sys.stderr, currname + " not found"
                #print >>sys.stderr, e
                pass
        '''
        #print >>sys.stderr, currbam
        for i, currread in enumerate(getbamrange(bamfile)): #this part takes a long time because I need to step through the entire file
            #continue
            #print >>sys.stderr, currread.name
            #if currread.name in seqnames:

            if currread.data["seq"]  in trnaseqcounts.names:
                currseqname = trnaseqcounts.names[currread.data["seq"]]
                
                bamline = currread.data["bamline"]
                bamline.qname = currseqname
                outfile.write(bamline)
        break

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
    parser.add_argument('--samplefile',
                       help='Sample file in format')
    parser.add_argument('--trnaloci',  nargs='+', default=list(),
                       help='bed file with tRNA features')
    parser.add_argument('--maturetrnas',  nargs='+', default=list(),
                       help='bed file with mature tRNA features')
    parser.add_argument('--trnatable',
                       help='table of tRNA features')
    
    
    args = parser.parse_args()
        
    #main(samplefile=args.samplefile, bedfile=args.bedfile, gtffile=args.bedfile, ensemblgtf=args.ensemblgtf, trnaloci=args.trnaloci, onlyfullpretrnas=args.onlyfullpretrnas,removepseudo=args.removepseudo,genetypefile=args.genetypefile,trnacounts=args.trnacounts,maturetrnas=args.maturetrnas,nofrag=args.nofrag,nomultimap=args.nomultimap,maxmismatches=args.maxmismatches)
    argvars = vars(args)
    #argvars["countfile"] = "stdout"
    main(**argvars)