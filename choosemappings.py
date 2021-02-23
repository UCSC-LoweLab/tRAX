#!/usr/bin/env python

import re
import sys
import os.path
import itertools
import pysam
import subprocess
import argparse

from trnasequtils import *


from collections import defaultdict

'''
Here is where I need to use the tRNA ontology between mature tRNAs and chromosomes


'''
defminnontrnasize = 20

maxmaps = 50





def isprimarymapping(mapping):
    return not (mapping.flag & 0x0100 > 0)

def getbesttrnamappings(trnafile, bamout = True, logfile = sys.stderr, progname = None, fqname = None, libname = None,setcountfile = None, extraseqfilename = None, minnontrnasize = defminnontrnasize):
    
    trnadata = transcriptfile(trnafile)
    trnatranscripts = set(trnadata.gettranscripts())
    extraseqs = set()
    seqdata = None
    if extraseqfilename is not None:
        seqdata = extraseqfile(extraseqfilename)
        extraseqs = set(trnadata.gettranscripts())
    readtargets = set()

    hitscores = dict()
    curreadname = None
    readlength = None
    totalmatch = re.compile(r"^(?:(?P<startclip>\d+)S)?(?P<matchlength>\d+)M(?:(?P<endclip>\d+)S)?$")
    #19S17M
    totalclips = 0
    totalreads = 0
    multimaps = 0
    duperemove = 0
    shortened = 0
    mapsremoved = 0
    totalmaps = 0
    getgithash
    
    scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"
    gitversion, githash = getgithash(scriptdir)
    #prints the sam header that samtools need to convert to bam
    #grab and iterate through all mappings of a single read
    #most mappers will ensure that read mappings are in the same order as they were in the fastq file, with all mappings of the same read together
    #once the file has been sorted using "samtools sort", this doesn't work anymore.  I can't detect that here, so no error message will be output
    trnareads = 0
    maxreads = 0
    diffreads = 0
    
    uniquenontrnas = 0
    nonuniquenontrnas = 0
    #gzsam = gzip.open(samfile, "rb")
    bamfile = pysam.Samfile("-", "r" )
    sort = True
    sortjob = None
    ambanticodon = 0
    ambamino = 0
    ambtrna = 0
    acsets = defaultdict(int)
    aminosets = defaultdict(int)
    trnasetcounts = defaultdict(int)
    newheader = bamfile.header.to_dict()
    #print >>sys.stderr, newheader
    newheader["RG"] = list()
    newheader["RG"].append(dict())
    imperfect = 0
    extraimperfect = 0
    if  progname is  not None:
        newheader["PG"].append({"PN" :progname, "ID": progname,"VN":gitversion})
    if fqname is not None:
        newheader["RG"][0]["ID"] = fqname 
    if libname is not None:
        newheader["RG"][0]["LB"] = libname
    
    if bamout:
        outfile = pysam.Samfile( "-", "wb", header = newheader )
    else:
        outfile = pysam.Samfile( "-", "w", header = newheader )
    for pairedname, allmaps in itertools.groupby(bamfile,lambda x: x.qname):
        allmaps = list(allmaps)
        if sum(curr.flag & 0x004 > 0 for curr in allmaps):
            continue
        totalreads += 1
        #print >>sys.stderr, "**"+pairedname
        readlength = None
        hitscores = dict()
        readtargets = set()
        clipsize = 50
        mappings = 0
        currscore = None
        newset = set()
        readlength = None
        #print >>sys.stderr, "**||"
        
        #iterate through all mappings of the current read
        #print >>sys.stderr, "**"+str(len(list(allmaps)))
        for currmap in allmaps:
            
            tagdict = dict()
            for curr in currmap.tags:
                tagdict[curr[0]] = curr[1]
            totalmaps += 1
            if currmap.tid is -1:
                continue
            #print >>sys.stderr, bamfile.getrname(currmap.tid)
            chromname = bamfile.getrname(currmap.tid)
            #sys.exit()
    
            readlength = len(currmap.seq)
            mappings += 1
            #if this is the best mapping of this read, discard any previous mappings
            
            
            #print >>sys.stderr, tagdict["AS"]
            #if the current score is worse than the new one
            if currscore is None or currscore < tagdict["AS"]:
                #print >>sys.stderr, str(currscore) +"<"+str( tagdict["AS"])
                newset = set()
                newset.add(currmap)
                currscore = tagdict["AS"]

            #if this mappings is the same as the previous best mapping, add it to the list
            elif currscore == tagdict["AS"]:
                newset.add(currmap)
                
            else:
                pass
            #print >>sys.stderr, currscore
        #here is where I count a bunch of thing so I can report them at the end
        if mappings > 1:
            multimaps += 1
        if len(newset) < mappings:
            #print  >>sys.stderr, pairedname
            #print >>sys.stderr, str(len(newset))+"/"+str(mappings)
            
            shortened += 1
        #print str(len(newset))+"\t"+str(readlength)
        #best mappings are printed out here
        if len(newset) >= 50:
            maxreads += 1
            #print >>sys.stderr, len(newset)
        finalset = list()
        #
        #print >>sys.stderr, ",".join(bamfile.getrname(curr.tid)  for curr in newset)
        #print >>sys.stderr, trnatranscripts
        if sum(bamfile.getrname(curr.tid) in trnatranscripts for curr in newset) > 0:
            
            trnareads += 1
            diff = len(newset) - sum(bamfile.getrname(curr.tid) in trnatranscripts for curr in newset)
            anticodons = frozenset(trnadata.getanticodon(bamfile.getrname(curr.tid)) for curr in newset if bamfile.getrname(curr.tid) in trnatranscripts)
            aminos = frozenset(trnadata.getamino(bamfile.getrname(curr.tid)) for curr in newset if bamfile.getrname(curr.tid) in trnatranscripts)
            trnamappings = list(curr for curr in newset if bamfile.getrname(curr.tid) in trnatranscripts)
            locusmaps = list(itertools.chain.from_iterable(trnadata.transcriptdict[bamfile.getrname(curr.tid)] for curr in trnamappings))
            
            if trnamappings[0].get_tag("XM")+trnamappings[0].get_tag("XO") > 0:
                imperfect += 1
            if trnamappings[0].get_tag("XM")+trnamappings[0].get_tag("XO") > 2:
                extraimperfect += 1
            readanticodon = "NNN"
            readamino = "Xxx"
            if len(anticodons - frozenset(['NNN'])) > 1:
                ambanticodon += 1
                acsets[anticodons] += 1
            
            if len(aminos - frozenset(['Und'])) > 1:
                ambamino += 1
                aminosets[aminos] += 1
                
            if diff > 0:
                diffreads += 1
            #finalset = trnareads
            if len(trnamappings) > 1:
                ambtrna += 1
            if setcountfile is not None:
                trnasetcounts[frozenset(bamfile.getrname(curr.tid) for curr in trnamappings)] += 1
            #tags = [("YA",len(anticodons))] + [("YM",len(aminos))]  + [("YR",len(trnamappings))]
            for currtrnamap in trnamappings:
                currtrnamap.tags = currtrnamap.tags + [("YA",len(anticodons))] + [("YM",len(aminos))]  + [("YR",len(trnamappings))] +  [("YL",len(locusmaps))]
            finalset = trnamappings
            
            
            for curr in newset:
                if bamfile.getrname(curr.tid) in trnatranscripts:
                    pass
                    #outfile.write(curr)
                    #curr.tags = curr.tags + readanticodon
                    #finalset.append(curr)
                    pass
                else:
                    #print curr.data["bamline"].rstrip()
                    pass
        elif sum(bamfile.getrname(curr.tid) in extraseqs for curr in newset) > 0:
            for currseqset in seqdata.seqlist():
                if sum(bamfile.getrname(curr.tid) in extraseqdict[currseqset] for curr in newset) > 0:
                    finalset = list(curr for curr in newset if bamfile.getrname(curr.tid) in extraseqdict[currseqset])
                    break
            
        else:
            #for non-tRNA, remove reads that are too small
            if readlength < minnontrnasize:
                continue
            #for non-tRNA, remove sets if mapped too many times
            

            for curr in newset:
                pass
            
                finalset.append(curr)
            #print >>sys.stderr, "len: "+str(len(finalset))    
            if len(finalset) > maxmaps:
                duperemove += 1
                #print >>sys.stderr, "**"
                continue
            if len(newset) > 1:
                nonuniquenontrnas += 1
            else:
                uniquenontrnas += 1
        #print >>sys.stderr,  sum(isprimarymapping(curr) for curr in finalset)
        #This bit is for ensuring that, if I remove the old primary mapping, a new one is chosen
        #Nesecarry for calculating read proportions
        mapsremoved += mappings - len(finalset)
        if sum(isprimarymapping(curr) for curr in finalset) < 1:
            
            
            for i, curr in enumerate(finalset):
                #This
                if i == 0:
                   
                    #print >>sys.stderr, "fixed "+str(len(finalset))+"/"+ str(len(newset))
                    curr.flag &= ~ 0x0100
                    outfile.write(curr)
                else:
                    outfile.write(curr)
        else:
            #print >>sys.stderr, "**"
            for curr in finalset:
                outfile.write(curr)
    '''       
    for curr in acsets.iterkeys():
        if ((1.*acsets[curr])/ambanticodon) > .1:
            print >>logfile, ",".join(curr) + ":"+str(acsets[curr])
    for curr in aminosets.iterkeys():
        if ((1.*aminosets[curr])/ambamino) > .1:
            print >>logfile, ",".join(curr) + ":"+str(aminosets[curr])
    '''
        
    if setcountfile is not None:
        setcounts = open(setcountfile, "w")
        for currset in trnasetcounts:
            print ",".join(currset) +"\t"+str(trnasetcounts[currset])
    #print >>logfile, str(diffreads)+"/"+str(trnareads)
    
    
    print >>logfile, "tRNA Reads with multiple transcripts:"+str(ambtrna)
    print >>logfile, "tRNA Reads with multiple anticodons:"+str(ambanticodon)
    print >>logfile, "tRNA Reads with multiple aminos:"+str(ambamino)
    print >>logfile, "Total tRNA Reads:"+str(trnareads)
    print >>logfile, "Single mapped non-tRNAs:"+str(uniquenontrnas)
    print >>logfile, "Multiply mapped non-tRNAs:"+str(nonuniquenontrnas)
    print >>logfile, "Imperfect matches:"+str(imperfect)+"/"+str(trnareads)
    #print >>logfile, "Extra Imperfect matches:"+str(extraimperfect)+"/"+str(trnareads)
    
    #print >>logfile, str(trnareads)+"/"+str(totalreads)
    #print >>logfile, str(maxreads)+"/"+str(totalreads)
    #print >>logfile, str(multimaps)+"/"+str(totalreads)
    #print >>logfile, str(shortened)+"/"+str(multimaps)
    #print >>logfile, "Mappings Removed:"+str(mapsremoved)+"/"+str(totalmaps)
    outfile.close()
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='get all best tRNA reads')
    parser.add_argument('trnaname',
                    help='name of tRNA database')
    parser.add_argument('--progname',
                       help='program name')
    parser.add_argument('--fqname',
                       help='fastq file name')
    parser.add_argument('--expname',
                       help='library name')

    parser.add_argument('--trnasetcounts',
                       help='Counts for all sets of tRNAs')
    parser.add_argument('--minnontrnasize',type=int,default=20,
                       help='Minimum read length for non-tRNAs')

    
    args = parser.parse_args()
    getbesttrnamappings(args.trnaname, progname = args.progname, fqname = args.fqname, libname = args.expname, setcountfile = args.trnasetcounts, minnontrnasize = args.minnontrnasize)

