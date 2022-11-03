#!/usr/bin/env python

import pysam
import sys
import argparse
import string
import itertools
from collections import defaultdict
import os.path
from trnasequtils import *
import subprocess

        
inputbam = sys.argv[1]
dbname = sys.argv[2]


    

bamfile = pysam.Samfile(inputbam, "rb" )

outfile = pysam.Samfile( "-", "w", template = bamfile )
#outfile = pysam.Samfile( outputbam, "wb", template = bamfile )

 
trnainfo = transcriptfile(dbname+"-trnatable.txt")

trnaloci = getnamedict(readbed(dbname+"-trnaloci.bed", includeintrons = True))


trnatranscripts = trnainfo.gettranscripts() 
npad = 20


bam_cins = 1
bam_cdel = 2


'''
566645416..
566645491

chr13	23413247	23413320	tRNA-Val-AAC-2-2	1000	-
chr13	23285400	23285473	tRNA-Val-AAC-2-1	1000	+


Doesnt quite work with multiple bam indels

CIGAR strings are relative to the genome strand

chr11	59318766	59318852	tRNA-Arg-TCT-3-2	1000	+	59318766	59318852	0	2	37,36	0,50
chr9	131102354	131102445	tRNA-Arg-TCT-3-1	1000	-	131102354	131102445	0	2	38,35	0,56


Actual from the site:
chr9	131102354	131102445	tRNA-Arg-TCT-3-1	1000	-	131102354	131102445	0	2	36,37,	0,54,



'''

def cigarlength(cigar):
    return sum(curr[1] for curr in cigar if curr[0] in set([0,bam_cins]))
    
def reverseintrons(introns, length):
    for currint in reversed(introns):
        yield tuple([length - currint[0] + 1, currint[1]])

cigarset = set()
uniquemode = False
for currmap in bamfile:
    chromname = bamfile.getrname(currmap.tid)
    readquals = currmap.query_qualities
    readseq = currmap.query_sequence
    readtags = currmap.get_tags() 
    origcigar = currmap.cigartuples
    if chromname in trnatranscripts:

        origstart = currmap.reference_start	
        for currlocusname in trnainfo.transcriptdict[chromname]:
            currlocusmap = currmap
            currlocus = trnaloci[currlocusname]
            currlocusmap.reference_id = bamfile.get_tid(currlocus.chrom)
            
            introns = list()
            if currlocus.data["blockcount"] > 1:
                lastblock = 0
                for i in range(currlocus.data["blockcount"]):
                    currblocksize = int(currlocus.data["blocksizes"][i])
                    currblockstart = int(currlocus.data["blockstarts"][i])
                    if lastblock != 0:
                        introns.append([lastblock,currblockstart - lastblock])
                    lastblock += currblocksize
            if len(origcigar) > 1:
                #continue 
                pass
            if currlocus.strand == '+':
                currmap.reference_start = origstart - npad + currlocus.start
                currmap.query_sequence = readseq
                currmap.query_qualities = readquals
                currmap.is_reverse = False
                
                currpoint = origstart - npad
                newcigar = list()
                for i in range(len(introns)): 
                    if currpoint >= introns[i][0]:
                        currmap.reference_start += introns[i][1]
                        currpoint += introns[i][1]
                for currcigar in origcigar:
                    #print >>sys.stderr, currcigar  
                    cigarset.add(currcigar[0])
                    if currcigar[0] in set([0,bam_cins]):
                        foundintron = False
                        for intronstart, intronlength in introns:
                            if currpoint < intronstart <  currpoint + currcigar[1]:
                                firseglength = intronstart  - currpoint
                                secseglength = currpoint + currcigar[1] - intronstart
                                newcigar.append(tuple([currcigar[0],firseglength]))
                                newcigar.append(tuple([bam_cdel,intronlength]))
                                newcigar.append(tuple([currcigar[0],secseglength]))
                                foundintron = True
                        if not foundintron:
                            newcigar.append(currcigar)
                        if currcigar[0] in set([0,bam_cdel]):
                            currpoint += currcigar[1]
                    else:
                        newcigar.append(currcigar)
                
                if cigarlength(newcigar) != cigarlength(origcigar):        
                    print(origcigar, file=sys.stderr)
                    print(newcigar, file=sys.stderr)
                    print("**||||", file=sys.stderr)
                    #sys.exit(1)
                    pass
                currmap.cigartuples = newcigar
            else:
                #continue
                currmap.reference_start =  currlocus.end - (origstart - npad + cigarlength(origcigar) + sum(curr[1] for curr in introns)) 
                currmap.query_sequence   = revcom(readseq)
                currmap.query_qualities = list(reversed(readquals))
                currmap.is_reverse = True
                
                currpoint = origstart - npad
                newcigar = list()
                        
                for intronstart, intronlength  in reverseintrons(introns, currlocus.length()): 
                    if currpoint + cigarlength(origcigar) < intronstart:
                        currmap.reference_start += intronlength
                        currpoint -= intronlength
                        pass
                for currcigar in reversed(origcigar):
                    cigarset.add(currcigar[0])
                    if currcigar[0] in set([0,bam_cins]):
                        foundintron = False
                        for intronstart, intronlength in introns:
                            if currpoint  < intronstart <  currpoint + currcigar[1]:
                                firseglength = intronstart  - currpoint + 1
                                secseglength = currpoint + currcigar[1] - intronstart - 1 
                                newcigar.append(tuple([currcigar[0],secseglength]))
                                newcigar.append(tuple([bam_cdel,intronlength]))
                                newcigar.append(tuple([currcigar[0],firseglength]))
                                foundintron = True
                                
                        if not foundintron:
                            newcigar.append(currcigar)
                        if currcigar[0] in set([0,bam_cdel]):
                            currpoint += currcigar[1]
                    else:
                        newcigar.append(currcigar)
                if cigarlength(newcigar) != cigarlength(origcigar):        
                    #print >>sys.stderr, currmap
                    print("neg", file=sys.stderr)
                    print(origcigar, file=sys.stderr)
                    print(newcigar, file=sys.stderr)
                    print("**||", file=sys.stderr)
                    #sys.exit(1)
                    pass
                currmap.cigartuples = newcigar
                
            currmap.set_tags(readtags)
            outfile.write(currmap)
        
    else:
        pass
        outfile.write(currmap)

