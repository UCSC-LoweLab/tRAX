#!/usr/bin/env python3

import pysam
import sys
import argparse
import os.path
from collections import defaultdict
from trnasequtils import *
import itertools
import time


parser = argparse.ArgumentParser(description='Process tRNA experiment.')

parser.add_argument('--samplefile',required=True,
                   help='sample file name')
parser.add_argument('--samplename',required=True,
                   help='new name for sample files')
parser.add_argument('--sizecutoff',
                   help='size cutoff')

args = parser.parse_args()

sizecutoff = 60
samplefilename = args.samplefile
newname = args.samplename
if args.sizecutoff is not None:
    sizecutoff = int(args.sizecutoff)


sampledata = samplefile(samplefilename)

longsamplelabels = list()
shortsamplelabels = list()
for currsample in sampledata.getsamples():
    currbam = sampledata.getbam(currsample)
    bamfile = pysam.Samfile(currbam, "r" )
    newheader = bamfile.header.to_dict()
    shortoutsample = currsample+"_lt"+str(sizecutoff)
    longoutsample = currsample+"_gt"+str(sizecutoff)
    
    
    shortsamplelabels.append(shortoutsample+"\t"+sampledata.getreplicatename(currsample)+"\t"+sampledata.getfastq(currsample))
    longsamplelabels.append(longoutsample+"\t"+sampledata.getreplicatename(currsample)+"\t"+sampledata.getfastq(currsample))

    shortoutfile = pysam.Samfile( shortoutsample+".bam", "wb", header = newheader )
    longoutfile = pysam.Samfile( longoutsample+".bam", "wb", header = newheader )
    

    for currline in bamfile:
        currlength = len(currline.get_forward_sequence())
        #print  (currline.get_forward_sequence(), file = sys.stderr)
        #print  (currlength, file = sys.stderr)
        if currlength > sizecutoff:
            longoutfile.write(currline)
        else:
            shortoutfile.write(currline)
        
longsamplefile = open(newname+"_long.txt", "w")
for currline in longsamplelabels:
    print(currline, file = longsamplefile)
shortsamplefile =  open(newname+"_short.txt", "w") 
for currline in shortsamplelabels:
    print(currline, file = shortsamplefile)
    


'''
samtools view -h "${NAME%%.*}_dedup.bam" | awk -v TS="${TDRSIZE}" 'length($10) <= TS || $1 ~ /^@/' | samtools view -bS - > "${NAME%%.*}_dedup_tdr_${TDRSIZE}.bam"
samtools view -h "${NAME%%.*}_dedup.bam" | awk -v TS="${TDRSIZE}" 'length($10) > TS || $1 ~ /^@/' | samtools view -bS - > "${NAME%%.*}_dedup_fl_${TDRSIZE}.bam"
'''

