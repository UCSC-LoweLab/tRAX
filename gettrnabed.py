#!/usr/bin/env python3

import re
import random
import math
import os
import subprocess
import tempfile
import sys
import itertools
from collections import defaultdict
import argparse
from parsetrnas import *


#This program gets the mature tRNA sequences

parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
parser.add_argument('--trnascan',  nargs='+', default=list(),
                   help='tRNAscan-SE file')
parser.add_argument('--rnacentral', nargs='+', default=list(),
                   help='RNAcentral tRNA file')
parser.add_argument('--genome', 
                   help='fasta sequence of genome')
parser.add_argument('--chromtranslate', 
                   help='translation table of chromosome names')

parser.add_argument('--tag', nargs='1',
                   help='tag to be added to tRNA name')

args = parser.parse_args()







alltrnas = list()
trnascantrnas = list()
for currfile in args.trnascan:
    trnascantrnas.extend(readtRNAscan(currfile, args.genome))
trnacentraltrnas = list()
for currfile in args.rnacentral:
    trnacentraltrnas.extend(readrnacentral(currfile, args.chromtranslate))
    
#alltrnas = list(getuniquetRNAs(trnascantrnas)) + trnacentraltrnas

alltrnas = trnascantrnas+trnacentraltrnas
                


anticodoncount = defaultdict(int)
trnanames = dict()
trnalist = list()
for currtrans in alltrnas:
    if currtrans.name is None:
        name = 'tRNA-'+currtrans.amino + currtrans.anticodon+ str(anticodoncount[currtrans.anticodon]+ 1)
        anticodoncount[currtrans.anticodon] += 1
    else:
        #print >>sys.stderr, "***"
        name = currtrans.name
    #trnanames[name] = currtrans
    trnalist.append(name)
    #print >>sys.stderr, name
    print(currtrans.loc.bedstring())

#sys.exit()
trnamods = dict()
allmods = set()
nomatches = 0
trnamismatches = dict()
