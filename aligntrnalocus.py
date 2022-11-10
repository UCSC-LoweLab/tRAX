#!/usr/bin/env python3

import pysam
import os
import sys
import argparse
import subprocess
from collections import defaultdict
from trnasequtils import *






def main(**args):
    args = defaultdict(lambda: None, args)
    stkfile = args["stkfile"]
    genomefile = os.path.expanduser(args["genomefile"])
    trnalocifile = args["trnaloci"]
    
    scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"        
    trnacmfile = args["cmmodel"]

    
    trnaloci = list(readbed(trnalocifile, orgdb = "genome", seqfile=genomefile))
    lociseqs = getseqdict(trnaloci, faifiles = {"genome":genomefile+".fai"})
    #print lociseqs
    #lociseqfile = tempmultifasta(lociseqs)
    devnull = open(os.devnull, 'w')
    seqfile = tempmultifasta(iter(lociseqs.items()))
    cmcommand = ['cmalign', "-o", stkfile,"--nonbanded","--notrunc", "-g",trnacmfile,seqfile.name]
    #print >>sys.stderr, " ".join(cmcommand)
    cmrun = subprocess.Popen(cmcommand, stdout = devnull)
    result = cmrun.wait()
    if result:
        print("Failure to align tRNAs", file=sys.stderr)
        sys.exit(1)
    seqfile.close()
    devnull.close()
    
    #trnaalign
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
    parser.add_argument('--trnaloci',
                       help='bed file of tRNA loci')      
    parser.add_argument('--genomefile',
                       help='fasta file of genome')
    parser.add_argument('--stkfile',
                       help='stockholm output file')
    parser.add_argument('--cmmodel', 
                       help='covariance model to use')
    args = vars(parser.parse_args())
    main(args)    

        