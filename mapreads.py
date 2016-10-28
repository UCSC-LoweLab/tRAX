#!/usr/bin/env python

import sys
import subprocess
import argparse
from tempfile import NamedTemporaryFile
import os
import os.path
import re
from trnasequtils import *

MAXMAPS = 100

'''
~/pythonsource/trnaseq/mapreads.py --samplefile=/projects/lowelab/users/holmes/pythonsource/trnatest/testcomp.txt --bowtiedb=hg19-trnagenome >counts.txt

~/pythonsource/trnaseq/mapreads.py --samplefile=agingshort.txt --trnafile=sacCer3-trnatable.txt --bowtiedb=sacCer3-trnagenome >aging-counts.txt

cat <(samtools view -H HBV10_1.bam) <(samtools view HBV10_1.bam | grep UNC11-SN627:301:C21RKACXX:3:1101:20489:57142) | ~/pythonsource/trnaseq/choosemappings.py hg19-trnatable.txt | samtools view -

'''

numcores = 4



def wrapbowtie2(bowtiedb, unpaired, outfile, scriptdir, trnafile, maxmaps = MAXMAPS,program = 'bowtie2', logfile = None, mapfile = None):
    '''
    I think the quals are irrelevant, and N should be scored as only slightly better than a mismatch
    this does both
    ignoring quals forces the mismatches to have a single score, which is necessary for precise N weighting
    
     --ignore-quals --np 5
     Very sensitive is necessary due to mismatch caused by modified base misreads
    '''

    bowtiecommand = program+' -x '+bowtiedb+' -k '+str(maxmaps)+' --very-sensitive --ignore-quals --np 5 --reorder -p '+str(numcores)+' -U '+unpaired

    #print >>sys.stderr, bowtiecommand
    
    bowtiecommand = bowtiecommand + ' | '+scriptdir+'choosemappings.py '+trnafile+' | samtools sort - '+outfile
    if logfile:
        print >>logfile,  bowtiecommand
    bowtierun = None
    logfile.flush()
    bowtierun = subprocess.Popen(bowtiecommand, shell = True, stderr = subprocess.PIPE)

    output = bowtierun.communicate()
    errinfo = output[1]
    if logfile is not None:
        print >>logfile, errinfo 
    logfile.flush()
    if bowtierun.returncode:
        print >>sys.stderr, "Failure to Bowtie2 map"
       
        sys.exit(1)

    rereadtotal = re.search(r'(\d+).*reads',errinfo )
    rereadunmap = re.search(r'\s*(\d+).*0 times',errinfo )
    rereadsingle = re.search(r'\s*(\d+).*exactly 1 time',errinfo )
    rereadmult = re.search(r'\s*(\d+).*>1 times',errinfo )
    if rereadtotal and rereadunmap and rereadsingle and rereadmult:
        totalreads = rereadtotal.group(1)
        unmappedreads = rereadunmap.group(1)
        singlemaps = rereadsingle.group(1)
        multmaps = rereadmult.group(1)
        return [unmappedreads,singlemaps,multmaps,totalreads]
        
    else:
        print >>sys.stderr, "Could not extract mapping information for "+unpaired
        print >>sys.stderr, "Exiting..."
        sys.exit(1)
        return None
    


#print >>sys.stderr, os.path.dirname(os.path.realpath(sys.argv[0]))
#print >>sys.stderr, os.path.abspath(__file__)
def main(**argdict):
    argdict = defaultdict(lambda: None, argdict)
    scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"
    sampledata = samplefile(argdict["samplefile"])
    trnafile = argdict["trnafile"]
    logfile = argdict["logfile"]
    mapfile = argdict["mapfile"]
    bowtiedb = argdict["bowtiedb"]
    lazycreate = argdict["lazy"]
    #sys.exit()
    
    workingdir = './'
    #samplefile = open(args.samplefile)
    
    samples = sampledata.getsamples()
    
    trnafile = trnafile
        
    if logfile:
        logfile = open(logfile,'a')
    else:
        logfile = sys.stderr

    unmaps = defaultdict(int)
    singlemaps = defaultdict(int)
    multimaps = defaultdict(int)
    totalreads = defaultdict(int)
    
    if not os.path.isfile(bowtiedb+".fa"):
        print >>sys.stderr, "No bowtie2 database "+bowtiedb
        sys.exit(1)
    for samplename in samples:
        bamfile = workingdir+samplename
        if lazycreate and os.path.isfile(bamfile+".bam"):
            print >>sys.stderr, "Skipping "+samplename
        else:
            print >>logfile, "Mapping "+samplename
            print >>sys.stderr, "Mapping "+samplename
            logfile.flush()
            mapresults = wrapbowtie2(bowtiedb, sampledata.getfastq(samplename),bamfile,scriptdir, trnafile,  logfile=logfile)

            
            if mapresults is not None:
                unmaps[samplename] = mapresults[0]
                singlemaps[samplename] = mapresults[1]
                multimaps[samplename] = mapresults[2]
                totalreads[samplename] = mapresults[3]
    
    if lazycreate:
        #here is where I might add stuff to read old files in lazy mode
        pass
    if mapfile is not None and not lazycreate:
        mapinfo = open(mapfile,'w')                
        print >>mapinfo, "\t".join(samples)
        print >>mapinfo, "unmap\t"+"\t".join(unmaps[currsample] for currsample in samples)
        print >>mapinfo, "single\t"+"\t".join(singlemaps[currsample] for currsample in samples)
        print >>mapinfo, "multi\t"+"\t".join(multimaps[currsample] for currsample in samples)
        #print >>mapinfo, "total\t"+"\t".join(totalreads[currsample] for currsample in samples)
        mapinfo.close()
    
        
        
        #print >>logfile, "Processing "+samplename +" mappings"
    logfile.close()
    
        
        
        #result = subprocess.call(scriptdir+'choosemappings.py '+trnafile+' <'+bamfile +' | samtools view -F 4 -b - | samtools sort - '+workingdir+samplename+'_sort', shell = True)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Map reads with bowtie2 and process mappings')
    
    parser.add_argument('--samplefile',
                       help='Sample file in format')
    parser.add_argument('--trnafile',
                       help='tRNA file in format')
    parser.add_argument('--logfile',
                       help='log file for error messages and mapping stats')
    parser.add_argument('--mapfile',
                       help='output table with mapping stats')
    parser.add_argument('--bowtiedb',
                       help='Location of Bowtie 2 database')
    parser.add_argument('--lazy', action="store_true", default=False,
                       help='do not remap if mapping results exist')
    
    args = parser.parse_args()
    main(samplefile = args.samplefile, trnafile= args.trnafile, logfile = args.logfile, bowtiedb = args.bowtiedb, lazy = args.lazy, mapfile = args.mapfile)


