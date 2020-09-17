#!/usr/bin/env python

import sys
import subprocess
import argparse
from tempfile import NamedTemporaryFile
import os
import os.path
import re
from trnasequtils import *
from multiprocessing import Pool, cpu_count
import time

MAXMAPS = 100




defaultminnontrnasize = 20

def wrapbowtie2(bowtiedb, unpaired, outfile, scriptdir, trnafile, maxmaps = MAXMAPS,program = 'bowtie2', logfile = None, mapfile = None, expname = None, samplename = None, minnontrnasize = defaultminnontrnasize, numcores = 1):
    '''
    I think the quals are irrelevant due to the RT step, and N should be scored as only slightly better than a mismatch
    this does both
    ignoring quals forces the mismatches to have a single score, which is necessary for precise N weighting
    
     --ignore-quals --np 5
     Very sensitive is necessary due to mismatch caused by modified base misreads
    '''

    bowtiecommand = program+' -x '+bowtiedb+' -k '+str(maxmaps)+' --very-sensitive --ignore-quals --np 5 --reorder -p '+str(numcores)+' -U '+unpaired

    #print >>sys.stderr, bowtiecommand
    temploc = os.path.basename(outfile)
    #bowtiecommand = bowtiecommand + ' | '+scriptdir+'choosemappings.py '+trnafile+' | samtools sort - '+outfile
    bowtiecommand = bowtiecommand + ' | '+scriptdir+'choosemappings.py '+trnafile+' --progname='+"TRAX"+ ' --fqname=' +unpaired+' --expname='+expname + ' --minnontrnasize='+str(minnontrnasize)+' | samtools sort -T '+tempfile.gettempdir()+"/"+temploc+'temp - -o '+outfile+'.bam'
    print >>sys.stderr,  bowtiecommand
    if logfile:
        print >>logfile,  bowtiecommand
        logfile.flush()
    bowtierun = None
    
    bowtierun = subprocess.Popen(bowtiecommand, shell = True, stderr = subprocess.PIPE)

    output = bowtierun.communicate()
    errinfo = output[1]
    if logfile is not None:
        print >>logfile, errinfo 
        logfile.flush()
    if bowtierun.returncode:
        return mapinfo(0,0,0,0, errinfo, samplename, failedrun = True, bowtiecommand = bowtiecommand)


    rereadtotal = re.search(r'(\d+).*reads',errinfo )
    rereadunmap = re.search(r'\s*(\d+).*0 times',errinfo )
    rereadsingle = re.search(r'\s*(\d+).*exactly 1 time',errinfo )
    rereadmult = re.search(r'\s*(\d+).*>1 times',errinfo )
    if rereadtotal and rereadunmap and rereadsingle and rereadmult:
        totalreads = rereadtotal.group(1)
        unmappedreads = rereadunmap.group(1)
        singlemaps = rereadsingle.group(1)
        multmaps = rereadmult.group(1)
        return mapinfo(singlemaps,multmaps,unmappedreads,totalreads, errinfo, samplename, bowtiecommand = bowtiecommand)
        
    else:
        print >>sys.stderr, "Could not map "+unpaired +", check mapstats file"
        print >>sys.stderr, "Exiting..."
        print >>sys.stderr, errinfo
        return mapinfo(0,0,0,0, errinfo, samplename, failedrun = True, bowtiecommand = bowtiecommand)
    

def checkheaders(bamname, fqname):
    try:
        bamfile = pysam.Samfile(bamname, "r" )
    except ValueError:
        return True
    except IOError as e:
        print >>sys.stderr, "Failed to read "+bamname
        print >>sys.stderr, e
        sys.exit(1)
    newheader = bamfile.header
    if len(newheader["PG"]) > 1 and newheader["PG"][1]["PN"] == "TRAX":
        
        if newheader["RG"][0]["ID"] != fqname:
            return False
    return True

#print >>sys.stderr, os.path.dirname(os.path.realpath(sys.argv[0]))
#print >>sys.stderr, os.path.abspath(__file__)

    
class mapinfo:
    def __init__(self, singlemap, multimap, unmap, totalreads, bowtietext, samplename, failedrun = False, bowtiecommand = None):
        self.unmaps = unmap
        self.singlemaps = singlemap
        self.multimaps = multimap
        self.totalreads = totalreads
        self.bowtietext = bowtietext
        self.samplename = samplename
        self.failedrun = failedrun
        self.bowtiecommand = bowtiecommand
        self.unmap = int(self.totalreads) - (int(self.multimaps) + int(self.singlemaps))
    def printbowtie(self, logfile = sys.stderr):
        print >>logfile, "******************************************************************"
        print >>logfile, self.bowtiecommand
        print >>logfile, self.bowtietext
def mapreads(*args, **kwargs):
    return wrapbowtie2(*args, **kwargs)
def mapreadspool(args):
    return mapreads(*args[0], **args[1])
    
def compressargs( *args, **kwargs):
    return tuple([args, kwargs])
def testmain(**argdict):
    argdict = defaultdict(lambda: None, argdict)
    scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"
    samplefilename = argdict["samplefile"]
    
    sampledata = samplefile(argdict["samplefile"])
    trnafile = argdict["trnafile"]
    logfile = argdict["logfile"]
    mapfile = argdict["mapfile"]
    bowtiedb = argdict["bowtiedb"]
    lazycreate = argdict["lazy"]
    minnontrnasize = argdict["minnontrnasize"]
    
    if "cores" in argdict:
        cores = int(argdict["cores"])
    else:
        cores = min(8,cpu_count())
    #sys.exit()
    print >>sys.stderr,"cores: "+str(cores) 
    workingdir = './'
    #samplefile = open(args.samplefile)
    
    samples = sampledata.getsamples()
    
    trnafile = trnafile
    print >>sys.stderr, "logging to "+logfile
    if logfile and lazycreate:
        logfile = open(logfile,'a')
        print >>logfile, "New mapping"
    elif logfile:
        logfile = open(logfile,'w')
    else:
        logfile = sys.stderr

    unmaps = defaultdict(int)
    singlemaps = defaultdict(int)
    multimaps = defaultdict(int)
    totalreads = defaultdict(int)
    
    if not os.path.isfile(bowtiedb+".fa"):
        print >>sys.stderr, "No bowtie2 database "+bowtiedb
        sys.exit(1)
    badsamples = list()
    for samplename in samples:
        bamfile = workingdir+samplename
        
        if lazycreate and os.path.isfile(bamfile+".bam"):   
            if not checkheaders(bamfile+".bam", sampledata.getfastq(samplename)):
                badsamples.append(bamfile+".bam")

                
            
        else:
            if os.path.isfile(bamfile+".bam"):

                if not checkheaders(bamfile+".bam", sampledata.getfastq(samplename)):
                    badsamples.append(bamfile+".bam")
    
    if len(badsamples) > 0:
        print >>sys.stderr, "Bam files "+",".join(badsamples)+" does not match fq files"
        print >>sys.stderr, "Aborting"
        sys.exit(1)               
    #'samtools sort -T '+tempfile.gettempdir()+"/"+outfile+'temp - -o '+outfile+'.bam'
    tempfilesover = list()
    missingfqfiles = list()
    for samplename in samples:
        #redundant but ensures compatibility
        bamfile = workingdir+samplename
        temploc = os.path.basename(bamfile)
        #print >>sys.stderr, "***"
        #print >>sys.stderr, samplename+'temp'
        
        for currfile in os.listdir(tempfile.gettempdir()):
            #
            if currfile.startswith(samplename+'temp'):
                tempfilesover.append(currfile)
        fqfile = sampledata.getfastq(samplename)
        if not os.path.isfile(fqfile):
            missingfqfiles.append(fqfile)
    if len(tempfilesover) > 0:
        for currfile in tempfilesover:
            print >>sys.stderr, tempfile.gettempdir() +"/"+ currfile + " temp bam files exists"
        print >>sys.stderr, "these files must be deleted to proceed"
        sys.exit(1)
    if len(missingfqfiles) > 0:
        print >>sys.stderr, ",".join(missingfqfiles) + " fastq files missing"
        sys.exit(1)
    mapresults = dict()
    multithreaded = True
    if multithreaded:
        mapargs = list()
        print >>sys.stderr, cores
        mappool = Pool(processes=cores)
        mapsamples = list()
        for samplename in samples:
            bamfile = workingdir+samplename
            
            if lazycreate and os.path.isfile(bamfile+".bam"):
                pass

                print >>sys.stderr, "Skipping "+samplename

            else:
                mapargs.append(compressargs(bowtiedb, sampledata.getfastq(samplename),bamfile,scriptdir, trnafile, expname = samplefilename, samplename = samplename, minnontrnasize = minnontrnasize))
                
                
                #mapresults[samplename] = mapreads(bowtiedb, sampledata.getfastq(samplename),bamfile,scriptdir, trnafile,  logfile=logfile, expname = samplefilename)
                mapsamples.append(samplename)
        #results = mappool.map(mapreadspool, mapargs)
        starttime = time.time()
        for currresult in mappool.imap_unordered(mapreadspool, mapargs):
            #print >>sys.stderr, "time "+currresult.samplename+": "+str(time.time() - starttime)
            if currresult.failedrun == True:
                print >>sys.stderr, "Failure to Bowtie2 map"
                #print >>sys.stderr, output[1]
                currresult.printbowtie(logfile)
                sys.exit(1)
            mapresults[currresult.samplename] = currresult
            currresult.printbowtie(logfile)
                
    else:
        for samplename in samples:
            bamfile = workingdir+samplename
            
            if lazycreate and os.path.isfile(bamfile+".bam"):
                pass
                    
                print >>sys.stderr, "Skipping "+samplename
                
            else:
        
                mapresults[samplename] = mapreads(bowtiedb, sampledata.getfastq(samplename),bamfile,scriptdir, trnafile,  logfile=logfile, expname = samplefilename, minnontrnasize = minnontrnasize)

    if lazycreate:
        #here is where I might add stuff to read old files in lazy mode
        pass
    if mapfile is not None and not lazycreate:
        mapinfo = open(mapfile,'w')                
        print >>mapinfo, "\t".join(samples)
        print >>mapinfo, "unmap\t"+"\t".join(mapresults[currsample].unmaps for currsample in samples)
        print >>mapinfo, "single\t"+"\t".join(mapresults[currsample].singlemaps for currsample in samples)
        print >>mapinfo, "multi\t"+"\t".join(mapresults[currsample].multimaps for currsample in samples)
        #print >>mapinfo, "total\t"+"\t".join(totalreads[currsample] for currsample in samples)
        mapinfo.close()
    
        
        
        #print >>logfile, "Processing "+samplename +" mappings"
    logfile.close()
    
        
        
        #result = subprocess.call(scriptdir+'choosemappings.py '+trnafile+' <'+bamfile +' | samtools view -F 4 -b - | samtools sort - '+workingdir+samplename+'_sort', shell = True)
def main(**argdict):
    argdict = defaultdict(lambda: None, argdict)
    scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"
    samplefilename = argdict["samplefile"]
    
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
    print >>sys.stderr, "logging to "+logfile
    if logfile and lazycreate:
        logfile = open(logfile,'a')
        print >>logfile, "New mapping"
    elif logfile:
        logfile = open(logfile,'w')
    else:
        logfile = sys.stderr

    unmaps = defaultdict(int)
    singlemaps = defaultdict(int)
    multimaps = defaultdict(int)
    totalreads = defaultdict(int)
    
    if not os.path.isfile(bowtiedb+".fa"):
        print >>sys.stderr, "No bowtie2 database "+bowtiedb
        sys.exit(1)
    badsamples = list()
    for samplename in samples:
        bamfile = workingdir+samplename
        
        if lazycreate and os.path.isfile(bamfile+".bam"):   
            if not checkheaders(bamfile+".bam", sampledata.getfastq(samplename)):
                badsamples.append(bamfile+".bam")

                
            
        else:
            if os.path.isfile(bamfile+".bam"):

                if not checkheaders(bamfile+".bam", sampledata.getfastq(samplename)):
                    badsamples.append(bamfile+".bam")
    
    if len(badsamples) > 0:
        print >>sys.stderr, "Bam files "+",".join(badsamples)+" does not match fq files"
        print >>sys.stderr, "Aborting"
        sys.exit(1)               
    #'samtools sort -T '+tempfile.gettempdir()+"/"+outfile+'temp - -o '+outfile+'.bam'
    tempfilesover = list()
    missingfqfiles = list()
    for samplename in samples:
        #redundant but ensures compatibility
        bamfile = workingdir+samplename
        temploc = os.path.basename(bamfile)
        #print >>sys.stderr, "***"
        #print >>sys.stderr, samplename+'temp'
        
        for currfile in os.listdir(tempfile.gettempdir()):
            #
            if currfile.startswith(samplename+'temp'):
                tempfilesover.append(currfile)
        fqfile = sampledata.getfastq(samplename)
        if not os.path.isfile(fqfile):
            missingfqfiles.append(fqfile)
    if len(tempfilesover) > 0:
        for currfile in tempfilesover:
            print >>sys.stderr, tempfile.gettempdir() +"/"+ currfile + " temp bam files exists"
        print >>sys.stderr, "these files must be deleted to proceed"
        sys.exit(1)
    if len(missingfqfiles) > 0:
        print >>sys.stderr, ",".join(missingfqfiles) + " fastq files missing"
        sys.exit(1)
    mapresults = dict()    
    for samplename in samples:
        bamfile = workingdir+samplename
        
        if lazycreate and os.path.isfile(bamfile+".bam"):
            pass
                
            print >>sys.stderr, "Skipping "+samplename
            
        else:
            if os.path.isfile(bamfile+".bam"):
                pass

            print >>logfile, "Mapping "+samplename
            print >>sys.stderr, "Mapping "+samplename
            logfile.flush()
            mapresults[samplename] = wrapbowtie2(bowtiedb, sampledata.getfastq(samplename),bamfile,scriptdir, trnafile,  logfile=logfile, expname = samplefilename, samplename = samplename)


    if lazycreate:
        #here is where I might add stuff to read old files in lazy mode
        pass
    if mapfile is not None and not lazycreate:
        mapinfo = open(mapfile,'w')                
        print >>mapinfo, "\t".join(samples)
        print >>mapinfo, "unmap\t"+"\t".join(mapresults[currsample].unmaps for currsample in samples)
        print >>mapinfo, "single\t"+"\t".join(mapresults[currsample].singlemaps for currsample in samples)
        print >>mapinfo, "multi\t"+"\t".join(mapresults[currsample].multimaps for currsample in samples)
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


