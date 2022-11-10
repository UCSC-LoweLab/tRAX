#!/usr/bin/env python3

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

class trnainfo:
    def __init__(self,multtrans,multac,multamino ,trna ,singlenon,multiplenon):
        
        self.multtrans    = int(multtrans)
        self.multac       = int(multac)
        self.multamino    = int(multamino)
        self.trna         = int(trna)
        self.singlenon    = int(singlenon)
        self.multiplenon  = int(multiplenon )
        self.multitrna = (self.multtrans + self.multac + self.multamino)
        self.singletrna = self.trna - self.multitrna
        
    def uniquereads(self):
       return  self.singletrna + self.singlenon
    def nonuniquereads(self):
       return self.multitrna + self.multiplenon
       
def wrapbowtie2(bowtiedb, unpaired, outfile, scriptdir, trnafile, maxmaps = MAXMAPS,program = 'bowtie2', logfile = None, mapfile = None, expname = None, samplename = None, minnontrnasize = defaultminnontrnasize, numcores = 1, local = False):
    '''
    I think the quals are irrelevant due to the RT step, and N should be scored as only slightly better than a mismatch
    this does both
    ignoring quals forces the mismatches to have a single score, which is necessary for precise N weighting
    
     --ignore-quals --np 5
     Very sensitive is necessary due to mismatch caused by modified base misreads
    ''' 
    localmode = " "
    if local:
        localmode = " --local "
    bowtiecommand = program+localmode+' -x '+bowtiedb+' -k '+str(maxmaps)+' --very-sensitive --ignore-quals --np 5 --reorder -p '+str(numcores)+' -U '+unpaired

    #print >>sys.stderr, bowtiecommand
    temploc = os.path.basename(outfile)
    #bowtiecommand = bowtiecommand + ' | '+scriptdir+'choosemappings.py '+trnafile+' | samtools sort - '+outfile
    bowtiecommand = bowtiecommand + ' | '+scriptdir+'choosemappings.py '+trnafile+' --progname='+"TRAX"+ ' --fqname=' +unpaired+' --expname='+expname + ' --minnontrnasize='+str(minnontrnasize)+' | samtools sort -T '+tempfile.gettempdir()+"/"+temploc+'temp - -o '+outfile+'.bam'
    print(bowtiecommand, file=sys.stderr)
    if logfile:
        print(bowtiecommand, file=logfile)
        logfile.flush()
    bowtierun = None
    
    bowtierun = subprocess.Popen(bowtiecommand, shell = True, stderr = subprocess.PIPE, universal_newlines=True)

    output = bowtierun.communicate()
    errinfo = output[1]
    if logfile is not None:
        print(errinfo, file=logfile) 
        logfile.flush()
    if bowtierun.returncode:
        return mapinfo(0,0,0,0, errinfo, samplename, failedrun = True, bowtiecommand = bowtiecommand)

    '''
    tRNA Reads with multiple transcripts:1637/1942451
    tRNA Reads with multiple anticodons:46/1942451
    tRNA Reads with multiple aminos:33/1942451
    Single mapped non-tRNAs:50865
    Multiply mapped non-tRNAs:313735
    Imperfect matches:1931867/1942451
    '''

    rereadmulttrans = re.search(r'tRNA Reads with multiple transcripts:(\d+)',errinfo )
    rereadmultac = re.search(r'tRNA Reads with multiple anticodons:(\d+)',errinfo )
    rereadmultamino = re.search(r'tRNA Reads with multiple aminos:(\d+)',errinfo )
    rereadtrna = re.search(r'Total tRNA Reads:(\d+)',errinfo )
    rereadsinglenon = re.search(r'Single mapped non-tRNAs:(\d+)',errinfo )
    rereadmultiplenon = re.search(r'Multiply mapped non-tRNAs:(\d+)',errinfo )
    trnamapinfo = None
    if rereadmulttrans and rereadmultac and rereadmultamino and rereadtrna and rereadsinglenon and rereadmultiplenon:
         trnamapinfo = trnainfo(rereadmulttrans.group(1),rereadmultac.group(1),rereadmultamino.group(1) ,rereadtrna.group(1) ,rereadsinglenon.group(1),rereadmultiplenon.group(1) ) 


    rereadtotal = re.search(r'(\d+).*reads',errinfo )
    rereadunmap = re.search(r'\s*(\d+).*0 times',errinfo )
    rereadsingle = re.search(r'\s*(\d+).*exactly 1 time',errinfo )
    rereadmult = re.search(r'\s*(\d+).*>1 times',errinfo )
    if rereadtotal and rereadunmap and rereadsingle and rereadmult:
        totalreads = rereadtotal.group(1)
        unmappedreads = rereadunmap.group(1)
        singlemaps = rereadsingle.group(1)
        multmaps = rereadmult.group(1)
        return mapinfo(singlemaps,multmaps,unmappedreads,totalreads, errinfo, samplename, bowtiecommand = bowtiecommand, trnamapinfo = trnamapinfo)
        
    else:
        print("Could not map "+unpaired +", check mapstats file", file=sys.stderr)
        print("Exiting...", file=sys.stderr)
        print(errinfo, file=sys.stderr)
        return mapinfo(0,0,0,0, errinfo, samplename, failedrun = True, bowtiecommand = bowtiecommand)
    

def checkheaders(bamname, fqname):
    try:
        bamfile = pysam.Samfile(bamname, "r" )
    except ValueError:
        return True
    except IOError as e:
        print("Failed to read "+bamname, file=sys.stderr)
        print(e, file=sys.stderr)
        sys.exit(1)
    newheader = bamfile.header
    if len(newheader["PG"]) > 1 and newheader["PG"][1]["PN"] == "TRAX":
        
        if newheader["RG"][0]["ID"] != fqname:
            return False
    return True

#print >>sys.stderr, os.path.dirname(os.path.realpath(sys.argv[0]))
#print >>sys.stderr, os.path.abspath(__file__)

    
class mapinfo:
    def __init__(self, singlemap, multimap, unmap, totalreads, bowtietext, samplename, failedrun = False, bowtiecommand = None, trnamapinfo = None):
        self.unmaps = unmap
        self.bowtiesinglemaps = singlemap
        self.bowtiemultimaps = multimap
        self.totalreads = totalreads
        self.bowtietext = bowtietext
        self.samplename = samplename
        self.failedrun = failedrun
        self.bowtiecommand = bowtiecommand
        
        
        self.trnamapinfo = trnamapinfo
        if trnamapinfo is not None:
            self.singlemaps = trnamapinfo.uniquereads()
            self.multimaps = trnamapinfo.nonuniquereads()
        else: #has to use the raw bowtie2 output if no tRNA data
            self.singlemaps = singlemap
            self.multimaps = multimap
        self.unmap = int(self.totalreads) - (int(self.multimaps) + int(self.singlemaps))
    def printbowtie(self, logfile = sys.stderr):
        print("******************************************************************", file=logfile)
        print(self.bowtiecommand, file=logfile)
        print(self.bowtietext, file=logfile)
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
    bamdir = argdict["bamdir"]
    local = argdict["local"]
    
    trnamapfile = argdict["trnamapfile"]
    if bamdir is None:
        bamdir = "./"
    
    if "cores" in argdict:
        cores = int(argdict["cores"])
    else:
        cores = min(8,cpu_count())
    #sys.exit()
    print("cores: "+str(cores), file=sys.stderr)
    workingdir = bamdir
    #samplefile = open(args.samplefile)
    
    samples = sampledata.getsamples()
    
    trnafile = trnafile
    print("logging to "+logfile, file=sys.stderr)
    if logfile and lazycreate:
        logfile = open(logfile,'a')
        print("New mapping", file=logfile)
    elif logfile:
        logfile = open(logfile,'w')
    else:
        logfile = sys.stderr

    unmaps = defaultdict(int)
    singlemaps = defaultdict(int)
    multimaps = defaultdict(int)
    totalreads = defaultdict(int)
    
    if not os.path.isfile(bowtiedb+".fa"):
        print("No bowtie2 database "+bowtiedb, file=sys.stderr)
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
        print("Bam files "+",".join(badsamples)+" does not match fq files", file=sys.stderr)
        print("Aborting", file=sys.stderr)
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
            print(tempfile.gettempdir() +"/"+ currfile + " temp bam files exists", file=sys.stderr)
        print("these files must be deleted to proceed", file=sys.stderr)
        sys.exit(1)
    if len(missingfqfiles) > 0:
        print(",".join(missingfqfiles) + " fastq files missing", file=sys.stderr)
        sys.exit(1)
    mapresults = dict()
    multithreaded = True
    if multithreaded:
        mapargs = list()
        print(cores, file=sys.stderr)
        mappool = Pool(processes=cores)
        mapsamples = list()
        for samplename in samples:
            bamfile = workingdir+samplename
            
            if lazycreate and os.path.isfile(bamfile+".bam"):
                pass

                print("Skipping "+samplename, file=sys.stderr)

            else:
                mapargs.append(compressargs(bowtiedb, sampledata.getfastq(samplename),bamfile,scriptdir, trnafile, expname = samplefilename, samplename = samplename, minnontrnasize = minnontrnasize, local = local))
                
                
                #mapresults[samplename] = mapreads(bowtiedb, sampledata.getfastq(samplename),bamfile,scriptdir, trnafile,  logfile=logfile, expname = samplefilename)
                mapsamples.append(samplename)
        #results = mappool.map(mapreadspool, mapargs)
        starttime = time.time()
        for currresult in mappool.imap_unordered(mapreadspool, mapargs):
            #print >>sys.stderr, "time "+currresult.samplename+": "+str(time.time() - starttime)
            if currresult.failedrun == True:
                print("Failure to Bowtie2 map", file=sys.stderr)
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
                    
                print("Skipping "+samplename, file=sys.stderr)
                
            else:
        
                mapresults[samplename] = mapreads(bowtiedb, sampledata.getfastq(samplename),bamfile,scriptdir, trnafile,  logfile=logfile, expname = samplefilename, minnontrnasize = minnontrnasize, local = local)

    if lazycreate:
        #here is where I might add stuff to read old files in lazy mode
        pass
    if mapfile is not None and not lazycreate:
        mapinfo = open(mapfile,'w')                
        print("\t".join(samples), file=mapinfo)
        print("unmap\t"+"\t".join(str(mapresults[currsample].unmaps) for currsample in samples), file=mapinfo)
        print("single\t"+"\t".join(str(mapresults[currsample].singlemaps) for currsample in samples), file=mapinfo)
        print("multi\t"+"\t".join(str(mapresults[currsample].multimaps) for currsample in samples), file=mapinfo)
        mapinfo.close()
        
    if trnamapfile is not None and not lazycreate:
        trnamapinfo = open(trnamapfile,'w')      
        
        print("\t".join(samples), file=trnamapinfo)
        print("multi_nontRNA\t"+"\t".join(str(mapresults[currsample].trnamapinfo.multiplenon) for currsample in samples), file=trnamapinfo)
        print("unique_nontRNA\t"+"\t".join(str(mapresults[currsample].trnamapinfo.singlenon) for currsample in samples), file=trnamapinfo)
        print("multi_amino\t"+"\t".join(str(mapresults[currsample].trnamapinfo.multamino) for currsample in samples), file=trnamapinfo)
        print("unique_amino\t"+"\t".join(str(mapresults[currsample].trnamapinfo.multac) for currsample in samples), file=trnamapinfo)
        print("unique_anticodon\t"+"\t".join(str(mapresults[currsample].trnamapinfo.multtrans) for currsample in samples), file=trnamapinfo)
        print("unique_tRNA\t"+"\t".join(str(mapresults[currsample].trnamapinfo.singletrna) for currsample in samples), file=trnamapinfo)
        


        #print >>mapinfo, "total\t"+"\t".join(totalreads[currsample] for currsample in samples)
        trnamapinfo.close()
        
        
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
    print("logging to "+logfile, file=sys.stderr)
    if logfile and lazycreate:
        logfile = open(logfile,'a')
        print("New mapping", file=logfile)
    elif logfile:
        logfile = open(logfile,'w')
    else:
        logfile = sys.stderr

    unmaps = defaultdict(int)
    singlemaps = defaultdict(int)
    multimaps = defaultdict(int)
    totalreads = defaultdict(int)
    
    if not os.path.isfile(bowtiedb+".fa"):
        print("No bowtie2 database "+bowtiedb, file=sys.stderr)
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
        print("Bam files "+",".join(badsamples)+" does not match fq files", file=sys.stderr)
        print("Aborting", file=sys.stderr)
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
            print(tempfile.gettempdir() +"/"+ currfile + " temp bam files exists", file=sys.stderr)
        print("these files must be deleted to proceed", file=sys.stderr)
        sys.exit(1)
    if len(missingfqfiles) > 0:
        print(",".join(missingfqfiles) + " fastq files missing", file=sys.stderr)
        sys.exit(1)
    mapresults = dict()    
    for samplename in samples:
        bamfile = workingdir+samplename
        
        if lazycreate and os.path.isfile(bamfile+".bam"):
            pass
                
            print("Skipping "+samplename, file=sys.stderr)
            
        else:
            if os.path.isfile(bamfile+".bam"):
                pass

            print("Mapping "+samplename, file=logfile)
            print("Mapping "+samplename, file=sys.stderr)
            logfile.flush()
            mapresults[samplename] = wrapbowtie2(bowtiedb, sampledata.getfastq(samplename),bamfile,scriptdir, trnafile,  logfile=logfile, expname = samplefilename, samplename = samplename)


    if lazycreate:
        #here is where I might add stuff to read old files in lazy mode
        pass
    if mapfile is not None and not lazycreate:
        mapinfo = open(mapfile,'w')                
        print("\t".join(samples), file=mapinfo)
        print("unmap\t"+"\t".join(mapresults[currsample].unmaps for currsample in samples), file=mapinfo)
        print("single\t"+"\t".join(mapresults[currsample].singlemaps for currsample in samples), file=mapinfo)
        print("multi\t"+"\t".join(mapresults[currsample].multimaps for currsample in samples), file=mapinfo)
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
    parser.add_argument('--trnamapfile',
                       help='output table with trna mapping stats')
    parser.add_argument('--bowtiedb',
                       help='Location of Bowtie 2 database')
    parser.add_argument('--lazy', action="store_true", default=False,
                       help='do not remap if mapping results exist')
    parser.add_argument('--local', action="store_true", default=False,
                       help='use local mapping')
    
    args = parser.parse_args()
    main(samplefile = args.samplefile, trnafile= args.trnafile, logfile = args.logfile, bowtiedb = args.bowtiedb, lazy = args.lazy, mapfile = args.mapfile)


