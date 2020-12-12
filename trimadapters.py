#!/usr/bin/env python

import pysam
import sys
import argparse
import string
import itertools
from collections import defaultdict
import os.path
import re
import subprocess
from multiprocessing import cpu_count
from multiprocessing.dummy import Pool


def updatemanifest(indexfilename, runname, runfile):
    filelocs = dict()
    try:
        indexfile = open(indexfilename)
        for currline in indexfile:
            fields = currline.split()
            if len(fields) > 1:
                filelocs[fields[0]] = fields[1]
        indexfile.close()
    except IOError as e:
        pass
    filelocs = filelocs
    filelocs[runname] = runfile
    outfile = open(indexfilename, "w") 
    for currname in filelocs.iterkeys():
           print >>outfile, currname+"\t"+filelocs[currname]
           

def runrscript(*script):
    print >>sys.stderr, "Rscript "+" ".join(script)

    retcode = subprocess.call("Rscript "+" ".join(script), shell=True,  stderr = subprocess.STDOUT)

    if retcode > 0:
        print >>sys.stderr, "R script "+script[0]+" failed"

        
        #sys.exit()
    return retcode
    
def subprocesspool(argnames):
    samplename = argnames[0]
    args = argnames[1]
    process = subprocess.Popen(*args[0], **args[1])
    process.wait()
    return samplename, process
def compressargs( *args, **kwargs):
    return tuple([args, kwargs])
    

def readseqprep(processoutput):
    seqprepcounts = dict()
    output = processoutput.communicate()
    errinfo = output[1]
    if processoutput.returncode != 0:
        print >>sys.stderr, "seqprep failed"
        print >>sys.stderr, errinfo
    
    #errinfo
    for line in errinfo.split("\n"):
    
        
        totalmatch = rereadtotal.match(line)
        mergematch = rereadmerge.match(line)
        discardmatch = rereaddiscard.match(line)
        if totalmatch:
            totalreads = int(totalmatch.group(1).replace(",",""))
        elif mergematch: 
            merged = int(mergematch.group(1).replace(",",""))
        elif discardmatch:
            discard = int(discardmatch.group(1).replace(",",""))
    
    seqprepcounts["merged"] = merged
    seqprepcounts["unmerged"] = totalreads - (merged + discard)
    seqprepcounts["discarded"] = discard
    return seqprepcounts, errinfo
    
def readcutadapt(processoutput):
    cutadaptcounts = dict()
    output = processoutput.communicate()
    errinfo = output[1]
    if processoutput.returncode != 0:
        print >>sys.stderr, "cutadapt failed"
        print >>sys.stderr, errinfo
    #print >>sys.stderr, errinfo
    
    for line in errinfo.split("\n"):
        totalmatch = recutadapttotal.match(line)
        trimmatch = recutadapttrimmed.match(line) 
        discardmatch = recutadaptshort.match(line) 
        writtenmatch = recutadaptwritten.match(line) 
        #trimmed = None
    
        if totalmatch:
            totalreads = int(totalmatch.group(1).replace(",",""))
        elif trimmatch: 
            trimmed = int(trimmatch.group(1).replace(",",""))
        elif discardmatch:
            discard = int(discardmatch.group(1).replace(",",""))
        elif writtenmatch:
            written = int(writtenmatch.group(1).replace(",",""))
           
    cutadaptcounts["trimmed"] = trimmed
    cutadaptcounts["untrimmed"] = totalreads - (trimmed + discard)
    cutadaptcounts["discarded"] = discard
    return cutadaptcounts, errinfo
scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"


parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
parser.add_argument('--runname',required=True,
                   help='run name to be used')
parser.add_argument('--runfile',required=True,
                   help='run file')
parser.add_argument('--firadapter',default ='AGATCGGAAGAGCACACGTC' ,
                   help='subset file')
parser.add_argument('--secadapter',default ='GATCGTCGGACTGTAGAACTC' ,
                   help='subset file')
parser.add_argument('--minlength',default ='15' ,
                   help='minimum length of sequence')
parser.add_argument('--singleend', action="store_true", default=False,
                   help='single-end mode (uses cutadapt)')
parser.add_argument('--cores',
                   help='number of processors to use')

args = parser.parse_args()

runname = args.runname
seqprepfile = args.runfile
firadapter = args.firadapter
secadapter =   args.secadapter
minlength = args.minlength
singleendmode = args.singleend
cores = cpu_count()
if args.cores is not None:
    cores = args.cores

#firadapter = 'AGATCGGAAGAGCACACGTC' 
#secadapter = 'GATCGTCGGACTGTAGAACTC'  


'''
Pairs Processed:	1207400
Pairs Merged:	1171370
Pairs With Adapters:	1190959
Pairs Discarded:	23902
CPU Time Used (Minutes):	2.748833

Total reads processed:                  10,000
Reads with adapters:                     9,926 (99.3%)
Reads that were too short:                  19 (0.2%)
Reads written (passing filters):         9,981 (99.8%)

'''     


rereadtotal = re.compile(r'Pairs Processed:\s+(\d+)')
rereadmerge = re.compile(r'Pairs Merged:\s+(\d+)' )
rereadapter = re.compile(r'Pairs With Adapters:\s+(\d+)')
rereaddiscard = re.compile(r'Pairs Discarded:\s+(\d+)' )     

recutadapttotal = re.compile(r'Total reads processed:\s+([\d\,]+)')
recutadapttrimmed = re.compile(r'Reads with adapters:\s+([\d\,]+)' )
recutadaptshort = re.compile(r'Reads that were too short:\s+([\d\,]+)')
recutadaptwritten = re.compile(r'Reads written (passing filters):\s+([\d\,]+)' )     


samplefiles = dict()
sampleorder = list()
cutadaptorder = list()
for currline in open(seqprepfile):
    
    fields = currline.split()
    if not singleendmode and len(fields) > 2:
        samplefiles[fields[0]] = tuple([fields[1],fields[2]])
        sampleorder.append(fields[0])
    elif singleendmode and len(fields) > 1:
        samplefiles[fields[0]] = tuple([fields[1]])
        cutadaptorder.append(fields[0])
        sampleorder.append(fields[0])


if len(sampleorder) < 1 and len(cutadaptorder) < 1:
    print >>sys.stderr, "Failed to read "+seqprepfile
    print >>sys.stderr, "Perhaps you failed to specify --singleend mode?"
    sys.exit(1)

    
totalreads = None
merged = None
discard = None
samplenum = 1
seqprepcounts = dict()
cutadaptcounts = dict()
allsamples = set()

minsize = 15
outputfiles = list()

prepout = ""
for currsample in samplefiles.iterkeys():
    for currfile in samplefiles[currsample]:
        if not os.path.isfile(currfile):
            print >>sys.stderr, currfile +" does not exist"
            sys.exit(1)
seqprepruns = dict()
cutadaptruns = dict()


#print >>sys.stderr, cores

trimpool = Pool(processes=int(cores))


for currsample in sampleorder:
    
    if not singleendmode:
        #seqprepcommmand = program+' -x '+bowtiedb+' -k '+str(maxmaps)+' --very-sensitive --ignore-quals --np 5 --reorder -p '+str(numcores)+' -U '+unpaired
        seqprepcommmand = 'SeqPrep -L '+str(minsize)+ ' -A '+firadapter+' -B '+secadapter +' -f '+samplefiles[currsample][0]+'  -r '+samplefiles[currsample][1]+' -1 '+currsample+'_left.fastq.gz     -2 '+currsample+'_right.fastq.gz   -s '+currsample+'_merge.fastq.gz'  
        outputfiles.append(currsample+'_merge.fastq.gz')
        print >>sys.stderr, seqprepcommmand
        #bowtiecommand = bowtiecommand + ' | '+scriptdir+'choosemappings.py '+trnafile+' | samtools sort - '+outfile
        
        seqprepruns[currsample] = None
        seqprepruns[currsample] = compressargs(seqprepcommmand, shell = True, stderr = subprocess.PIPE)

    else:
        #seqprepcommmand = program+' -x '+bowtiedb+' -k '+str(maxmaps)+' --very-sensitive --ignore-quals --np 5 --reorder -p '+str(numcores)+' -U '+unpaired
        #cutadapt -m 15 --adapter='TGGAATTCTCGGGTGCCAAGG'  Testicular_sperm/Fraction4/Fraction4_S2_L002_R1_001.fastq.gz                           | gzip -c > trimmed/Fraction4_S2_L002_R1_001_TRIM.fastq.gz                      

        cutadaptcommand = 'cutadapt -m '+str(minsize)+ ' --adapter='+firadapter+' '+samplefiles[currsample][0]  +' | gzip -c >'+ currsample+'_trimmed.fastq.gz'
        outputfiles.append(currsample+'_trimmed.fastq.gz')
        #print >>sys.stderr, cutadaptcommand
        #bowtiecommand = bowtiecommand + ' | '+scriptdir+'choosemappings.py '+trnafile+' | samtools sort - '+outfile
        
        cutadaptruns[currsample] = None
        cutadaptruns[currsample] = compressargs(cutadaptcommand, shell = True, stderr = subprocess.PIPE)
        


if not singleendmode:
    results = trimpool.imap_unordered(subprocesspool, list(tuple([currsample, seqprepruns[currsample]]) for currsample in sampleorder))
    for samplename, spoutput in results:
        print >>sys.stderr, samplename +" merged"
        #print >>sys.stderr, spoutput
        seqprepcounts[samplename], errinfo = readseqprep(spoutput)
        prepout += errinfo
else:
    results = trimpool.imap_unordered(subprocesspool, list(tuple([currsample, cutadaptruns[currsample]]) for currsample in sampleorder))
    for samplename, caoutput in results:
        
        print >>sys.stderr, samplename +" trimmed"
        cutadaptcounts[samplename], errinfo = readcutadapt(caoutput)
        prepout += errinfo
        
#print >>sys.stderr,  cutadaptcounts.keys()
#sys.exit()
''' 
for currsample in sampleorder:
    if not singleendmode:
        output = seqprepruns[currsample].communicate()
        errinfo = output[1]
        if seqprepruns[currsample].returncode != 0:
            print >>sys.stderr, "seqprep failed"
            print >>sys.stderr, errinfo
        
        prepout += errinfo
        for line in errinfo.split("\n"):
        
            
            totalmatch = rereadtotal.match(line)
            mergematch = rereadmerge.match(line)
            discardmatch = rereaddiscard.match(line)
            if totalmatch:
                totalreads = int(totalmatch.group(1).replace(",",""))
            elif mergematch: 
                merged = int(mergematch.group(1).replace(",",""))
            elif discardmatch:
                discard = int(discardmatch.group(1).replace(",",""))
        
        seqprepcounts[currsample]["merged"] = merged
        seqprepcounts[currsample]["unmerged"] = totalreads - (merged + discard)
        seqprepcounts[currsample]["discarded"] = discard
    else:
        output = cutadaptruns[currsample].communicate()
        errinfo = output[1]
        if cutadaptruns[currsample].returncode != 0:
            print >>sys.stderr, "seqprep failed"
            print >>sys.stderr, errinfo
        print >>sys.stderr, errinfo
        prepout += errinfo
        for line in errinfo.split("\n"):
            totalmatch = recutadapttotal.match(line)
            trimmatch = recutadapttrimmed.match(line) 
            discardmatch = recutadaptshort.match(line) 
            writtenmatch = recutadaptwritten.match(line) 
            #trimmed = None

            if totalmatch:
                totalreads = int(totalmatch.group(1).replace(",",""))
            elif trimmatch: 
                trimmed = int(trimmatch.group(1).replace(",",""))
            elif discardmatch:
                discard = int(discardmatch.group(1).replace(",",""))
            elif writtenmatch:
                written = int(writtenmatch.group(1).replace(",",""))
               
        cutadaptcounts[currsample]["trimmed"] = trimmed
        cutadaptcounts[currsample]["untrimmed"] = totalreads - (trimmed + discard)
        cutadaptcounts[currsample]["discarded"] = discard
'''
        
#print >>sys.stderr, cutadaptcounts

if not singleendmode:
    samplefile = open(runname+"_sp.txt", "w")
    logfile = open(runname+"_log.txt", "w")
    replicatefile = open(runname+"_manifest.txt", "w")
    for i, curr in enumerate(sampleorder):
        print >>replicatefile, curr+"\t"+outputfiles[i]
    print >>samplefile,"\t".join(sampleorder)
    for currtype in ["merged","unmerged","discarded"]:
        print >>samplefile,currtype+"\t"+"\t".join(str(seqprepcounts[currsample][currtype]) for currsample in sampleorder)
        
    samplefile.close()
    print >>logfile ,"samplefile: "+seqprepfile
    print >>logfile ,"first adapter: "+ firadapter
    print >>logfile ,"second adapter: "+ secadapter
    print >>logfile ,"output files: "+ ",".join(outputfiles)
    print >>logfile ," ".join(sys.argv)
    print >>logfile ,"************************"
    print >>logfile ,prepout


    logfile.close()    
    runrscript(scriptdir+"/featuretypesreal.R",runname+"_sp.txt",runname+"_sp.pdf")
else:
    samplefile = open(runname+"_ca.txt", "w")
    logfile = open(runname+"_log.txt", "w")
    replicatefile = open(runname+"_manifest.txt", "w")
    for i, curr in enumerate(cutadaptorder):
        print >>replicatefile, curr+"\t"+outputfiles[i]
    print >>samplefile,"\t".join(cutadaptorder)
    for currtype in ["trimmed","untrimmed","discarded"]:
        print >>samplefile,currtype+"\t"+"\t".join(str(cutadaptcounts[currsample][currtype]) for currsample in cutadaptorder)
        
    samplefile.close()
    print >>logfile ,"samplefile: "+seqprepfile
    print >>logfile ,"adapter: "+ firadapter
    print >>logfile ,"output files: "+ ",".join(outputfiles)
    print >>logfile ," ".join(sys.argv)
    print >>logfile ,"************************"
    print >>logfile ,prepout

    logfile.close()
    runrscript(scriptdir+"/featuretypesreal.R",runname+"_ca.txt",runname+"_ca.pdf")


updatemanifest("trimindex.txt",runname,seqprepfile)

