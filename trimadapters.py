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
    for currname in filelocs.keys():
           print(currname+"\t"+filelocs[currname], file=outfile)
           

def runrscript(*script):
    print("Rscript "+" ".join(script), file=sys.stderr)

    retcode = subprocess.call("Rscript "+" ".join(script), shell=True,  stderr = subprocess.STDOUT, universal_newlines=True)

    if retcode > 0:
        print("R script "+script[0]+" failed", file=sys.stderr)

        
        #sys.exit()
    return retcode
    
def subprocesspool(argnames):
    samplename = argnames[0]
    args = argnames[1]
    process = subprocess.Popen(*args[0], **args[1], universal_newlines=True)
    process.wait()
    return samplename," ".join(args[0]),  process
def compressargs( *args, **kwargs):
    return tuple([args, kwargs])
    

def readseqprep(processoutput):
    seqprepcounts = dict()
    output = processoutput.communicate()
    errinfo = output[1]
    if processoutput.returncode != 0:
        print("seqprep failed", file=sys.stderr)
        print(errinfo, file=sys.stderr)
    
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
        print("cutadapt failed", file=sys.stderr)
        print(errinfo, file=sys.stderr)
    print(errinfo, file=sys.stderr)
    trimmed = 0 
    discard = 0
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
           
    cutadaptcounts["trimmed"] = trimmed - discard
    cutadaptcounts["untrimmed"] = totalreads - (trimmed)
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
parser.add_argument('--umilength',default ='0',
                   help='length of UMI (uses umi_tools to extract if present)')
parser.add_argument('--umithreeprime', action="store_true", default=False,
                   help='umi is at the three prime end')
parser.add_argument('--cores',
                   help='number of processors to use')

args = parser.parse_args()

runname = args.runname
seqprepfile = args.runfile
firadapter = args.firadapter
secadapter =   args.secadapter
minlength = args.minlength
singleendmode = args.singleend
threeprimeumi = args.umithreeprime
umilength = int(args.umilength)
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
    print("Failed to read "+seqprepfile, file=sys.stderr)
    print("Perhaps you failed to specify --singleend mode?", file=sys.stderr)
    sys.exit(1)

    
totalreads = None
merged = None
discard = None
samplenum = 1
seqprepcounts = dict()
cutadaptcounts = dict()
allsamples = set()

minsize = int(minlength) + int(umilength)
outputfiles = list()

prepout = ""
for currsample in samplefiles.keys():
    for currfile in samplefiles[currsample]:
        if not os.path.isfile(currfile):
            print(currfile +" does not exist", file=sys.stderr)
            sys.exit(1)
seqprepruns = dict()
cutadaptruns = dict()


#print >>sys.stderr, cores

trimpool = Pool(processes=int(cores))
umiargs = ""
if threeprimeumi:
    umiargs += " --3prime "

for currsample in sampleorder:
    
    if not singleendmode:
        seqprepcommmand = ""
        if umilength  > 0:
            seqprepcommmand = 'SeqPrep -L '+str(minsize)+ ' -A '+firadapter+' -B '+secadapter +' -f '+samplefiles[currsample][0]+'  -r '+samplefiles[currsample][1]+' -1 '+currsample+'_left.fastq.gz     -2 '+currsample+'_right.fastq.gz   -s '+currsample+'_m.fastq.gz '
            seqprepcommmand += ";umi_tools extract -v 0 --stdin="+currsample+"_m.fastq.gz "+umiargs+"--bc-pattern="+("N"*int(umilength))+" --stdout="+currsample+'_merge.fastq.gz -L '+currsample+'_umilog.txt'
        else:
            seqprepcommmand = 'SeqPrep -L '+str(minsize)+ ' -A '+firadapter+' -B '+secadapter +' -f '+samplefiles[currsample][0]+'  -r '+samplefiles[currsample][1]+' -1 '+currsample+'_left.fastq.gz     -2 '+currsample+'_right.fastq.gz   -s '+currsample+'_merge.fastq.gz'


        outputfiles.append(currsample+'_merge.fastq.gz')
        print(seqprepcommmand, file=sys.stderr)
        #bowtiecommand = bowtiecommand + ' | '+scriptdir+'choosemappings.py '+trnafile+' | samtools sort - '+outfile
        
        seqprepruns[currsample] = None
        seqprepruns[currsample] = compressargs(seqprepcommmand, shell = True, stderr = subprocess.PIPE, universal_newlines=True)

    else:
        cutadaptcommand = ""
        #seqprepcommmand = program+' -x '+bowtiedb+' -k '+str(maxmaps)+' --very-sensitive --ignore-quals --np 5 --reorder -p '+str(numcores)+' -U '+unpaired
        #cutadapt -m 15 --adapter='TGGAATTCTCGGGTGCCAAGG'  Testicular_sperm/Fraction4/Fraction4_S2_L002_R1_001.fastq.gz                           | gzip -c > trimmed/Fraction4_S2_L002_R1_001_TRIM.fastq.gz                      
        if umilength > 0:
            cutadaptcommand = 'cutadapt -m '+str(minsize)+ ' --adapter='+firadapter+' '+samplefiles[currsample][0]  +' | gzip -c >'+ currsample+'_t.fastq.gz'
            cutadaptcommand += ";umi_tools extract --stdin="+currsample+"_t.fastq.gz "+umiargs+"--bc-pattern="+("N"*int(umilength))+" --stdout="+currsample+'_trimmed.fastq.gz -L '+currsample+'_umilog.txt'
            
        else:
            cutadaptcommand = 'cutadapt -m '+str(minsize)+ ' --adapter='+firadapter+' '+samplefiles[currsample][0]  +' | gzip -c >'+ currsample+'_trimmed.fastq.gz'


        outputfiles.append(currsample+'_trimmed.fastq.gz')
        print(cutadaptcommand, file=sys.stderr)
        #bowtiecommand = bowtiecommand + ' | '+scriptdir+'choosemappings.py '+trnafile+' | samtools sort - '+outfile
        
        cutadaptruns[currsample] = None
        cutadaptruns[currsample] = compressargs(cutadaptcommand, shell = True, stderr = subprocess.PIPE, universal_newlines=True)
        
print("***", file=sys.stderr)
'''

umi_tools extract --bc-pattern=NNNNNN --log=umi.log --stdout=output.fastq.gz 
'''
if not singleendmode:
    results = trimpool.imap_unordered(subprocesspool, list(tuple([currsample, seqprepruns[currsample]]) for currsample in sampleorder))
    print("***||", file=sys.stderr)
    for samplename, command, spoutput in results:
        print(samplename +" merged", file=sys.stderr)
        #print >>sys.stderr, spoutput
        seqprepcounts[samplename], errinfo = readseqprep(spoutput)
        prepout += samplename +"\n"
        prepout += command+"\n"
        prepout += errinfo+"\n"
else:
    results = trimpool.imap_unordered(subprocesspool, list(tuple([currsample, cutadaptruns[currsample]]) for currsample in sampleorder))
    for samplename,  command, caoutput in results:
        
        print(samplename +" trimmed", file=sys.stderr)
        cutadaptcounts[samplename], errinfo = readcutadapt(caoutput)
        prepout += samplename+"\n"
        prepout += command+"\n"
        prepout += errinfo+"\n"
        
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
        print(curr+"\t"+outputfiles[i], file=replicatefile)
    print("\t".join(sampleorder), file=samplefile)
    for currtype in ["merged","unmerged","discarded"]:
        print(currtype+"\t"+"\t".join(str(seqprepcounts[currsample][currtype]) for currsample in sampleorder), file=samplefile)
        
    samplefile.close()
    print("samplefile: "+seqprepfile, file=logfile)
    print("first adapter: "+ firadapter, file=logfile)
    print("second adapter: "+ secadapter, file=logfile)
    print("output files: "+ ",".join(outputfiles), file=logfile)
    print(" ".join(sys.argv), file=logfile)
    print("************************", file=logfile)
    print(prepout, file=logfile)


    logfile.close()    
    runrscript(scriptdir+"/featuretypesreal.R",runname+"_sp.txt",runname+"_sp.pdf")
else:
    samplefile = open(runname+"_ca.txt", "w")
    logfile = open(runname+"_log.txt", "w")
    replicatefile = open(runname+"_manifest.txt", "w")
    for i, curr in enumerate(cutadaptorder):
        print(curr+"\t"+outputfiles[i], file=replicatefile)
    print("\t".join(cutadaptorder), file=samplefile)
    for currtype in ["trimmed","untrimmed","discarded"]:
        print(currtype+"\t"+"\t".join(str(cutadaptcounts[currsample][currtype]) for currsample in cutadaptorder), file=samplefile)
        
    samplefile.close()
    print("samplefile: "+seqprepfile, file=logfile)
    print("adapter: "+ firadapter, file=logfile)
    print("output files: "+ ",".join(outputfiles), file=logfile)
    print(" ".join(sys.argv), file=logfile)
    print("************************", file=logfile)
    print(prepout, file=logfile)

    logfile.close()
    runrscript(scriptdir+"/featuretypesreal.R",runname+"_ca.txt",runname+"_ca.pdf")


updatemanifest("trimindex.txt",runname,seqprepfile)

