#!/usr/bin/env python3

import pysam
import sys
import argparse
import string
import itertools
from collections import defaultdict
import os.path
from trnasequtils import *
import subprocess
import tempfile
import os

from multiprocessing import Process, Queue, Pool
import time

from distutils.spawn import find_executable



def get_location(program, allowfail = False):
    progloc = find_executable(program)
    if find_executable(program) is None and not allowfail:
        print("Could not find "+program+" in path", file=sys.stderr)
        print("Aborting", file=sys.stderr)
        sys.exit(1)
    else:
        return progloc
        
        
def convertbam(dbname,inputbam, outputbam, scriptdir, force = False, logfile = sys.stderr):
    if not os.path.isfile(outputbam) or force:
   
        #print >>sys.stderr, tempfile.gettempprefix()
        tempprefix = inputbam.split(".")[0]+"_"+str(os.getpid())
        
        bamconvertcommand =scriptdir+'/convertbam.py '+inputbam+' '+dbname+' | samtools sort -T '+tempfile.gettempdir()+'/convert'+tempprefix+' - -o '+outputbam
        #print >>sys.stderr, bamconvertcommand
        print(bamconvertcommand, file=sys.stderr)
        #sys.exit(1)
        if logfile:
            #print >>logfile,  bamconvertcommand
            pass
        bamconvertrun = None
        logfile.flush()
        bamconvertrun = subprocess.Popen(bamconvertcommand, shell = True, stderr = subprocess.PIPE,universal_newlines=True )
        
        output = bamconvertrun.communicate()
        errinfo = output[1]
        if logfile is not None:
            pass
            print(errinfo, file=logfile) 
        logfile.flush()
        if bamconvertrun.returncode:
            print("Failure to convert bam to genome space", file=sys.stderr)
            print("check logfile", file=sys.stderr)
            logfile.close()
            sys.exit(1)

def samtoolsmerge(bamfiles, outbam, force = False):
    samtoolsloc = get_location("samtools")
    if not os.path.isfile(outbam) or force:
        samcommand = [samtoolsloc,"merge","-f", outbam]
        samcommand.extend(bamfiles)
        
        samtoolsjob = subprocess.Popen(samcommand,stdout = subprocess.PIPE,stderr = subprocess.STDOUT,universal_newlines=True  )
        print(" ".join(samcommand), file=sys.stderr)
        samtoolsresults = samtoolsjob.communicate()[0]
        print(samtoolsresults, file=sys.stderr)
    

	
	
def createmultiwigtrackdb(sampledata, expname,trackfile, shortlabel = "", longlabel = "",suffix = '', startpriority = 3.0, stacked = False):
    trackcolors = list(['0,217,47','47,142,248','220,21,235','264,115,6','95,238,230'])
    #trackcolors = list(['55,128,128','204,0,0','0,204,0'])    #'120,235,204'
    currpriority = startpriority
    trackdb = trackfile
    
    print("track "+expname+suffix+"tracks            ", file=trackdb)
    #print  >>trackdb, "compositeTrack on                     "
    print("superTrack on show", file=trackdb)
    
    print("shortLabel "+expname+" "+shortlabel, file=trackdb)
    print("longLabel Data from "+expname+" "+longlabel, file=trackdb)
    print("visibility full", file=trackdb)
    #print  >>trackdb, "type bigWig                           "
    #print  >>trackdb, "dragAndDrop on                        "
    #print  >>trackdb, "autoScale on                          "
    #print  >>trackdb, "alwaysZero on                         "
    #print  >>trackdb, "maxHeightPixels 256:100:32              "
    print("\n", file=trackdb)
    
    #print >>sys.stderr, sampledata.samplelist
    for currrep in sampledata.allreplicates():
        #print >>sys.stderr, sampledata.getrepsamples(currrep)
        for currstrand in ['Plus','Minus']:
            print("\ttrack "+currrep+suffix+'_'+currstrand+"tracks", file=trackdb)
            print("\tcontainer multiWig", file=trackdb)
            print("\tshortLabel "+currrep+suffix+" "+currstrand+" Strand", file=trackdb)
            print("\tlongLabel Data from "+expname+" "+currrep+suffix+" "+currstrand+" Strand", file=trackdb)
            print("\ttype bigWig", file=trackdb)
            print("\tparent "+expname+suffix+"tracks on", file=trackdb)
            print("\tdragAndDrop on", file=trackdb)
            if stacked:
                print("\taggregate solidOverlay", file=trackdb)
                
            else:
                print("\taggregate transparentOverlay", file=trackdb)
            print("\tshowSubtrackColorOnUi on", file=trackdb)
            print("\tautoScale on", file=trackdb)
            print("\talwaysZero on", file=trackdb)
            print("\tpriority "+str(currpriority + .1)+"  ", file=trackdb)
            print("\tmaxHeightPixels 256:100:32", file=trackdb)
            print("\tvisibility full", file=trackdb)
            print("\n", file=trackdb)
            currpriority += .2
            repsamples = sampledata.getrepsamples(currrep)
            for i, currsample in enumerate(repsamples):
                print("\t\ttrack "+currsample+suffix+'_'+currstrand+"track", file=trackdb)
                print("\t\ttype bigWig", file=trackdb)
                print("\t\tparent "+currrep+suffix+'_'+currstrand+"tracks", file=trackdb)
                print("\t\tshortLabel "+currsample+suffix+" "+currstrand+" Strand", file=trackdb)
                print("\t\tlongLabel Data from "+expname+" "+currsample+suffix+" "+currstrand+" Strand", file=trackdb)
                print("\t\tcolor "+trackcolors[i % len(trackcolors)]+"", file=trackdb)
                print("\t\tbigDataUrl "+currsample+suffix+"."+currstrand+".bw", file=trackdb)
                print("\t\tvisibility full", file=trackdb)
                
                print("\n", file=trackdb)


    
def createtrackdb(allreps, expname):

    trackdb = open (expname+"/trackhub/trackdb.txt", "w")
    currpriority = 2.3
    
    print("track "+expname+"tracks            ", file=trackdb)
    print("compositeTrack on                     ", file=trackdb)
    print("shortLabel Data from "+expname+"   ", file=trackdb)
    print("longLabel Data from "+expname+"    ", file=trackdb)
    print("type bigWig                           ", file=trackdb)
    print("dragAndDrop on                        ", file=trackdb)
    print("autoScale on                          ", file=trackdb)
    print("alwaysZero on                         ", file=trackdb)
    print("maxHeightPixels 100:32:8              ", file=trackdb)
    print("\n", file=trackdb)
    
    for currrep in allreps:
        

        print("track "+currrep+"plus                                             ", file=trackdb)
        print("parent "+expname+"tracks                                   ", file=trackdb)
        print("bigDataUrl "+currrep+".Plus.bw                                 ", file=trackdb)
        print("shortLabel Plus "+currrep+"                                    ", file=trackdb)
        print("longLabel Plus strand coverage "+currrep+" all mapped reads    ", file=trackdb)
        print("color 220,148,44                                             ", file=trackdb)
        print("type bigWig                                                  ", file=trackdb)
        print("priority "+str(currpriority)+"                                       ", file=trackdb)
        print("\n", file=trackdb)
        
        
        print("track "+currrep+"Minus                                             ", file=trackdb)
        print("parent "+expname+"tracks                                    ", file=trackdb)
        print("bigDataUrl "+currrep+".Minus.bw                                 ", file=trackdb)
        print("shortLabel Minus "+currrep+"                                    ", file=trackdb)
        print("longLabel Minus strand coverage "+currrep+" all mapped reads    ", file=trackdb)
        print("color 112,73,18                                       ", file=trackdb)
        print("type bigWig                                                  ", file=trackdb)
        print("priority "+str(currpriority + .1)+"                                  ", file=trackdb)
        print("\n\n\n", file=trackdb)

        currpriority += .2
        
'''
chr9
  37244 chr1_KI270762v1_alt     141679  141702  1.1772
  37245 chr1_KI270766v1_alt     92568   92634   1.1772
  sort -k1,1 -k2,2n" with LC_COLLATE=C
'''

def makebigwigs(bamfile, repname, faifile, directory,scriptdir,filterloci = False, suffix = '',scalefactor = 1):
    #print >>sys.stderr, 'zsh -c "bedGraphToBigWig =(samtools view -b -F 0x10 '+bamfile+' | genomeCoverageBed -bg -ibam stdin -g '+faifile+') '+faifile+' '+directory+"/"+repname+'.Plus.bw"'
    
    #print >>sys.stderr, 'zsh -c "bedGraphToBigWig =(samtools view -b -F 0x10 '+bamfile+' | genomeCoverageBed -scale '+' -bg -ibam stdin -g '+faifile+') ' +faifile+' '+directory+"/"+repname+'.Plus.bw"'
    filtercommand = ''
    if filterloci:
        filtercommand = scriptdir+'/filterunique.py --uniqloci | '
    #print >>sys.stderr, 'zsh -c "bedGraphToBigWig =(samtools view -b -F 0x10 ' +bamfile+' | '+filtercommand+' genomeCoverageBed -scale '+str(1./scalefactor)+' -bg -ibam stdin -g '+faifile+' | sort -k1,1 -k2,2n) '+faifile+' '+directory+"/"+repname+suffix+'.Plus.bw"'
    print(faifile, file=sys.stderr )
    print('zsh -c "bedGraphToBigWig =(samtools view -b -F 0x10 ' +bamfile+' | '+filtercommand+' genomeCoverageBed -scale '+str(1./scalefactor)+' -bg -ibam stdin -g '+faifile+' | sort -k1,1 -k2,2n) '+faifile+' '+directory+"/"+repname+suffix+'.Plus.bw"', file=sys.stderr)
    plusjob = subprocess.Popen('zsh -c "bedGraphToBigWig =(samtools view -b -F 0x10 ' +bamfile+' | '+filtercommand+' genomeCoverageBed -scale '+str(1./scalefactor)+' -bg -ibam stdin -g '+faifile+' | sort -k1,1 -k2,2n) '+faifile+' '+directory+"/"+repname+suffix+'.Plus.bw"', shell = True,universal_newlines=True)
    minusjob = subprocess.Popen('zsh -c "bedGraphToBigWig =(samtools view -b -f 0x10 '+bamfile+' | '+filtercommand+' genomeCoverageBed -scale '+str(1./scalefactor)+' -bg -ibam stdin -g '+faifile+' | sort -k1,1 -k2,2n) '+faifile+' '+directory+"/"+repname+suffix+'.Minus.bw"', shell = True,universal_newlines=True)
    plusjob.wait()
    minusjob.wait()
    if plusjob.returncode != 0 or minusjob.returncode != 0:
        print("conversion to bigwig of "+bamfile+" failed", file=sys.stderr)
        sys.exit(1)
    pass




def maketracks(dbname, currbam, genomebam, currsample, scriptdir,trackdir, sizefactors = None):
    if sizefactors is None:
        sizefactors = defaultdict(int)
    convertbam(dbname, currbam, genomebam, scriptdir, force = True)
    makebigwigs(genomebam, currsample, dbname+"-tRNAgenome.fa.fai",trackdir,scriptdir, scalefactor =  sizefactors[currsample])
    #makebigwigs(genomebam, currsample, dbname+"-tRNAgenome.fa.fai",trackdir,scriptdir, filterloci = True, suffix = 'uniqloci', scalefactor =  sizefactors[currsample])	
    return genomebam
    
def maketracksspool(args):
    return maketracks(*args[0], **args[1])
def compressargs( *args, **kwargs):
    return tuple([args, kwargs])
def main(**args):
    dbname = args["genomedatabase"]
    samplefilename = args["samplefile"]
    sampledata = samplefile(args["samplefile"])
    expname = args["expname"]
    trackdir = expname+"/trackhub"
    scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"
    sizefactors = getsizefactors( expname+"/"+expname+"-SizeFactors.txt")
    if not os.path.exists(trackdir):
        os.makedirs(trackdir)
    allsamples = sampledata.getsamples()
    faidxjob = subprocess.Popen("samtools faidx "+dbname+"-tRNAgenome.fa",shell = True,universal_newlines=True )
    faidxjob.wait()
    convetbampool = Pool(processes=8)
    trackargs = list()
    starttime = time.time()
    #print >>sys.stderr, starttime
    #print >>sys.stderr, "**||"
    threadmode = True
    for currsample in allsamples:
        currbam = sampledata.getbam(currsample)
        genomebam = currsample+"-genome.bam"
        if not threadmode:
            maketracks(dbname, currbam, genomebam, currsample,scriptdir,trackdir, sizefactors = sizefactors)
        else:
            trackargs.append(compressargs(dbname, currbam, genomebam, currsample,scriptdir,trackdir, sizefactors = sizefactors))
            
            #convertbam(dbname, currbam, genomebam, scriptdir, force = True)
        #makebigwigs(genomebam, currsample, dbname+"-tRNAgenome.fa.fai",trackdir,scriptdir, scalefactor =  sizefactors[currsample])
        #makebigwigs(genomebam, currsample, dbname+"-tRNAgenome.fa.fai",trackdir,scriptdir, filterloci = True, suffix = 'uniqloci', scalefactor =  sizefactors[currsample])
    if threadmode:
        for  currresult in convetbampool.imap_unordered(maketracksspool, trackargs):
            print(currresult+":"+str(time.time() - starttime), file=sys.stderr)
            pass
        

    trackfile = open(expname+"/trackhub/trackdb.txt", "w")
    createmultiwigtrackdb(sampledata,expname, trackfile, shortlabel = "all", longlabel = "all")
    #createmultiwigtrackdb(sampledata,expname,trackfile, shortlabel = "unique mapping only", longlabel = "uniquely mapping only", suffix = 'uniqloci', startpriority = 8.0)
    '''


    for currrep in sampledata.allreplicates():
        repsamples = sampledata.getrepsamples(currrep)
        samtoolsmerge(list(curr+"-genome.bam" for curr in repsamples), currrep+"-mergegenome.bam", True)
        pysam.index(currrep+"-mergegenome.bam")
        makebigwigs(currrep+"-mergegenome.bam", currrep, dbname+"-tRNAgenome.fa.fai",trackdir)
    
    createtrackdb(sampledata.allreplicates(),expname)
    '''
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert TRAX bam file into genome bam file.')
    parser.add_argument('--genomedatabase', 
                       help='fasta sequence of genome')
    parser.add_argument('--samplefile', 
                       help='sample file')
    parser.add_argument('--expname', 
                       help='experiment name')

    
    args = vars(parser.parse_args())
    main(args)    
