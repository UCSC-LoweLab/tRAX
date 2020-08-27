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
import tempfile
import os

from multiprocessing import Process, Queue, Pool
import time

from distutils.spawn import find_executable



def get_location(program, allowfail = False):
    progloc = find_executable(program)
    if find_executable(program) is None and not allowfail:
        print >>sys.stderr, "Could not find "+program+" in path"
        print >>sys.stderr, "Aborting"
        sys.exit(1)
    else:
        return progloc
        
        
def convertbam(dbname,inputbam, outputbam, scriptdir, force = False, logfile = sys.stderr):
    if not os.path.isfile(outputbam) or force:
   
        #print >>sys.stderr, tempfile.gettempprefix()
        tempprefix = inputbam.split(".")[0]+"_"+str(os.getpid())
        bamconvertcommand =scriptdir+'/convertbam.py '+inputbam+' '+dbname+' | samtools sort -T '+tempfile.gettempdir()+'/convert'+tempprefix+' - -o '+outputbam
        #print >>sys.stderr, bamconvertcommand
        #sys.exit(1)
        if logfile:
            #print >>logfile,  bamconvertcommand
            pass
        bamconvertrun = None
        logfile.flush()
        bamconvertrun = subprocess.Popen(bamconvertcommand, shell = True, stderr = subprocess.PIPE)
        
        output = bamconvertrun.communicate()
        errinfo = output[1]
        if logfile is not None:
            pass
            print >>logfile, errinfo 
        logfile.flush()
        if bamconvertrun.returncode:
            print >>sys.stderr, "Failure to convert bam to genome space"
            print >>sys.stderr, "check logfile"
            logfile.close()
            sys.exit(1)

def samtoolsmerge(bamfiles, outbam, force = False):
    samtoolsloc = get_location("samtools")
    if not os.path.isfile(outbam) or force:
        samcommand = [samtoolsloc,"merge","-f", outbam]
        samcommand.extend(bamfiles)
        
        samtoolsjob = subprocess.Popen(samcommand,stdout = subprocess.PIPE,stderr = subprocess.STDOUT )
        print >>sys.stderr, " ".join(samcommand)
        samtoolsresults = samtoolsjob.communicate()[0]
        print >>sys.stderr, samtoolsresults
    

	
	
def createmultiwigtrackdb(sampledata, expname,trackfile, shortlabel = "", longlabel = "",suffix = '', startpriority = 3.0, stacked = False):
    trackcolors = list(['0,217,47','47,142,248','220,21,235','264,115,6','95,238,230'])
    #trackcolors = list(['55,128,128','204,0,0','0,204,0'])    #'120,235,204'
    currpriority = startpriority
    trackdb = trackfile
    
    print  >>trackdb, "track "+expname+suffix+"tracks            "
    #print  >>trackdb, "compositeTrack on                     "
    print  >>trackdb, "superTrack on show"
    
    print  >>trackdb, "shortLabel "+expname+" "+shortlabel
    print  >>trackdb, "longLabel Data from "+expname+" "+longlabel
    print  >>trackdb, "visibility full"
    #print  >>trackdb, "type bigWig                           "
    #print  >>trackdb, "dragAndDrop on                        "
    #print  >>trackdb, "autoScale on                          "
    #print  >>trackdb, "alwaysZero on                         "
    #print  >>trackdb, "maxHeightPixels 256:100:32              "
    print  >>trackdb, "\n"
    
    #print >>sys.stderr, sampledata.samplelist
    for currrep in sampledata.allreplicates():
        #print >>sys.stderr, sampledata.getrepsamples(currrep)
        for currstrand in ['Plus','Minus']:
            print  >>trackdb, "\ttrack "+currrep+suffix+'_'+currstrand+"tracks"
            print  >>trackdb, "\tcontainer multiWig"
            print  >>trackdb, "\tshortLabel "+currrep+suffix+" "+currstrand+" Strand"
            print  >>trackdb, "\tlongLabel Data from "+expname+" "+currrep+suffix+" "+currstrand+" Strand"
            print  >>trackdb, "\ttype bigWig"
            print  >>trackdb, "\tparent "+expname+suffix+"tracks on"
            print  >>trackdb, "\tdragAndDrop on"
            if stacked:
                print  >>trackdb, "\taggregate solidOverlay"
                
            else:
                print  >>trackdb, "\taggregate transparentOverlay"
            print  >>trackdb, "\tshowSubtrackColorOnUi on"
            print  >>trackdb, "\tautoScale on"
            print  >>trackdb, "\talwaysZero on"
            print  >>trackdb, "\tpriority "+str(currpriority + .1)+"  "
            print  >>trackdb, "\tmaxHeightPixels 256:100:32"
            print  >>trackdb, "\tvisibility full"
            print  >>trackdb, "\n"
            currpriority += .2
            repsamples = sampledata.getrepsamples(currrep)
            for i, currsample in enumerate(repsamples):
                print  >>trackdb, "\t\ttrack "+currsample+suffix+'_'+currstrand+"track"
                print  >>trackdb, "\t\ttype bigWig"
                print  >>trackdb, "\t\tparent "+currrep+suffix+'_'+currstrand+"tracks"
                print  >>trackdb, "\t\tshortLabel "+currsample+suffix+" "+currstrand+" Strand"
                print  >>trackdb, "\t\tlongLabel Data from "+expname+" "+currsample+suffix+" "+currstrand+" Strand"
                print  >>trackdb, "\t\tcolor "+trackcolors[i % len(trackcolors)]+""
                print  >>trackdb, "\t\tbigDataUrl "+currsample+suffix+"."+currstrand+".bw"
                print  >>trackdb, "\t\tvisibility full"
                
                print  >>trackdb, "\n"


    
def createtrackdb(allreps, expname):

    trackdb = open (expname+"/trackhub/trackdb.txt", "w")
    currpriority = 2.3
    
    print  >>trackdb, "track "+expname+"tracks            "
    print  >>trackdb, "compositeTrack on                     "
    print  >>trackdb, "shortLabel Data from "+expname+"   "
    print  >>trackdb, "longLabel Data from "+expname+"    "
    print  >>trackdb, "type bigWig                           "
    print  >>trackdb, "dragAndDrop on                        "
    print  >>trackdb, "autoScale on                          "
    print  >>trackdb, "alwaysZero on                         "
    print  >>trackdb, "maxHeightPixels 100:32:8              "
    print  >>trackdb, "\n"
    
    for currrep in allreps:
        

        print >>trackdb, "track "+currrep+"plus                                             "
        print  >>trackdb, "parent "+expname+"tracks                                   "
        print  >>trackdb, "bigDataUrl "+currrep+".Plus.bw                                 "
        print  >>trackdb, "shortLabel Plus "+currrep+"                                    "
        print  >>trackdb, "longLabel Plus strand coverage "+currrep+" all mapped reads    "
        print  >>trackdb, "color 220,148,44                                             "
        print  >>trackdb, "type bigWig                                                  "
        print  >>trackdb, "priority "+str(currpriority)+"                                       "
        print  >>trackdb, "\n"
        
        
        print  >>trackdb, "track "+currrep+"Minus                                             "
        print  >>trackdb, "parent "+expname+"tracks                                    "
        print  >>trackdb, "bigDataUrl "+currrep+".Minus.bw                                 "
        print  >>trackdb, "shortLabel Minus "+currrep+"                                    "
        print  >>trackdb, "longLabel Minus strand coverage "+currrep+" all mapped reads    "
        print  >>trackdb, "color 112,73,18                                       "
        print  >>trackdb, "type bigWig                                                  "
        print  >>trackdb, "priority "+str(currpriority + .1)+"                                  "
        print  >>trackdb, "\n\n\n"

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
 
    plusjob = subprocess.Popen('zsh -c "bedGraphToBigWig =(samtools view -b -F 0x10 ' +bamfile+' | '+filtercommand+' genomeCoverageBed -scale '+str(1./scalefactor)+' -bg -ibam stdin -g '+faifile+' | sort -k1,1 -k2,2n) '+faifile+' '+directory+"/"+repname+suffix+'.Plus.bw"', shell = True)
    minusjob = subprocess.Popen('zsh -c "bedGraphToBigWig =(samtools view -b -f 0x10 '+bamfile+' | '+filtercommand+' genomeCoverageBed -scale '+str(1./scalefactor)+' -bg -ibam stdin -g '+faifile+' | sort -k1,1 -k2,2n) '+faifile+' '+directory+"/"+repname+suffix+'.Minus.bw"', shell = True)
    plusjob.wait()
    minusjob.wait()
    if plusjob.returncode != 0 or minusjob.returncode != 0:
        print >>sys.stderr, "conversion to bigwig failed"
        sys.exit(1)
    pass




def maketracks(dbname, currbam, genomebam, currsample, scriptdir,trackdir, sizefactors = None):
    if sizefactors is None:
        sizefactors = defaultdict(int)
    convertbam(dbname, currbam, genomebam, scriptdir, force = False)
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
    faidxjob = subprocess.Popen("samtools faidx "+dbname+"-tRNAgenome.fa",shell = True)
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
            print >>sys.stderr, currresult+":"+str(time.time() - starttime)
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
