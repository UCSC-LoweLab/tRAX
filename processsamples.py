#!/usr/bin/env python

import sys
import argparse
import os.path
from collections import defaultdict
from trnasequtils import *
import itertools
import subprocess
import os
import time

import mapreads
import countreads
import getcoverage
import getends
import countreadtypes
import maketrackhub
import traxqc
from distutils.version import LooseVersion, StrictVersion
from multiprocessing import Pool, cpu_count


 


#expname is experiment name
#dbname is database name
#samplefile is sample file
#$4 is bed feature for other sRNAs




parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
parser.add_argument('--experimentname',required=True,
                   help='experiment name to be used')
parser.add_argument('--databasename',required=True,
                   help='name of the tRNA database')
parser.add_argument('--samplefile',required=True,
                   help='sample file')
parser.add_argument('--ensemblgtf',
                   help='The ensembl gene list for that species')

parser.add_argument('--exppairs',
                   help='List of sample pairs to compare')
parser.add_argument('--bedfile',  nargs='*',
                   help='Additional bed files for feature list')
parser.add_argument('--lazyremap', action="store_true", default=False,
                   help='Skip mapping reads if bam files exit')
parser.add_argument('--nofrag', action="store_true", default=False,
                   help='Omit fragment determination (Used for TGIRT mapping)')
parser.add_argument('--olddeseq', action="store_true", default=False,
                   help='Use old DESeq1 for analysis')
parser.add_argument('--nosizefactors', action="store_true", default=False,
                   help='Don\'t use Deseq size factors in plotting')
'''
parser.add_argument('--maxmismatch',
                   help='Maximum allowed mismatches')
'''
parser.add_argument('--minnontrnasize',type=int,default=20,
                   help='Minimum read length for non-tRNAs')
parser.add_argument('--paironly', action="store_true", default=False,
                   help='Generate only pair files (for adding a pair file after initial processing)')
parser.add_argument('--makehub', action="store_true", default=False,
                   help='make a track hub')
parser.add_argument('--hubonly', action="store_true", default=False,
                   help='Only make the track hub')
parser.add_argument('--maponly', action="store_true", default=False,
                   help='Only do the mapping step')
#parser.add_argument('--maketdr', action="store_true", default=False,
#                   help='create tdrs')
#parser.add_argument('--makeall', action="store_true", default=False,
#                   help='make both track hub and tdrs')
#parser.add_argument('--splittypecounts', action="store_true", default=False,
#                   help='Split type counts into tRNA types')
parser.add_argument('--dumpother', action="store_true", default=False,
                   help='Dump "other" features when counting gene types')
parser.add_argument('--cores',
                   help='number of cores to use')

parser.add_argument('--bamdir',
                   help='directory for placing bam files (default current working directory)')


rlogname = "Rlog.txt"
rlogfile = open(rlogname, "w")

def runrscript(*script):
    print >>sys.stderr, "Rscript "+" ".join(script)
    print >>rlogfile, "*******************************************************" 
    print >>rlogfile, "Rscript "+" ".join(script)
    rlogfile.flush()
    retcode = subprocess.call("Rscript "+" ".join(script), shell=True, stdout = rlogfile, stderr = subprocess.STDOUT)

    if retcode > 0:
        print >>rlogfile, script[0]+" failed"
        print >>sys.stderr, "R script "+script[0]+" failed"
        print >>sys.stderr, "Check "+rlogname+" for details"
        
        #sys.exit()
    return retcode
    

class trnadatabase:
    def __init__(self, dbname):
        self.dbname = dbname
        self.trnatable = dbname+"-trnatable.txt"
        self.bowtiedb = dbname+"-tRNAgenome"
        self.locifile = dbname+"-trnaloci.bed"
        self.maturetrnas=dbname+"-maturetRNAs.bed"
        self.trnaalign = dbname+"-trnaalign.stk"
        self.locialign = dbname+"-trnaloci.stk"
        self.trnafasta = dbname+"-maturetRNAs.fa"
        self.modomics = dbname+"-modomics.txt"
        self.otherseqs = dbname+"-otherseqs.txt"
        self.dbinfo = dbname+"-dbinfo.txt"
    def test(self):
        bowtie2job = subprocess.Popen(["bowtie2","-x",self.bowtiedb, "-U", scriptdir+"test.fq"],stdout = subprocess.PIPE,stderr = subprocess.STDOUT )
        rstatsresults = rstatsjob.communicate()[0]
        if bowtie2job.returncode  != 0:
                print >>sys.stderr, "bowtie2 failed to run"
    def getorgtype(self):
        orgtype = "euk"
        for currline in open(self.dbinfo):
            fields = currline.split()
            if fields[0] == "orgmode":
                orgtype = fields[1]
        return orgtype
            
class expdatabase:
    def __init__(self, expname):
        self.expname = expname
        self.mapinfo = expname+"/"+expname+"-mapinfo.txt"
        self.mapplot = expname+"/"+expname+"-mapinfo.pdf"

        self.trnamapfile = expname+"/"+expname+"-trnamapinfo.txt"
        self.trnamapplot = expname+"/"+expname+"-trnamapinfo.pdf"
        
        
        self.maplog = expname+"/"+expname+"-mapstats.txt"
        self.genetypes = expname+"/"+expname+"-genetypes.txt"
        self.genecounts = expname+"/"+expname+"-readcounts.txt"
        self.trnacounts = expname+"/"+expname+"-trnacounts.txt"
        
        self.normalizedcounts = expname+"/"+expname+"-normalizedreadcounts.txt"
        self.sizefactors = expname+"/"+expname+"-SizeFactors.txt"

        self.genetypecounts=expname+"/"+expname+"-typecounts.txt"
        self.genetypeplot=expname+"/"+expname+"-typecounts.pdf"

        self.genetyperealcounts=expname+"/"+expname+"-typerealcounts.txt"
        self.genetyperealplot=expname+"/"+expname+"-typerealcounts.pdf"
        
        self.trnaaminofile=expname+"/"+expname+"-aminocounts.txt"
        self.trnaaminoplot=expname+"/"+expname+"-aminocounts.pdf"
        self.trnaaminorealplot=expname+"/"+expname+"-aminorealcounts.pdf"
        
        
        self.trnaanticodonfile=expname+"/"+expname+"-anticodoncounts.txt"
        
        self.trnalengthfile=expname+"/"+expname+"-readlengths.txt"
        self.trnalengthplot=expname+"/"+expname+"-readlengths.pdf"
        
        self.trnacoveragefile=expname+"/"+expname+"-coverage.txt"
        self.trnacoverageplot=expname+"/"+expname+"-coverage.pdf"
        self.trnacombinecoverageplot=expname+"/"+expname+"-combinecoverage.pdf"

        self.locicoveragefile=expname+"/pretRNAs/"+expname+"-pretRNAcoverage.txt"
        self.locicoverageplot=expname+"/pretRNAs/"+expname+"-pretRNAcoverage.pdf"
        self.locicombinecoverageplot=expname+"/pretRNAs/"+expname+"-pretRNAcombinecoverage.pdf"
        
        self.trnamismatchfile = expname+"/mismatch/"+expname+"-mismatchcoverage.txt"
        self.trnamismatchplot = expname+"/mismatch/"+expname+"-mismatchcoverage.pdf"
        
        self.trnadeletefile = expname+"/mismatch/"+expname+"-deletecoverage.txt"
        self.trnadeleteplot = expname+"/mismatch/"+expname+"-deletecoverage.pdf"
        
        self.trnamismatchreport = expname+"/mismatch/"+expname+"-mismatchreport.txt"
        self.trnauniquefile=expname+"/unique/"+expname+"-trnauniquecounts.txt"
        self.trnaendfile=expname+"/"+expname+"-trnaendcounts.txt"
        
        
        self.pcaplot = expname+"/"+expname+"-pca.pdf"
        self.pcatrnaplot = expname+"/"+expname+"-pcatrna.pdf"
        
        self.qaoutputname = expname+"/"+expname+"-qa.html"
        
        
        

def mapsamples(samplefile, trnainfo,expinfo, lazyremap, bamdir = "./",  cores = 8, minnontrnasize = 20):
    mapreads.testmain(samplefile=samplefile, trnafile=trnainfo.trnatable,bowtiedb=trnainfo.bowtiedb, bamdir = bamdir, otherseqs = trnainfo.otherseqs,logfile=expinfo.maplog,mapfile=expinfo.mapinfo,trnamapfile = expinfo.trnamapfile,lazy=lazyremap, cores = cores,minnontrnasize = minnontrnasize)
def countfeatures(samplefile, trnainfo,expinfo, ensgtf, bedfiles,  bamdir = "./", cores = 8):
    countreads.testmain(samplefile=samplefile,ensemblgtf=ensgtf,maturetrnas=[trnainfo.maturetrnas], bamdir = bamdir, otherseqs = trnainfo.otherseqs,trnaloci=[trnainfo.locifile],removepseudo=True,genetypefile=expinfo.genetypes,trnatable=trnainfo.trnatable,countfile=expinfo.genecounts,bedfile=bedfiles, trnacounts = expinfo.trnacounts,trnaends = expinfo.trnaendfile,trnauniquecounts = expinfo.trnauniquefile,nofrag=nofrag, cores = cores)
    runrscript(scriptdir+"/pcareadcounts.R",expinfo.genecounts,samplefile,expinfo.pcaplot)
    runrscript(scriptdir+"/pcareadcounts.R",expinfo.trnacounts,samplefile,expinfo.pcatrnaplot)

def counttypes(samplefile, trnainfo,expinfo, ensgtf, bedfiles, bamdir = "./",  ignoresizefactors = False, countfrags = False, bamnofeature = False, cores = 8):
    if not ignoresizefactors:
         
        countreadtypes.testmain(sizefactors=expinfo.sizefactors,combinereps= True , bamdir = bamdir, otherseqs = trnainfo.otherseqs, samplefile=samplefile,maturetrnas=[trnainfo.maturetrnas],trnatable=trnainfo.trnatable,trnaaminofile=expinfo.trnaaminofile,trnaanticodonfile = expinfo.trnaanticodonfile,ensemblgtf=ensgtf,trnaloci=[trnainfo.locifile],countfile=expinfo.genetypecounts,realcountfile=expinfo.genetyperealcounts,bedfile= bedfiles,readlengthfile =  expinfo.trnalengthfile ,countfrags=countfrags, bamnofeature = bamnofeature, cores = cores)
        #Plot reads by gene type and tRNAs by amino acid
        runrscript(scriptdir+"/featuretypes.R",expinfo.genetypecounts,expinfo.genetypeplot)
        runrscript(scriptdir+"/featuretypes.R",expinfo.trnaaminofile,expinfo.trnaaminoplot)
        runrscript(scriptdir+"/featuretypesreal.R",expinfo.trnaaminofile,expinfo.trnaaminorealplot)
        
        runrscript(scriptdir+"/featuretypesreal.R",expinfo.genetyperealcounts,expinfo.genetyperealplot)

        runrscript(scriptdir+"/readlengthhistogram.R",expinfo.trnalengthfile,expinfo.trnalengthplot)
    else:
        countreadtypes.testmain(combinereps= True ,samplefile=samplefile,maturetrnas=[trnainfo.maturetrnas],otherseqs = expinfo.otherseqs, bamdir = bamdir, trnatable=trnainfo.trnatable,trnaaminofile=expinfo.trnaaminofile,trnaanticodonfile = expinfo.trnaanticodonfile,ensemblgtf=ensgtf,trnaloci=[trnainfo.locifile],countfile=expinfo.genetypecounts,realcountfile=expinfo.genetyperealcounts,bedfile= bedfiles,readlengthfile =  expinfo.trnalengthfile,countfrags=countfrags, cores = cores)
        #Plot reads by gene type and tRNAs by amino acid
        runrscript(scriptdir+"/featuretypes.R",expinfo.genetypecounts,expinfo.genetypeplot)
        runrscript(scriptdir+"/featuretypes.R",expinfo.trnaaminofile,expinfo.trnaaminoplot)
        runrscript(scriptdir+"/featuretypesreal.R",expinfo.trnaaminofile,expinfo.trnaaminorealplot)
        runrscript(scriptdir+"/featuretypesreal.R",expinfo.genetyperealcounts,expinfo.genetyperealplot)

        runrscript(scriptdir+"/readlengthhistogram.R",expinfo.trnalengthfile,expinfo.trnalengthplot)
        
def gettrnacoverage(samplefile, trnainfo,expinfo, bamdir = "./",  orgtype = "euk",ignoresizefactors = False, cores = 8):
    #print >>sys.stderr, orgtype
    if not ignoresizefactors:
        getcoverage.testmain(samplefile=samplefile,bedfile=[trnainfo.maturetrnas],locibed=[trnainfo.locifile],locistk=trnainfo.locialign, bamdir = bamdir, lociedgemargin=30,sizefactors=expinfo.sizefactors,orgtype = orgtype,locicoverage=expinfo.locicoveragefile,stkfile=trnainfo.trnaalign, allcoverage=expinfo.trnacoveragefile,trnafasta = trnainfo.trnafasta, cores = cores)
        runrscript(scriptdir+"/newcoverageplots.R","--cov="+expinfo.trnacoveragefile,"--locicov="+expinfo.locicoveragefile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--allcov="+expinfo.trnacoverageplot,"--runname="+expname,"--modomics="+trnainfo.modomics,"--combinecov="+expinfo.trnacombinecoverageplot,"--directory="+expname)
        runrscript(scriptdir+"/boxplotmismatches.R","--runname="+expinfo.expname,"--mismatch="+expinfo.trnacoveragefile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--directory="+expname+"/mismatch/")
    else:
        getcoverage.testmain(samplefile=samplefile,bedfile=[trnainfo.maturetrnas],stkfile=trnainfo.trnaalign,uniquename=expname+"/"+expname,orgtype = orgtype, bamdir = bamdir, allcoverage=expinfo.trnacoveragefile,trnafasta = trnainfo.trnafasta, cores = cores)
        runrscript(scriptdir+"/newcoverageplots.R","--cov="+expinfo.trnacoveragefile,"--locicov="+expinfo.locicoveragefile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--allcov="+expinfo.trnacoverageplot,"--runname="+expname,"--modomics="+trnainfo.modomics,"--combinecov="+expinfo.trnacombinecoverageplot,"--directory="+expname)
        
        runrscript(scriptdir+"/boxplotmismatches.R","--runname="+expinfo.expname,"--mismatch="+expinfo.trnacoveragefile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--directory="+expname+"/mismatch/")
'''
def getendscoverage(samplefile, trnainfo,expinfo, ignoresizefactors = False):
    if not ignoresizefactors:
        getends.main(samplefile=samplefile,bedfile=[trnainfo.maturetrnas],sizefactors=expinfo.sizefactors,stkfile=trnainfo.trnaalign,uniquename=expname+"/mismatch/"+expname, allmismatch=expinfo.trnamismatchfile,trnafasta = trnainfo.trnafasta,mismatchfile=expinfo.trnamismatchfile,mismatchreport=expinfo.trnamismatchreport, indelfile=expinfo.trnadeletefile)
        runrscript(scriptdir+"/endplots.R","--cov="+expinfo.trnamismatchfile,"--mismatchcov="+expinfo.trnamismatchfile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--allcov="+expinfo.trnamismatchplot,"--uniquename="+expname+"/mismatch/"+expname,"--modomics="+trnainfo.modomics,"--directory="+expname+"/mismatch/")
        runrscript(scriptdir+"/boxplotmismatches.R","--mismatch="+expinfo.trnamismatchreport,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--directory="+expname+"/mismatch/")
    else:
        getends.main(samplefile=samplefile,bedfile=[trnainfo.maturetrnas],stkfile=trnainfo.trnaalign,uniquename=expname+"/mismatch/"+expname, allmismatch=expinfo.trnamismatchfile,trnafasta = trnainfo.trnafasta,mismatchfile=expinfo.trnamismatchfile,mismatchreport=expinfo.trnamismatchreport )
        runrscript(scriptdir+"/endplots.R","--cov="+expinfo.trnamismatchfile,"--mismatchcov="+expinfo.trnamismatchfile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--allcov="+expinfo.trnamismatchplot,"--uniquename="+expname+"/mismatch/mismatch/"+expname,"--modomics="+trnainfo.modomics,"--directory="+expname+"/mismatch/")
        runrscript(scriptdir+"/boxplotmismatches.R","--mismatch="+expinfo.trnamismatchreport,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--directory="+expname+"/mismatch/")
def getlocuscoverage(samplefile, trnainfo,expinfo, ignoresizefactors = False):
    if not ignoresizefactors:
        getcoverage.main(samplefile=samplefile ,bedfile=[trnainfo.locifile],sizefactors=expinfo.sizefactors,stkfile=trnainfo.locialign,edgemargin=30, uniquegenome=expname+"/"+expname+"loci",allcoverage=expinfo.locicoveragefile) #removed minextend = 5
        runrscript(scriptdir+"/locuscoverage.R", "--cov="+expinfo.locicoveragefile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--allcov="+expinfo.locicoverageplot,"--combinecov="+expinfo.locicombinecoverageplot,"--directory="+expname)
    else:
        getcoverage.main(samplefile=samplefile ,bedfile=[trnainfo.locifile],stkfile=trnainfo.locialign,edgemargin=30, uniquegenome=expname+"/"+expname+"loci",allcoverage=expinfo.locicoveragefile)
        runrscript(scriptdir+"/locuscoverage.R", "--cov="+expinfo.locicoveragefile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--allcov="+expinfo.locicoverageplot,"--combinecov="+expinfo.locicombinecoverageplot,"--directory="+expname)
'''
def gettdrinfo(samplefile, dbname,expname):
    
    tdrcommand = " ".join(["bash",scriptdir+"/"+"tdrtrax.bash", samplefile, dbname,expname+"/"+expname+"tdrs"])
    print >>sys.stderr, tdrcommand
    tdrjob = subprocess.Popen(tdrcommand,stdout = subprocess.PIPE,stderr = subprocess.STDOUT, shell=True )
    print >>sys.stderr, tdrjob.communicate()[0]
    
def createtrackhub(samplefile, trnainfo,expinfo):
    maketrackhub.main(genomedatabase=trnainfo, samplefile=samplefile,expname=expinfo.expname)
def gettraxqc(samplefile, trnainfo,expinfo, tgirtmode = False):
    traxqc.main(samplefile=samplefile,databasename=trnainfo.dbname,experimentname=expinfo.expname,tgirt = tgirtmode, output=expinfo.qaoutputname)


        

args = parser.parse_args()
dbname = args.databasename
expname = args.experimentname
pairfile =  args.exppairs
ensgtf = args.ensemblgtf
samplefilename = args.samplefile
lazyremap = args.lazyremap
bedfiles= args.bedfile
nofrag= args.nofrag
nosizefactors = args.nosizefactors
olddeseq = args.olddeseq
bamdir = args.bamdir
maponly = args.maponly

mismatch = False
paironly= args.paironly
splittypecounts = False
bamnofeature = args.dumpother

minnontrnasize = args.minnontrnasize

if args.cores is None:
    cores = min(8,cpu_count())
else:
    cores = int(args.cores)

hubonly = args.hubonly

makehubs = args.makehub 
maketdrs=  False # args.maketdr

'''
if args.makeall:
    makehubs = True
    maketdrs = True
'''

scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"



def testsamtools(): #Version: 1.6
    samversionre = re.compile(r"Version\:\s*([\.\d]+)")
    samtoolsloc = get_location("samtools")
    if samtoolsloc is None:
            print >>sys.stderr, "Cannot find samtools in path"
            print >>sys.stderr, "Make sure samtools is installed"
    samtoolsjob = subprocess.Popen([samtoolsloc,"--help"],stdout = subprocess.PIPE,stderr = subprocess.STDOUT )
    samtoolsresults = samtoolsjob.communicate()[0]
    if samtoolsjob.returncode  != 0:
            print >>sys.stderr, "Samtools failed to run"
            print >>sys.stderr, "Make sure samtools is functioning" 
    samtoolsres = samversionre.search(samtoolsresults)
    if samtoolsres:
        if LooseVersion(samtoolsres.group(1)) < LooseVersion("1.0.0"):
            print >>sys.stderr, "Old samtools version "+samtoolsres.group(1)+" found"
            print >>sys.stderr, "Upgrade to latest version"
            sys.exit(1)
    else:
        print >>sys.stderr, "Could not find samtools version number"
        
def testrstats():
    rstatsversionre = re.compile(r"R\s+version\s+((\d+)\.(\d+)\.(\d+))")
    rstatsloc = get_location("R")
    if rstatsloc is None:
            print >>sys.stderr, "Cannot find R in path"
            print >>sys.stderr, "Make sure R is installed"
            sys.exit(1)
    rstatsjob = subprocess.Popen([rstatsloc, "--version"],stdout = subprocess.PIPE,stderr = subprocess.STDOUT )
    rstatsresults = rstatsjob.communicate()[0]
    if rstatsjob.returncode  != 0:
            print >>sys.stderr, "R failed to run"
            print >>sys.stderr, "Make sure R is functioning" 
    rstatsres = rstatsversionre.search(rstatsresults)
    if rstatsres:
        if LooseVersion(rstatsres.group(1)) < LooseVersion("3.1.2"):
            print >>sys.stderr, "Old R version "+rstatsres.group(1)+" found"
            print >>sys.stderr, "Upgrade to latest version"
            sys.exit(1)
    else:
        print >>sys.stderr, "Could not find R version number"


        
testrstats()
get_location("Rscript")

testsamtools()
get_location("bowtie2")

gitversion, gitversionhash = getgithash(scriptdir)

#trnainfo.test(trnainfo)


sampledata = samplefile(samplefilename)
samples = sampledata.getsamples()
if len(samples) < len(set(samples)):
    print >>sys.stderr, "duplicate sample names in sample file"
    sys.exit(1)
for currsample in samples:
    if '-' in currsample:
        print >>sys.stderr, "Sample names containing '-' character are not allowed"
        sys.exit(1)
    if currsample[0].isdigit():
        print >>sys.stderr, "Sample names starting with digits are not allowed"
        sys.exit(1)
replicates = sampledata.allreplicates()
for currsample in replicates:
    
    if '-' in currsample:
        print >>sys.stderr, "Sample names containing '-' character are not allowed"
        sys.exit(1)
    if currsample[0].isdigit():
        print >>sys.stderr, "Sample names starting with digits are not allowed"
        sys.exit(1)
        
replicates = set(replicates)        
       

if pairfile is not None:
    
    missingnames = set()
    for fir, sec in getpairfile(pairfile):
        #print >>sys.stderr, "**"
        if fir not in replicates:
            missingnames.add(fir)
        if sec not in replicates:
            missingnames.add(sec)
    if len(missingnames) > 0:
        print >>sys.stderr, "Pair names "+",".join(missingnames)+" are not present in sample file"
        sys.exit(1)
        
#sys.exit(0)
deseqversion = "DESeq2"
if olddeseq:
    deseqversion = "DESeq"
if runrscript(scriptdir+"checkRmodules.R",deseqversion) > 0:
    print >>sys.stderr, "Not all R modules needed are installed"
    print >>sys.stderr, "check README for needed R modules"
    sys.exit(1)
    





#mkdir -p expname
if not os.path.exists(expname):
    os.makedirs(expname)
'''
if not os.path.exists(expname+"/indiv"):
    os.makedirs(expname+"/indiv")
'''
if not os.path.exists(expname+"/mismatch"):
    os.makedirs(expname+"/mismatch")
if not os.path.exists(expname+"/pretRNAs"):
    os.makedirs(expname+"/pretRNAs")
if not os.path.exists(expname+"/unique"):
    os.makedirs(expname+"/unique")

    
    
if bedfiles is None:   
    bedfiles = list()
dbname = os.path.expanduser(dbname)
if ensgtf is not None:
    ensgtf = os.path.expanduser(ensgtf)
if bedfiles is not None:    
    bedfiles = list(os.path.expanduser(curr) for curr in bedfiles)


trnainfo = trnadatabase(dbname)
orgtype = trnainfo.getorgtype()
#print >>sys.stderr, orgtype
expinfo = expdatabase(expname)
getsamples = samplefile(samplefilename)
if len(getsamples.getsamples()) == 1:
    nosizefactors = True
    

#if only pairfile
if pairfile and paironly:
    if olddeseq:
        deseqret = runrscript(scriptdir+"/deseq1.R",expname,expinfo.genecounts,samplefilename)
        if deseqret == 2:
            print >>sys.stderr, "Deseq analysis failed, cannot continue"
            sys.exit(1)    
    else:
        print >>sys.stderr, scriptdir+"/analyzecounts.R",expname,expinfo.genecounts,samplefilename

        deseqret = runrscript(scriptdir+"/analyzecounts.R",expname,expinfo.genecounts,samplefilename, pairfile)

        if deseqret == 2:
            print >>sys.stderr, "Deseq analysis failed, cannot continue"
            sys.exit(1)
    
    runrscript(scriptdir+"/makescatter.R",expname,expinfo.normalizedcounts,trnainfo.trnatable,expinfo.genetypes,samplefilename,pairfile)
    sys.exit(0)
elif paironly:
    print >>sys.stderr, "pair only mode used but no --pairfile used"
    sys.exit(1)

if hubonly:
    print >>sys.stderr, "Creating trackhub"      

    createtrackhub(samplefilename, dbname,expinfo)
    sys.exit(0)
#getendscoverage(samplefilename, trnainfo,expinfo, nosizefactors)
#
#
#gettdrinfo(samplefilename, dbname,expname)
        #tdrtrax.bash samplefile.txt traxdb outputname
        
#coverage plot of tRNAs
      

#Map the reads
runtime = time.time()
loctime = time.localtime(runtime)
print >>sys.stderr, "Mapping Reads"
#need to check here for names with dashes
mapsamples(samplefilename, trnainfo,expinfo, lazyremap, bamdir = bamdir, cores = cores, minnontrnasize = minnontrnasize)

runinfoname = expname+"/"+expname+"-runinfo.txt"
dbinfo = None
if not lazyremap:
    dbinfo = open(runinfoname,"w")
    print >>dbinfo, "Starting"

else:
    dbinfo = open(runinfoname,"a")
    print >>dbinfo, "---------------------------------------------------------"
    print >>dbinfo, "redoing"
    
print >>dbinfo, "expname\t"+expname
print >>dbinfo, "time\t"+str(runtime)+" ("+str(loctime[1])+"/"+str(loctime[2])+"/"+str(loctime[0])+")"
print >>dbinfo, "samplefile\t"+os.path.realpath(samplefilename)
print >>dbinfo, "dbname\t"+os.path.realpath(dbname)
print >>dbinfo, "git version\t"+gitversion

print >>dbinfo, "git version hash\t"+gitversionhash

print >>dbinfo, "command\t"+" ".join(sys.argv)
dbinfo.close()

runrscript(scriptdir+"/featuretypes.R",expinfo.mapinfo,expinfo.mapplot)
runrscript(scriptdir+"/featuretypes.R",expinfo.trnamapfile,expinfo.trnamapplot)

#print >>sys.stderr, "Counting Read Types"
#counttypes(samplefilename, trnainfo,expinfo, ensgtf, bedfiles, ignoresizefactors = nosizefactors)

if maponly:
    sys.exit()
#Count the reads for DEseq2 and scatter plots
print >>sys.stderr, "Counting Reads"
countfeatures(samplefilename, trnainfo,expinfo, ensgtf, bedfiles, bamdir = bamdir,  cores = cores)
#Create a plot of mapped reads                                
print >>sys.stderr, "Analyzing counts"



#Analyze counts and create scatter plots if pair file is provided
if pairfile:
    if olddeseq:
        deseqret = runrscript(scriptdir+"/deseq1.R",expname,expinfo.genecounts,samplefilename)
        if deseqret >  0:
            print >>sys.stderr, "Deseq analysis failed, cannot continue"
            sys.exit(1)    
    else:
        deseqret = runrscript(scriptdir+"/analyzecounts.R",expname,expinfo.genecounts,samplefilename, pairfile)
        print >>sys.stderr, scriptdir+"/analyzecounts.R",expname,expinfo.genecounts,samplefilename, pairfile

        if deseqret > 0:
            print >>sys.stderr, "Deseq analysis failed, cannot continue"
            print >>sys.stderr, "Likely that a sample did not contain enough reads"
            sys.exit(1)
    
    runrscript(scriptdir+"/makescatter.R",expname,expinfo.normalizedcounts,trnainfo.trnatable,expinfo.genetypes,samplefilename,pairfile)
    runrscript(scriptdir+"/endplot.R","--ends="+expinfo.trnaendfile,"--trna="+trnainfo.trnatable, "--samples="+samplefilename,"--directory="+expname+"/mismatch/")
    

elif not nosizefactors:
    if olddeseq:
        deseqret = runrscript(scriptdir+"/deseq1.R",expname,expinfo.genecounts,samplefilename)
        if deseqret == 2:
            print >>sys.stderr, "Deseq analysis failed, cannot continue"
            sys.exit(1)    
    else:
        deseqret = runrscript(scriptdir+"/analyzecounts.R",expname,expinfo.genecounts,samplefilename) 
        print >>sys.stderr, scriptdir+"/analyzecounts.R",expname,expinfo.genecounts,samplefilename
        if deseqret == 2:
            print >>sys.stderr, "Deseq analysis failed, cannot continue"
            sys.exit(1)
#Count the reads by gene type
print >>sys.stderr, "Counting Read Types"
counttypes(samplefilename, trnainfo,expinfo, ensgtf, bedfiles,  bamdir = bamdir, ignoresizefactors = nosizefactors,countfrags =  splittypecounts, bamnofeature = bamnofeature, cores = cores)


#coverage plot of tRNAs
print >>sys.stderr, "Generating Read Coverage plots"      
#print >>sys.stderr, orgtype
gettrnacoverage(samplefilename, trnainfo,expinfo, bamdir = bamdir, orgtype = orgtype, ignoresizefactors = nosizefactors, cores = cores)

#coverage plot of pre-tRNAs
#getlocuscoverage(samplefilename, trnainfo,expinfo, nosizefactors)

#coverage plot of mismatches
#getmismatchcoverage(samplefilename, trnainfo,expinfo, nosizefactors)
#print >>sys.stderr, "Counting mismatches"      

#getendscoverage(samplefilename, trnainfo,expinfo, nosizefactors)

gettraxqc(samplefilename, trnainfo, expinfo, tgirtmode = nofrag)

if makehubs:
    print >>sys.stderr, "Creating trackhub"      

    createtrackhub(samplefilename, dbname,expinfo)


if (os.path.isfile(scriptdir+"/"+"tdrtrax.bash") and maketdrs):
    print >>sys.stderr, "Creating tdrs"      

    gettdrinfo(samplefilename, dbname,expname)
        #tdrtrax.bash samplefile.txt traxdb outputname


