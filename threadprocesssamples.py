#!/usr/bin/env python3

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
parser.add_argument('--bedfile',  nargs='+', default=list(),
                   help='Additional bed files for feature list')
parser.add_argument('--lazyremap', action="store_true", default=False,
                   help='Skip mapping reads if bam files exit')
parser.add_argument('--nofrag', action="store_true", default=False,
                   help='Omit fragment determination (Used for TGIRT mapping)')
parser.add_argument('--olddeseq', action="store_true", default=False,
                   help='Use old DESeq1 for analysis')
parser.add_argument('--nosizefactors', action="store_true", default=False,
                   help='Don\'t use Deseq size factors in plotting')
parser.add_argument('--maxmismatch',
                   help='Maximum allowed mismatches')
parser.add_argument('--mismatch',
                   help='Generate (experimental) mismatch coverage charts')

parser.add_argument('--paironly', action="store_true", default=False,
                   help='Generate only pair files (for adding a pair file after initial processing)')
parser.add_argument('--makehub', action="store_true", default=False,
                   help='make a track hub')
parser.add_argument('--hubonly', action="store_true", default=False,
                   help='make only the track hub')
parser.add_argument('--maketdr', action="store_true", default=False,
                   help='create tdrs')
parser.add_argument('--makeall', action="store_true", default=False,
                   help='make both track hub and tdrs')
parser.add_argument('--splittypecounts', action="store_true", default=False,
                   help='Split type counts into tRNA types')
parser.add_argument('--cores',
                   help='number of cores to use')


rlogname = "Rlog.txt"
rlogfile = open(rlogname, "w")

def runrscript(*script):
    print("Rscript "+" ".join(script), file=sys.stderr)
    print("*******************************************************", file=rlogfile) 
    print("Rscript "+" ".join(script), file=rlogfile)
    rlogfile.flush()
    retcode = subprocess.call("Rscript "+" ".join(script), shell=True, stdout = rlogfile, stderr = subprocess.STDOUT)

    if retcode > 0:
        print(script[0]+" failed", file=rlogfile)
        print("R script "+script[0]+" failed", file=sys.stderr)
        print("Check "+rlogname+" for details", file=sys.stderr)
        
        #sys.exit()
    return retcode
    

class trnadatabase:
    def __init__(self, dbname):
        self.trnatable = dbname+"-trnatable.txt"
        self.bowtiedb = dbname+"-tRNAgenome"
        self.locifile = dbname+"-trnaloci.bed"
        self.maturetrnas=dbname+"-maturetRNAs.bed"
        self.trnaalign = dbname+"-trnaalign.stk"
        self.locialign = dbname+"-trnaloci.stk"
        self.trnafasta = dbname+"-maturetRNAs.fa"
        self.modomics = dbname+"-modomics.txt"
        self.otherseqs = dbname+"-otherseqs.txt"
    def test(self):
        bowtie2job = subprocess.Popen(["bowtie2","-x",self.bowtiedb, "-U", scriptdir+"test.fq"],stdout = subprocess.PIPE,stderr = subprocess.STDOUT )
        rstatsresults = rstatsjob.communicate()[0]
        if bowtie2job.returncode  != 0:
                print("bowtie2 failed to run", file=sys.stderr)

        
class expdatabase:
    def __init__(self, expname):
        self.mapinfo = expname+"/"+expname+"-mapinfo.txt"
        self.mapplot = expname+"/"+expname+"-mapinfo.pdf"
        
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
        
        self.trnalengthfile=expname+"/"+expname+"-readlengths.txt"
        self.trnalengthplot=expname+"/"+expname+"-readlengths.pdf"
        
        self.trnacoveragefile=expname+"/"+expname+"-coverage.txt"
        self.trnacoverageplot=expname+"/"+expname+"-coverage.pdf"
        self.trnacombinecoverageplot=expname+"/"+expname+"-combinecoverage.pdf"

        self.locicoveragefile=expname+"/"+expname+"-pretRNAcoverage.txt"
        self.locicoverageplot=expname+"/pretRNAs/"+expname+"-pretRNAcoverage.pdf"
        self.locicombinecoverageplot=expname+"/pretRNAs/"+expname+"-pretRNAcombinecoverage.pdf"
        
        self.trnamismatchfile = expname+"/mismatch/"+expname+"-mismatchcoverage.txt"
        self.trnamismatchplot = expname+"/mismatch/"+expname+"-mismatchcoverage.pdf"
        
        self.trnadeletefile = expname+"/mismatch/"+expname+"-deletecoverage.txt"
        self.trnadeleteplot = expname+"/mismatch/"+expname+"-deletecoverage.pdf"
        
        self.trnamismatchreport = expname+"/mismatch/"+expname+"-mismatchreport.txt"
        self.trnauniquefile=expname+"/"+expname+"-trnauniquecounts.txt"
        
        self.pcaplot = expname+"/"+expname+"-pca.pdf"
        self.pcatrnaplot = expname+"/"+expname+"-pcatrna.pdf"
        
        
        

def mapsamples(samplefile, trnainfo,expinfo, lazyremap, cores = 8):
    mapreads.testmain(samplefile=samplefile, trnafile=trnainfo.trnatable,bowtiedb=trnainfo.bowtiedb,otherseqs = trnainfo.otherseqs,logfile=expinfo.maplog,mapfile=expinfo.mapinfo, lazy=lazyremap, cores = cores)
def countfeatures(samplefile, trnainfo,expinfo, ensgtf, bedfiles, cores = 8):
    countreads.testmain(samplefile=samplefile,ensemblgtf=ensgtf,maturetrnas=[trnainfo.maturetrnas],otherseqs = trnainfo.otherseqs,trnaloci=[trnainfo.locifile],removepseudo=True,genetypefile=expinfo.genetypes,trnatable=trnainfo.trnatable,countfile=expinfo.genecounts,bedfile=bedfiles, trnacounts = expinfo.trnacounts,trnauniquecounts = expinfo.trnauniquefile,nofrag=nofrag, cores = cores)
    runrscript(scriptdir+"/pcareadcounts.R",expinfo.genecounts,samplefile,expinfo.pcaplot)
    runrscript(scriptdir+"/pcareadcounts.R",expinfo.trnacounts,samplefile,expinfo.pcatrnaplot)

def counttypes(samplefile, trnainfo,expinfo, ensgtf, bedfiles, ignoresizefactors = False, countfrags = False, cores = 8):
    if not ignoresizefactors:
        
        countreadtypes.testmain(sizefactors=expinfo.sizefactors,combinereps= True ,otherseqs = trnainfo.otherseqs, samplefile=samplefile,maturetrnas=[trnainfo.maturetrnas],trnatable=trnainfo.trnatable,trnaaminofile=expinfo.trnaaminofile,ensemblgtf=ensgtf,trnaloci=[trnainfo.locifile],countfile=expinfo.genetypecounts,realcountfile=expinfo.genetyperealcounts,bedfile= bedfiles,readlengthfile =  expinfo.trnalengthfile ,countfrags=countfrags, cores = cores)
        #Plot reads by gene type and tRNAs by amino acid
        runrscript(scriptdir+"/featuretypes.R",expinfo.genetypecounts,expinfo.genetypeplot)
        runrscript(scriptdir+"/featuretypes.R",expinfo.trnaaminofile,expinfo.trnaaminoplot)
        runrscript(scriptdir+"/featuretypesreal.R",expinfo.trnaaminofile,expinfo.trnaaminoplot)

        runrscript(scriptdir+"/readlengthhistogram.R",expinfo.trnalengthfile,expinfo.trnalengthplot)
    else:
        countreadtypes.main(combinereps= True ,samplefile=samplefile,maturetrnas=[trnainfo.maturetrnas],otherseqs = expinfo.otherseqs,trnatable=trnainfo.trnatable,trnaaminofile=expinfo.trnaaminofile,ensemblgtf=ensgtf,trnaloci=[trnainfo.locifile],countfile=expinfo.genetypecounts,realcountfile=expinfo.genetyperealcounts,bedfile= bedfiles,readlengthfile =  expinfo.trnalengthfile,countfrags=countfrags, cores = cores)
        #Plot reads by gene type and tRNAs by amino acid
        runrscript(scriptdir+"/featuretypes.R",expinfo.genetypecounts,expinfo.genetypeplot)
        runrscript(scriptdir+"/featuretypes.R",expinfo.trnaaminofile,expinfo.trnaaminoplot)
        runrscript(scriptdir+"/featuretypesreal.R",expinfo.trnaaminofile,expinfo.trnaaminoplot)

        runrscript(scriptdir+"/readlengthhistogram.R",expinfo.trnalengthfile,expinfo.trnalengthplot)
        
def gettrnacoverage(samplefile, trnainfo,expinfo, ignoresizefactors = False, cores = 8):
    if not ignoresizefactors:
        getcoverage.testmain(samplefile=samplefile,bedfile=[trnainfo.maturetrnas],locibed=[trnainfo.locifile],locistk=trnainfo.locialign,lociedgemargin=30,sizefactors=expinfo.sizefactors,locicoverage=expinfo.locicoveragefile,stkfile=trnainfo.trnaalign, allcoverage=expinfo.trnacoveragefile,trnafasta = trnainfo.trnafasta, cores = cores)
        runrscript(scriptdir+"/newcoverageplots.R","--cov="+expinfo.trnacoveragefile,"--locicov="+expinfo.locicoveragefile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--allcov="+expinfo.trnacoverageplot,"--runname="+expname,"--modomics="+trnainfo.modomics,"--combinecov="+expinfo.trnacombinecoverageplot,"--directory="+expname)
        runrscript(scriptdir+"/boxplotmismatches.R","--mismatch="+expinfo.trnacoveragefile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--directory="+expname+"/mismatch/")
    else:
        getcoverage.testmain(samplefile=samplefile,bedfile=[trnainfo.maturetrnas],stkfile=trnainfo.trnaalign,uniquename=expname+"/"+expname, allcoverage=expinfo.trnacoveragefile,trnafasta = trnainfo.trnafasta, cores = cores)
        runrscript(scriptdir+"/newcoverageplots.R","--cov="+expinfo.trnacoveragefile,"--locicov="+expinfo.locicoveragefile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--allcov="+expinfo.trnacoverageplot,"--runname="+expname,"--modomics="+trnainfo.modomics,"--combinecov="+expinfo.trnacombinecoverageplot,"--directory="+expname)
        
        runrscript(scriptdir+"/boxplotmismatches.R","--mismatch="+expinfo.trnacoveragefile,"--trna="+trnainfo.trnatable,"--samples="+samplefile,"--directory="+expname+"/mismatch/")
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
    print(tdrcommand, file=sys.stderr)
    tdrjob = subprocess.Popen(tdrcommand,stdout = subprocess.PIPE,stderr = subprocess.STDOUT, shell=True )
    print(tdrjob.communicate()[0], file=sys.stderr)
    
def createtrackhub(samplefile, trnainfo,expinfo):
    maketrackhub.main(genomedatabase=trnainfo, samplefile=samplefile,expname=expname)

        

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
mismatch = args.mismatch 
paironly= args.paironly
splittypecounts = args.splittypecounts

cores = args.cores
if cores is None:
    cores = min(8,cpu_count())

hubonly = args.hubonly

makehubs = args.makehub 
maketdrs= args.maketdr

if args.makeall:
    makehubs = True
    maketdrs = True


scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"



def testsamtools(): #Version: 1.6
    samversionre = re.compile(r"Version\:\s*([\.\d]+)")
    samtoolsloc = get_location("samtools")
    if samtoolsloc is None:
            print("Cannot find samtools in path", file=sys.stderr)
            print("Make sure samtools is installed", file=sys.stderr)
    samtoolsjob = subprocess.Popen([samtoolsloc,"--help"],stdout = subprocess.PIPE,stderr = subprocess.STDOUT )
    samtoolsresults = samtoolsjob.communicate()[0]
    if samtoolsjob.returncode  != 0:
            print("Samtools failed to run", file=sys.stderr)
            print("Make sure samtools is functioning", file=sys.stderr) 
    samtoolsres = samversionre.search(samtoolsresults)
    if samtoolsres:
        if LooseVersion(samtoolsres.group(1)) < LooseVersion("1.0.0"):
            print("Old samtools version "+samtoolsres.group(1)+" found", file=sys.stderr)
            print("Upgrade to latest version", file=sys.stderr)
            sys.exit(1)
    else:
        print("Could not find samtools version number", file=sys.stderr)
        
def testrstats():
    rstatsversionre = re.compile(r"R\s+version\s+((\d+)\.(\d+)\.(\d+))")
    rstatsloc = get_location("R")
    if rstatsloc is None:
            print("Cannot find R in path", file=sys.stderr)
            print("Make sure R is installed", file=sys.stderr)
            sys.exit(1)
    rstatsjob = subprocess.Popen([rstatsloc, "--version"],stdout = subprocess.PIPE,stderr = subprocess.STDOUT )
    rstatsresults = rstatsjob.communicate()[0]
    if rstatsjob.returncode  != 0:
            print("R failed to run", file=sys.stderr)
            print("Make sure R is functioning", file=sys.stderr) 
    rstatsres = rstatsversionre.search(rstatsresults)
    if rstatsres:
        if LooseVersion(rstatsres.group(1)) < LooseVersion("3.1.2"):
            print("Old R version "+rstatsres.group(1)+" found", file=sys.stderr)
            print("Upgrade to latest version", file=sys.stderr)
            sys.exit(1)
    else:
        print("Could not find R version number", file=sys.stderr)


        
testrstats()
get_location("Rscript")

testsamtools()
get_location("bowtie2")

gitversion, gitversionhash = getgithash(scriptdir)

#trnainfo.test(trnainfo)


sampledata = samplefile(samplefilename)
samples = sampledata.getsamples()
for currsample in samples:
    if '-' in currsample:
        print("Sample names containing '-' character are not allowed", file=sys.stderr)
        sys.exit(1)
    if currsample[0].isdigit():
        print("Sample names starting with digits are not allowed", file=sys.stderr)
        sys.exit(1)
replicates = sampledata.allreplicates()
for currsample in replicates:
    
    if '-' in currsample:
        print("Sample names containing '-' character are not allowed", file=sys.stderr)
        sys.exit(1)
    if currsample[0].isdigit():
        print("Sample names starting with digits are not allowed", file=sys.stderr)
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
        print("Pair names "+",".join(missingnames)+" are not present in sample file", file=sys.stderr)
        sys.exit(1)
        
#sys.exit(0)
deseqversion = "DESeq2"
if olddeseq:
    deseqversion = "DESeq"
if runrscript(scriptdir+"checkRmodules.R",deseqversion) > 0:
    print("Not all R modules needed are installed", file=sys.stderr)
    print("check README for needed R modules", file=sys.stderr)
    sys.exit(1)
    





#mkdir -p expname
if not os.path.exists(expname):
    os.makedirs(expname)
if not os.path.exists(expname+"/indiv"):
    os.makedirs(expname+"/indiv")
if not os.path.exists(expname+"/mismatch"):
    os.makedirs(expname+"/mismatch")
if not os.path.exists(expname+"/pretRNAs"):
    os.makedirs(expname+"/pretRNAs")

    
    
    
dbname = os.path.expanduser(dbname)
if ensgtf is not None:
    ensgtf = os.path.expanduser(ensgtf)
bedfiles = list(os.path.expanduser(curr) for curr in bedfiles)


trnainfo = trnadatabase(dbname)
expinfo = expdatabase(expname)
getsamples = samplefile(samplefilename)
if len(getsamples.getsamples()) == 1:
    nosizefactors = True
    

#if only pairfile
if pairfile and paironly:
    if olddeseq:
        deseqret = runrscript(scriptdir+"/deseq1.R",expname,expinfo.genecounts,samplefilename)
        if deseqret == 2:
            print("Deseq analysis failed, cannot continue", file=sys.stderr)
            sys.exit(1)    
    else:
        print(scriptdir+"/analyzecounts.R",expname,expinfo.genecounts,samplefilename, file=sys.stderr)

        deseqret = runrscript(scriptdir+"/analyzecounts.R",expname,expinfo.genecounts,samplefilename, pairfile)

        if deseqret == 2:
            print("Deseq analysis failed, cannot continue", file=sys.stderr)
            sys.exit(1)
    
    runrscript(scriptdir+"/makescatter.R",expname,expinfo.normalizedcounts,trnainfo.trnatable,expinfo.genetypes,samplefilename,pairfile)
    sys.exit(0)
elif paironly:
    print("pair only mode used but no --pairfile used", file=sys.stderr)
    sys.exit(1)

if hubonly:
    print("Creating trackhub", file=sys.stderr)      

    createtrackhub(samplefilename, dbname,expname)
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
print("Mapping Reads", file=sys.stderr)
#need to check here for names with dashes
mapsamples(samplefilename, trnainfo,expinfo, lazyremap)

runinfoname = expname+"/"+expname+"-runinfo.txt"
dbinfo = None
if not lazyremap:
    dbinfo = open(runinfoname,"w")
    print("Starting", file=dbinfo)

else:
    dbinfo = open(runinfoname,"a")
    print("---------------------------------------------------------", file=dbinfo)
    print("redoing", file=dbinfo)
    
print("expname\t"+expname, file=dbinfo)
print("time\t"+str(runtime)+" ("+str(loctime[1])+"/"+str(loctime[2])+"/"+str(loctime[0])+")", file=dbinfo)
print("samplefile\t"+os.path.realpath(samplefilename), file=dbinfo)
print("dbname\t"+os.path.realpath(dbname), file=dbinfo)
print("git version\t"+gitversion, file=dbinfo)

print("git version hash\t"+gitversionhash, file=dbinfo)
dbinfo.close()

runrscript(scriptdir+"/featuretypes.R",expinfo.mapinfo,expinfo.mapplot)


#print >>sys.stderr, "Counting Read Types"
#counttypes(samplefilename, trnainfo,expinfo, ensgtf, bedfiles, ignoresizefactors = nosizefactors)


#Count the reads for DEseq2 and scatter plots
print("Counting Reads", file=sys.stderr)
countfeatures(samplefilename, trnainfo,expinfo, ensgtf, bedfiles, cores = cores)
#Create a plot of mapped reads                                
print("Analyzing counts", file=sys.stderr)



#Analyze counts and create scatter plots if pair file is provided
if pairfile:
    if olddeseq:
        deseqret = runrscript(scriptdir+"/deseq1.R",expname,expinfo.genecounts,samplefilename)
        if deseqret == 2:
            print("Deseq analysis failed, cannot continue", file=sys.stderr)
            sys.exit(1)    
    else:
        deseqret = runrscript(scriptdir+"/analyzecounts.R",expname,expinfo.genecounts,samplefilename, pairfile)
        print(scriptdir+"/analyzecounts.R",expname,expinfo.genecounts,samplefilename, pairfile, file=sys.stderr)

        if deseqret == 2:
            print("Deseq analysis failed, cannot continue", file=sys.stderr)
            sys.exit(1)
    
    runrscript(scriptdir+"/makescatter.R",expname,expinfo.normalizedcounts,trnainfo.trnatable,expinfo.genetypes,samplefilename,pairfile)
elif not nosizefactors:
    if olddeseq:
        deseqret = runrscript(scriptdir+"/deseq1.R",expname,expinfo.genecounts,samplefilename)
        if deseqret == 2:
            print("Deseq analysis failed, cannot continue", file=sys.stderr)
            sys.exit(1)    
    else:
        deseqret = runrscript(scriptdir+"/analyzecounts.R",expname,expinfo.genecounts,samplefilename)
        print(scriptdir+"/analyzecounts.R",expname,expinfo.genecounts,samplefilename, file=sys.stderr)
        if deseqret == 2:
            print("Deseq analysis failed, cannot continue", file=sys.stderr)
            sys.exit(1)
#Count the reads by gene type
print("Counting Read Types", file=sys.stderr)
counttypes(samplefilename, trnainfo,expinfo, ensgtf, bedfiles, ignoresizefactors = nosizefactors,countfrags =  splittypecounts, cores = cores)


#coverage plot of tRNAs
print("Generating Read Coverage plots", file=sys.stderr)      
gettrnacoverage(samplefilename, trnainfo,expinfo, ignoresizefactors = nosizefactors, cores = cores)

#coverage plot of pre-tRNAs
#getlocuscoverage(samplefilename, trnainfo,expinfo, nosizefactors)

#coverage plot of mismatches
#getmismatchcoverage(samplefilename, trnainfo,expinfo, nosizefactors)
#print >>sys.stderr, "Counting mismatches"      

#getendscoverage(samplefilename, trnainfo,expinfo, nosizefactors)

if makehubs:
    print("Creating trackhub", file=sys.stderr)      

    createtrackhub(samplefilename, dbname,expname)


if (os.path.isfile(scriptdir+"/"+"tdrtrax.bash") and maketdrs):
    print("Creating tdrs", file=sys.stderr)      

    gettdrinfo(samplefilename, dbname,expname)
        #tdrtrax.bash samplefile.txt traxdb outputname


