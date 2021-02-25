#!/usr/bin/env python

import sys
import argparse
import os.path
from collections import defaultdict
from trnasequtils import *
import itertools
import subprocess
import getmaturetrnas
import aligntrnalocus
from distutils.spawn import find_executable
from distutils.version import LooseVersion, StrictVersion
import time
from multiprocessing import Pool, cpu_count


def get_location(program, allowfail = False):
    progloc = find_executable(program)
    if find_executable(program) is None and not allowfail:
        print >>sys.stderr, "Could not find "+program+" in path"
        print >>sys.stderr, "Aborting"
        sys.exit(1)
    else:
        return progloc
        

eukpositions = list([-1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,'17a',18,19,20,'20a','20b',21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14','e15','e16','e17','e18','e19',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76])
archpositions = list([-1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,'17a',18,19,20,'20a','20b',21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76])
bactpositions = list([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14','e15','e16','e17',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76])

mitopositions = list([-1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,'e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14','e15','e16','e17',46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76])



#print >>sys.stderr, len(positions)
#this gets the tRNA numbers by the sprinzel numbering system
def gettrnanums(trnaalign, margin = 0, orgtype = "euk", mode = "transcript"):
    trnanum = list()
    currcount = 0
    enum = 1
    gapnum = 1
    intronnum = 1
    positions = eukpositions 
    if orgtype == "arch":
    	positions = archpositions
    elif orgtype == "mito":
        positions = mitopositions
    elif orgtype == "bact":
        positions = bactpositions 
    for i in range(margin):
        trnanum.append('head'+str(margin - i))
    for i, struct in enumerate(trnaalign.consensus):
        if currcount >= len(positions):
            trnanum.append('gap'+str(gapnum))
            gapnum += 1
            currcount += 1
        elif struct in  set("+=*"):
            #special case to account for differences between loci/transcripts
            if currcount == 0 and struct == '=' and orgtype != "bact":
                currcount = 1
                gapnum = 1
            if positions[currcount] == 'e':
                trnanum.append('e'+str(enum))
                enum += 1
                currcount += 1
                gapnum = 1
            elif positions[currcount] == '-':
                trnanum.append(str(currcount)+'.gap'+str(gapnum))
                gapnum += 1
                currcount += 1
            else:
                trnanum.append(str(positions[currcount]))
                currcount += 1
                gapnum = 1
        else:
            #if intron

            if positions[currcount] == 38:
                trnanum.append('intron'+str(intronnum))
                intronnum += 1
            else:
                
                trnanum.append(str(currcount)+'.gap'+str(gapnum))
                gapnum += 1
    for i in range(margin):
        trnanum.append('tail'+str(i+1))
    return trnanum

parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
parser.add_argument('--databasename',required=True,
                   help='database name to be used')
parser.add_argument('--genomefile',required=True,
                   help='Fasta file containing genome sequence')
parser.add_argument('--trnascanfile',required=True,
                   help='output from tRNAscan-SE run')
parser.add_argument('--addseqs',
                   help='file with additional sets of sequence transcripts')
parser.add_argument('--gtrnafafile',
                   help='Fasta file of tRNA sequences from gtRNAdb')
parser.add_argument('--namemapfile',
                   help='Name mapping from gtRNAdb')
parser.add_argument('--orgmode',
                   help='organism mode (euk/arch/bact/mito, default euk)')
parser.add_argument('--forcecca', action="store_true", default=False,
                   help='force the addition of a CCA tail')

def shellcall(shellcommand,failquit = False):
    retcode = subprocess.call(shellcommand, shell=True)
    if retcode > 0 and failquit:
        print >>sys.stderr, "Command failed:"
        print >>sys.stderr, shellcall
        print >>sys.stderr, "quiting program..."
        sys.exit(1)
    elif retcode > 0:
        print >>sys.stderr, "Command failed:"
        print >>sys.stderr, shellcall
    return retcode
        
        
args = parser.parse_args()
dbname = args.databasename
scanfile = args.trnascanfile
genomefile = args.genomefile
gtrnafafile = args.gtrnafafile
namemapfile = args.namemapfile
forcecca = args.forcecca
if namemapfile is None:
    print >>sys.stderr, "Name map file currently needed for tRNA database creation"
    sys.exit(1)
orgmode = args.orgmode
addseqs = args.addseqs
dbdirectory = os.path.dirname(dbname) + "/"
if dbdirectory == "/":
    dbdirectory = ""
dbname = os.path.basename(dbname)


runtime = time.time()
loctime = time.localtime(runtime)


#$1 is database name
#trnascanfile is trnascan file
#genomefile is fasta file of genome

#test command line programs


    
get_location("samtools")
get_location("bowtie2-build")

if not os.path.isfile(genomefile+".fai"):
    shellcall("samtools faidx "+genomefile)
    
scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"

''''
cmbuild --enone --hand -F

'''
prokmode = False
if orgmode is None:
    orgmode = "euk"

if orgmode == "euk":
    maturemodel =  scriptdir+'trnamature-euk.cm'
    trnamodel =  scriptdir+'TRNAinf-euk.cm'
elif orgmode == "arch":
    maturemodel =  scriptdir+'trnamature-arch.cm'
    trnamodel =  scriptdir+'TRNAinf-arch.cm'
    prokmode = True
elif orgmode == "mito":
    maturemodel =  scriptdir+'TRNAMatureMitoinf.cm'
    trnamodel =  scriptdir+'TRNAinf.cm'
    prokmode = False
elif orgmode == "bact":
    maturemodel =  scriptdir+'trnamature-bact.cm'
    trnamodel =  scriptdir+'TRNAinf-bact.cm'
    prokmode = True
    
if forcecca:
    prokmode = False

transcriptstk = dbdirectory+dbname+"-trnaalign.stk"
locusstk = dbdirectory+dbname+"-trnaloci.stk"
getmaturetrnas.main(trnascan=[scanfile], genome=genomefile,gtrnafa=gtrnafafile,namemap=namemapfile, bedfile=dbdirectory+dbname+"-maturetRNAs.bed",maturetrnatable=dbdirectory+dbname+"-trnatable.txt",trnaalignment=transcriptstk,locibed=dbdirectory+dbname+"-trnaloci.bed",maturetrnafa=dbdirectory+dbname+"-maturetRNAs.fa",cmmodel = maturemodel, prokmode = prokmode)
aligntrnalocus.main(genomefile=genomefile,stkfile=locusstk,trnaloci=dbdirectory+dbname+"-trnaloci.bed", cmmodel = trnamodel)



transalign = list(readrnastk(open(transcriptstk, "r")))[0]
transnums = gettrnanums(transalign, margin = 0, orgtype = orgmode, mode = "transcript")
alignnumfile = open(dbdirectory+dbname+"-alignnum.txt", "w")
for i, currbase in enumerate(transalign.consensus):
    print >>alignnumfile, "\t".join([str(i), transalign.consensus[i], transalign.currstruct[i], transnums[i]])
alignnumfile.close()
locialign = list(readrnastk(open(locusstk, "r")))[0]
locinums = gettrnanums(locialign, margin = 0, orgtype = orgmode, mode = "locus")
locinumfile = open(dbdirectory+dbname+"-locusnum.txt", "w")
for i, currbase in enumerate(locialign.consensus):
    print >>locinumfile, "\t".join([str(i), locialign.consensus[i], locialign.currstruct[i], locinums[i]])


locinumfile.close()








seqfiles = dict()
newseqs = dict()
seqfastaname = ""
if addseqs:
    seqfastaname = dbdirectory+dbname+"-additionals.fa"
    seqfasta = open(seqfastaname, "w")
    seqcounts = defaultdict(int)
    otherseqs = open(dbdirectory+dbname+"-otherseqs.txt", "w")
    for currline in open(addseqs):
        
        fields = currline.split()
        if len(fields) < 2:
            continue
        seqsname = fields[0]   
        seqsfile = fields[1]
        print >>otherseqs, seqsname + "\t"+dbname+"-"+seqsname+"_seq.fa"+"\t"+dbname+"-"+seqsname+"_seq.bed"
        seqfiles[fields[0]] = fields[1]

        seqbed = open(dbdirectory+dbname+"-"+seqsname+"_seq.bed", "w")
        
        for name, seq in readmultifasta(fields[1]):
            
            if name not in seqcounts:
                
                print >>seqfasta, ">"+name
                seqcounts[name] += 1
            else:
                
                print >>seqfasta, ">"+name +"."+str(seqcounts[name])
                seqcounts[name] += 1
            print >>seqfasta, (20 *"N")+seq.upper()+(20 *"N")
        
            print >>seqbed, "\t".join([name, str(20), str(len(seq) + 20), name,"1000", "+"])
        
        seqbed.close()
    otherseqs.close() 
    #print >>sys.stderr, "cat "+" ".join(newseqs.values())+" >"+dbdirectory+dbname+"-additionals.fa"
    #shellcall("cat "+" ".join(newseqs.values())+" >"+dbdirectory+dbname+"-additionals.fa", failquit = True)
    seqfasta.close()
    
#shellcall("cat "+dbdirectory+dbname+"-maturetRNAs.fa "+genomefile+" "+" ".join(seqfiles.values())+" >"+dbdirectory+dbname+"-tRNAgenome.fa", failquit = True)

shellcall("cat "+dbdirectory+dbname+"-maturetRNAs.fa "+genomefile+" "+seqfastaname+" >"+dbdirectory+dbname+"-tRNAgenome.fa", failquit = True)
#sys.exit(1)    
scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"
gitversion, githash = getgithash(scriptdir)
    
    
dbinfo = open(dbdirectory+dbname+ "-dbinfo.txt","w")
print >>dbinfo, "time\t"+str(runtime)+"("+str(loctime[1])+"/"+str(loctime[2])+"/"+str(loctime[0])+")"
print >>dbinfo, "creation\t"+" ".join(sys.argv)
print >>dbinfo, "genomefile\t"+str(genomefile)
print >>dbinfo, "trnascanfile\t"+str(scanfile)
print >>dbinfo, "orgmode\t"+str(orgmode)
print >>dbinfo, "forcecca\t"+str(forcecca)
print >>dbinfo, "git version\t"+str(gitversion)

if addseqs:
    print >>dbinfo, "additional transcripts\t"+" ".join(name+":"+source for name, source in seqfiles.iteritems())

dbinfo.close()
cores = None
if cores is None:
    cores = min(8,cpu_count())
indexoption = ""
if True:
    indexoption = "--large-index"
shellcall("bowtie2-build "+dbdirectory+dbname+"-tRNAgenome.fa "+dbdirectory+dbname+"-tRNAgenome "+indexoption+" -p "+str(cores)+"", failquit = True)




