#!/usr/bin/env python

import re
import os
import sys
import itertools

from collections import defaultdict
import argparse
from parsetrnas import *
from trnasequtils import *



    
#This program gets the mature tRNA sequences

def main(**args):
    args = defaultdict(lambda: None, args)
    #>Saccharomyces_cerevisiae_tRNA-Ala-AGC-1-2 (tRNAscan-SE ID: chrVI.trna6) chrVI:204924-204996 (-) Ala (AGC) 73 bp Sc: 69.6
                                #nmt-tRNA-Gln-TTG-13-1
    #trnanamere = re.compile(r"^\>\S+_(tRNA\-\w+\-\w+\-\d+\-\d+)\s+")#\((S+)\)")
    genomefile = os.path.expanduser(args["genome"])
    if "maturetrnafa" in args:
        maturetrnafa = open(args["maturetrnafa"], "w")
    else:
        maturetrnafa=sys.stdout
    gtrnatrans = None
    prokmode = args["prokmode"]
    if args["namemap"]:
        gtrnafa = open(args["namemap"], "r")
        gtrnatrans = dict()
        for currline in gtrnafa:
            fields = currline.rstrip().split()


            if "tRNAscan-SE_id" == fields[0] or len(fields) < 2 or len(fields[1].split("-")) < 4:
                continue
            shortname = fields[0].split("-")[0]

            gtrnatrans[shortname] = fields[1] 

        #sys.exit()
    elif args["gtrnafa"]:
        trnanamere = re.compile(r"^\>\S+_((?:\w+\-)?tRNA\-\w+\-[\w\?]+\-\d+\-\d+)\s+\((?:tRNAscan\-SE\s+ID:\s+)?(\S+)\)")

        gtrnafa = open(args["gtrnafa"], "r")
        gtrnatrans = dict()
        for currline in gtrnafa:
            trnamatch = trnanamere.match(currline)
            #print >>sys.stderr, currline
            if trnamatch:

                #print >>sys.stderr, trnamatch.group(1)+":"+trnamatch.group(2)
                gtrnatrans[trnamatch.group(2)] = trnamatch.group(1)
            elif currline.startswith(">"):
                pass
                #print >>sys.stderr, currline
        if len(gtrnatrans.keys()) == 0:
            print >>sys.stderr, "Could not extract names from gtrnadb fasta file"
            print >>sys.stderr, "must have names in format '>Saccharomyces_cerevisiae_tRNA-Ala-AGC-1-10 (tRNAscan-SE ID: chrXIII.trna9)'"
            sys.exit()
            
    
        
    alltrnas = list()
    trnascantrnas = list()
    trnadbtrnas = list()
    trnacentraltrnas = list()
    for currfile in args["trnascan"]:
        if gtrnatrans:
            trnadbtrnas.extend(readtRNAdb(currfile, genomefile, gtrnatrans))
        else:
            trnascantrnas.extend(readtRNAscan(currfile, genomefile))
            #print >>sys.stderr, len(trnascantrnas)
            
    #print >>sys.stderr, "||**"
    '''
    for currfile in args["trnascan"]:
        trnascantrnas.extend(readtRNAscan(currfile, args["genome"]))
    trnacentraltrnas = list()
    for currfile in args["rnacentral"]:
        trnacentraltrnas.extend(readrnacentral(currfile,args.chromtranslate,mode = 'transcript'))
    '''    
    extratrnas = list()
    if args["addtrna"]:
        #addtrnafa = open(args["addtrna"], "r")
        for currname, currseq in readmultifasta(args["addtrna"]):
            #I'm making sure that the name fits this tRNA convention I'm establishing
            namefields = currname.split("-")
            #print >>sys.stderr,currname
            #print >>sys.stderr,currseq
            if len(namefields) != 4: #or namefields[0] != "newtRNA":
                print >>sys.stderr, "additional tRNA "+currname+" from "+args["addtrna"] +" does not use a valid tRNA name"
                print >>sys.stderr, "valid named in the format 'newtRNA-xxx-NNN-1-1'"
                print >>sys.stderr, len(namefields) 
                sys.exit(1)
                pass
            
            
            extratrnas.append(tRNAtranscript(currseq.replace("U","T"), None, namefields[1], namefields[2], [], None, name = currname, artificialtrna = True))
    alltrnas = list(getuniquetRNAs(trnascantrnas)) + trnacentraltrnas + trnadbtrnas + extratrnas
    trnabed = None
    if args["bedfile"]:
        trnabed = open(args["bedfile"], "w")
    
    trnatable = None
    if args["maturetrnatable"]:
        trnatable = open(args["maturetrnatable"], "w")
    
    
    
    
    def readmultistk(struct):
        currrecord = ""
        structs = list()
        for line in struct.split("\n"):
            currrecord += line+"\n"
            if line == "//":
                yield currrecord
                currrecord = ""
                
    
    margin = 20
    anticodoncount = defaultdict(int)
    trnanames = dict()
    trnalist = list()
    for currtrans in alltrnas:
        if currtrans.name is None:
            name = 'tRNA-'+currtrans.amino + currtrans.anticodon+ str(anticodoncount[currtrans.anticodon]+ 1)
            currtrans.name = name
            anticodoncount[currtrans.anticodon] += 1
            
            
    scriptdir = os.path.dirname(os.path.realpath(sys.argv[0]))+"/"
    trnacmfile = args["cmmodel"]

        
    stkfile = args["trnaalignment"]
    if args["trnaalignment"]:
        devnull = open(os.devnull, 'w')
        seqfile = tempmultifasta(((currtrans.name, currtrans.getmatureseq(addcca = not prokmode)) for currtrans in alltrnas))
        cmcommand = ['cmalign', "-o", stkfile,"--nonbanded", "--notrunc","-g",trnacmfile,seqfile.name]
        
        cmrun = subprocess.Popen(cmcommand, stdout = devnull)
        result = cmrun.wait()
        if result:
            print >>sys.stderr, " ".join(cmcommand)
            print >>sys.stderr, "Failure to align tRNAs"

            sys.exit(1)
        seqfile.close()
        #print >>sys.stderr, " ".join(cmcommand)

        
    if args["locibed"]:
        locibed = open(args["locibed"],"w")
        
    for currtrans in alltrnas:
        name = currtrans.name
        #trnanames[name] = currtrans
        trnalist.append(name)
        #print >>sys.stderr, name
        print >>maturetrnafa, ">"+name
        matureseq = currtrans.getmatureseq(addcca = not prokmode)
        print >>maturetrnafa, str("N" * margin) +matureseq+str("N" * margin)
        locinames = ",".join(currlocus.name for currlocus in currtrans.loci)
        if locinames == "":
            locinames = "NA"
        if trnatable is not None:
            print >>trnatable, "\t".join([name,locinames,currtrans.amino,currtrans.anticodon])
        if trnabed is not None:
            transcriptrange = GenomeRange("genome", name, margin, margin + len(matureseq), strand = "+", name = name)
            print >>trnabed, transcriptrange.bedstring()
        if args["locibed"]:
            itemrgb = "0"
            for currlocus in currtrans.loci:
                pass
                trnalength = currlocus.loc.end - currlocus.loc.start
                if currlocus.intron is None:
                    print >>locibed, currlocus.loc.bedstring()+"\t"+"\t".join([str(currlocus.loc.start), str(currlocus.loc.end),itemrgb,"1",str(trnalength),"0"])
                else:
                    
                    blockcounts = 2
                    blocksizes = str(currlocus.intron[0] + 1)+","+str(trnalength - currlocus.intron[1] - 1)
                    blockstarts = "0,"+str(currlocus.intron[1] + 1)
                    print >>locibed, currlocus.loc.bedstring()+"\t"+"\t".join([str(currlocus.loc.start), str(currlocus.loc.end),itemrgb,str(blockcounts),str(blocksizes),str(blockstarts)])

    
    
    #sys.exit()
    trnamods = dict()
    allmods = set()
    nomatches = 0
    trnamismatches = dict()
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate fasta file containing mature tRNA sequences.')
    parser.add_argument('--trnascan',  nargs='+', default=list(),
                       help='tRNAscan-SE file')
    parser.add_argument('--rnacentral', nargs='+', default=list(),
                       help='RNAcentral tRNA file')
    parser.add_argument('--genome', 
                       help='fasta sequence of genome')
    parser.add_argument('--chromtranslate', 
                       help='translation table of chromosome names')
    parser.add_argument('--maturetrnafa', 
                       help='output file for mature tRNA sequences')
    parser.add_argument('--bedfile', 
                       help='Output bedfile of coordinates for mature tRNA')
    parser.add_argument('--tag', nargs='1',
                       help='tag to be added to tRNA name')
    parser.add_argument('--maturetrnatable',
                       help='Output table of mature tRNAs')
    parser.add_argument('--locibed',
                       help='Output BED format trna loci')
    parser.add_argument('--trnaalignment',
                       help='Output stockholm format alignment of mature tRNAs')
    parser.add_argument('--gtrnafa', 
                       help='fasta sequence tRNAs from tRNAscan-SE database')
    parser.add_argument('--cmmodel', 
                       help='covariance model to use')
    parser.add_argument('--prokmode', action="store_true", default=False,
                       help='skip adding CCA to mature tRNAs')
    
    args = vars(parser.parse_args())
    main(args)    

