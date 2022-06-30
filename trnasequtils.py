#!/usr/bin/env python

import pysam
import sys
import tempfile
import re
import gzip
import subprocess
import os.path
from distutils.spawn import find_executable

from collections import defaultdict

allaminos = ('Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly','His','Ile','Ile2','Leu','Lys','Met','iMet','fMet','Phe','Pro','Ser','Thr','Trp','Tyr','Val','SeC','Sup','Undet')

def readmultifasta(fafile):
    #print chrom+":"+ chromstart+"-"+ chromendre
    if fafile == "stdin":
        fafile = sys.stdin
    elif fafile.endswith(".gz"):
        fafile = gzip.open(fafile, "rb")
    else:
        fafile = open(fafile, "r")
    currloc = 0
    currseq = None
    sequence = ""
    reheader = re.compile(r"\>([^\s\,]+)")
    for line in fafile:
        line = line.rstrip("\n")
        currheader = reheader.match(line)
        if currheader and sequence != "" and currseq is not None:
            yield currseq, sequence
            currseq = currheader.groups(1)[0]
            sequence = ""
        elif currheader:
            currseq = currheader.groups(1)[0]
        else:
            sequence += line
    if currseq is not None:
        yield currseq, sequence

#returns location of the genome file
def readmultifastq(fqfile, getnameline = False):
    #print chrom+":"+ chromstart+"-"+ chromend
    if fqfile == "stdin":
        fqfile = sys.stdin
    elif fqfile.endswith(".gz"):
        fqfile = gzip.open(fqfile, "rb")
    else:
        fqfile = open(fqfile, "r")
    currloc = 0
    currseq = None
    sequence = ""
    quality = ""
    if getnameline:
        reheader = re.compile(r"\@(.+)")
    else:
        reheader = re.compile(r"\@([^\s\,]+)")
    qualheader = re.compile(r"\+([^\s\,]*)")
    readqual = False
    for line in fqfile:
        #print line
        line = line.rstrip("\n")
        seqheader = reheader.match(line)
        qheader = qualheader.match(line)
        if readqual:
            quality += line
            if len(quality) == len(sequence):
                yield currseq, sequence, quality
                readqual = False
        elif seqheader:
            currseq = seqheader.groups(1)[0]
            sequence = ""
            quality = ""
        elif qheader and readqual == False:
            readqual = True
            pass
        else:
            sequence += line
                    
def fastadict(fafile):
    seqdict = dict()
    for name, seq in readmultifasta(fafile):
        seqdict[name] = seq
    return seqdict
        
def tempmultifasta(allseqs):
    fafile = tempfile.NamedTemporaryFile(suffix=".fa")


    for seqname, seq in allseqs:
        fafile.write(">"+seqname+"\n")
        fafile.write(seq+"\n")
    fafile.flush()
    return fafile
def invertstrand(strand):
    if strand == "+":
        return "-"
    elif strand == "-":
        return "+"

class alignment:
    def __init__(self, sequences):
        if len(sequences.keys()) < 1:
            raise EmptyAlignException() 
        self.aligns = sequences

        if max(len(curr) for curr in self.aligns.itervalues()) != min(len(curr) for curr in self.aligns.itervalues()):
            print >>sys.stderr, "Non-matching sequence lengths in multiple alignment"
        #self.alignlength = len(self.aligns[self.aligns.keys()[0]])
        self.alignlength = max(len(curr) for curr in self.aligns.itervalues())
    def getseqlength(self, seqname):
        return len(self.aligns[seqname].replace("-",""))
    def numseqs(self):
        return len(self.aligns.keys())
    def getseqnames(self):
        return self.aligns.keys()
    def getalignseq(self ,seqname):
        return self.aligns[seqname]
    def toseqdict(self):
        seqs = dict()
        for currname, currseq in self.aligns.iteritems():
            seqs[currname] = string.translate(currseq, None, "-.~")
        return seqs
    def getsubset(self, subset):
        newaligns = dict()
        for curr in subset:
            newaligns[curr] = ""
        for i in range(self.alignlength):
            if any(self.aligns[curr][i] not in gapchars for curr in subset):
                for curr in subset:
                    newaligns[curr] += self.aligns[curr][i]
        return alignment(newaligns)
    def prottocodonalignment(self, nucseqs):
        #I don't do any checking to make sure that the nucleotide sequence is equivalent here
        newalign = defaultdict(str)
        currpoint = defaultdict(int)
        for i in range(self.alignlength):
            for currgene in nucseqs.iterkeys():
                if self.aligns[currgene][i] != "-":
                    newalign[currgene] += nucseqs[currgene][currpoint[currgene]:currpoint[currgene] + 3]
                    currpoint[currgene] += 3
                else:
                    newalign[currgene] += "---"
         
        #print "|".join(str(curr) for curr in currpoint.values())
        #print "|".join(str(len(curr)) for curr in nucseqs.values())
        #print "|".join(curr for curr in newalign.values())
        
        #print str(self.alignlength) pyroIsla1  
        return alignment(newalign, aligned = True)
    def tempmultifasta(self):
        return tempmultifasta(self.aligns)
    def fastaformat(self, prefix = ""):
        output = ""
        for seqname, seq in self.aligns.iteritems():
            output += ">"+prefix+seqname+"\n"+seq+"\n"
        return output
    def getsegment(self,start, end):
        newalign = dict()
        for name, seq in self.aligns.iteritems():
            newalign[name] = seq[start:end]
        return alignment(newalign)
    def getsubsets(self,windowsize, stepsize):
        for curr in range(0, self.alignlength - windowsize, stepsize):
            start = curr
            end = curr + windowsize
            yield self.getsegment(start,end)
    def removeemptyseqs(self):
        newalign = dict()
        for name, seq in self.aligns.iteritems():
            if len(string.translate(seq, None, "-.")) != 0:
                newalign[name] = seq
        return alignment(newalign)
    #http://bioinformatics.oxfordjournals.org/content/25/5/668.full

    def phylipformat(self):
        #print ",".join(str(curr)for curr in self.aligns.keys())
        output = str(len(self.aligns.keys()))+ " "+str(self.alignlength)+"\n"
        for currgene in self.aligns.iterkeys():
            #sys.stderr.write("**"+str(currgene).ljust( 11)+"**\n")
            output += str(currgene).ljust( 14)
            output +=  self.aligns[currgene]+"\n"
        #print output
        return output


    def nexusformat(self):
        output = "#NEXUS\nBegin data;\nDimensions ntax="+str(len(self.aligns.keys()))+" nchar="+str(self.alignlength)+";\n"
        output += "Format datatype=dna symbols=\""+"ATCG"+"\" missing=? gap=-;\n"
        output += "Matrix\n"
        for currgene in self.aligns.iterkeys():
            output += currgene.ljust( 11)
            output +=  self.aligns[currgene]+"\n"
        output+=";\nEnd;"    
    def getrealcoordinate(self, alignname, coord):
        currreal = 0
        #print self.aligns[alignname]
        #print coord
        for i in range(0, coord):
            if self.aligns[alignname][i] not in  set("-.~"):
                currreal += 1
        return min([currreal, len(list(curr for curr in self.aligns[alignname] if curr not in set("-.~"))) - 1])
    def getseqrange(self, alignname, start, end):
        return self.getsegment(self.getaligncoordinate(alignname, start),self.getaligncoordinate(alignname, end))

            
    def printstk(self, name = None):
        print "# STOCKHOLM 1.0"
        if name is not None:
            print "#=GF ID "+name
        padlength = max(len(curr) for curr in self.aligns.iterkeys()) + 4
        for currname, currseq in self.aligns.iteritems():
            print string.ljust(currname,padlength ) +currseq
        print "//"
    def printhtml(self, name = None):
        print "<CODE>"
        for currname, currseq in self.aligns.iteritems():
            print currname + "\t"+currseq +"<BR/>"
        print "</CODE>"
    def clustalwformat(self):
        output = "CLUSTAL W 2.1 multiple sequence alignment\n\n"
        conservestring = ""
        for i in range(0, self.alignlength):
            conservestring += ":"                     
        for currpos in range(0, self.alignlength, 60):
            for seqname, seq in self.aligns.iteritems():
                output += seqname+"\t"+seq[currpos:min([currpos+60,self.alignlength])]+"\n"
            output += "\t"+conservestring[currpos:min([currpos+60,self.alignlength])] + "\n\n"
            
        return output
        

class RnaAlignment(alignment):
    def __init__(self,alignseqs, structure, consensus = None,energies = None):
        self.aligns = alignseqs
        self.currstruct = structure
        self.energies = energies 
        self.consensus = consensus
        self.alignlength = max(len(curr) for curr in alignseqs.values())
    def addupstream(self, seqs, struct = None):
        newseqs = dict()
        #print >>sys.stderr, seqs.keys()
        for curr in self.aligns.iterkeys():
            newseqs[curr] = seqs[curr] + self.aligns[curr]
        if struct is None:
            newstruct = (max(len(curr) for curr in seqs.itervalues()) * ":") + self.currstruct
        else:
            newstruct = struct + self.currstruct
        return RnaAlignment(newseqs, newstruct)
    def adddownstream(self, seqs, struct = None):
        newseqs = dict()
        for curr in self.aligns.iterkeys():
            newseqs[curr] =  self.aligns[curr]+  seqs[curr]
        if struct is None:
            newstruct = self.currstruct + (max(len(curr) for curr in seqs.itervalues()) * ":")
        else:
            newstruct = self.currstruct + struct
        return RnaAlignment(newseqs, newstruct)
    def addmargin(self, length):
        newseqs = dict()
        for curr in self.aligns.iterkeys():
            newseqs[curr] =  length*"N" + self.aligns[curr]+ length*"N"
        newstruct = length * ":" + self.currstruct + length * ":"
        if self.consensus is None:
            newconsensus = None
        else:
            newconsensus  =   length*"N" + self.consensus+ length*"N"
        return RnaAlignment(newseqs, newstruct, consensus = newconsensus)

    def viennaformat(self):
        output = ""
        for currseq in self.aligns.iterkeys():
            output += ">"+currseq+"\n"
            output += self.aligns[currseq]+"\n"
            output += self.currstruct+"\n"
        return output
    def viennatempfile(self):
        viennafile = tempfile.NamedTemporaryFile()
        viennafile.write(self.viennaformat())
        viennafile.flush()
        return viennafile
    def printstk(self, name = None):
        print "# STOCKHOLM 1.0"
        if name is not None:
            print "#=GF ID "+name
        for currname, currseq in self.aligns.iteritems():
            print currname + "\t"+currseq
        structline = ""
        structpos = 0
        secpos = 0
        print "#=GC SS_cons\t"+self.currstruct
        print "//"
    def printhtml(self, name = None):
        for currname, currseq in self.aligns.iteritems():
            print currname + "\t"+currseq +"</BR>"
        structline = ""
        structpos = 0
        secpos = 0
        print "#=GC SS_cons\t"+self.currstruct+"</BR>"
        print "//"+"</BR>"
def convertmaturealign(rnaalign):
    newseqs = dict()
    newalign = rnaalign.struct
    for name, seq in rnaalign.aligns.iteritems():
        newseqs[name]
    
def readrnastk(stk):
    seqs = defaultdict(str)
    struct = ""
    consensus = ""
    energyscore = None
    for line in stk:
        line = line.rstrip()
        if line.startswith("//"):
            #print >>sys.stderr, "||||||||||||||||||"
            #print >>sys.stderr, consensus
            if consensus == "":
                consensus = None
            #print >>sys.stderr, consensus
            yield RnaAlignment(seqs, struct, consensus = consensus)
            seqs = defaultdict(str)
            struct = ""
            consensus = ""
            energyscore = None
        elif not line.startswith("#") and len(line.split()) > 1:
            currname = line.split()[0]
            currseq = line.split()[1]
            seqs[currname] += currseq
        elif line.startswith("#=GC SS_cons"):
            struct += line.split()[2]
        elif line.startswith("#=GC RF"):
            #print >>sys.stderr, "********************"
            #
            #print >>sys.stderr, line.split()[2]
            #print >>sys.stderr, "#############################"
            consensus += line.split()[2]
            
def uniqueorder(inp):
    
    alldata = set()
    for curr in inp:
        if curr not in alldata:
            yield curr
            alldata.add(curr)
class transcriptfile:
    def __init__(self, trnafilename):
        trnafile = open(trnafilename)
        locustranscript = dict()
        trnatranscripts = list()
        loci = list()
        amino = dict()
        anticodon = dict()
        transcriptdict = defaultdict(set)
        aminoorder = list()
        anticodonorder = list()
        for i, line in enumerate(trnafile):
            fields = line.split()
            if len(fields) < 2:
                continue
            trnatranscripts.append(fields[0])
            amino[fields[0]] = fields[2]
            anticodon[fields[0]] = fields[3]
            aminoorder.append(fields[2])
            anticodonorder.append(fields[3])
            for currlocus in fields[1].split(','):
                locustranscript[currlocus] = fields[0]
                loci.append(currlocus)
                transcriptdict[fields[0]].add(currlocus)

        
        self.locustranscript = locustranscript
        self.transcripts = trnatranscripts
        self.amino = amino
        self.anticodon = anticodon
        self.transcriptdict = transcriptdict
        self.loci = loci
        
        self.aminoorder = tuple(uniqueorder(aminoorder))             
        self.anticodonorder = tuple(uniqueorder(anticodonorder))
    def gettranscripts(self):
        return set(self.transcripts)
    def getlocustranscript(self, locus):
        return  self.locustranscript[locus]
    def getloci(self):
        return  self.loci
    def getamino(self, trna):
        return  self.amino[trna]
    def getanticodon(self, trna):
        return  self.anticodon[trna]
        
    def allaminos(self):
        return  self.aminoorder
    def allanticodons(self):
        return  self.anticodonorder
    def getaminotranscripts(self, trnaamino):
        return  set(curr for curr in self.transcripts if trnaamino == self.amino[curr])
    def getanticodontranscripts(self, trnaanticodon):
        return   set(curr for curr in self.transcripts if trnaanticodon == self.anticodon[curr])
        
        
def getpairfile(pairfilename):
    pairfile = open(pairfilename)
    for currline in pairfile:
        
        fields = currline.split()
        if len(fields) > 1:
            yield fields[0], fields[1]
    
class extraseqfile:
    def __init__(self, extraseqfilename):
        try:
            extrafile = open(extraseqfilename)
            seqlist = list()
            seqfasta = dict()
            seqbed = dict()
            directory = os.path.dirname(extraseqfilename)
            for i, line in enumerate(extrafile):
                fields = line.split()
                if len(fields) < 2:
                    continue
                seqfasta[fields[0]] = directory+"/"+fields[1]
                seqbed[fields[0]] = directory+"/"+fields[2]
                
                seqlist.append(fields[0])

            
            #bamlist = list(curr + "_sort.bam" for curr in samplefiles.iterkeys())
            self.seqlist = seqlist
            self.seqfasta = seqfasta
            self.seqbed = seqbed
            #self.bamlist = list(curr+ "_sort.bam" for curr in samplelist)
        except IOError as e:
            self.seqlist = list()
            self.seqfasta = dict()
            self.seqbed = dict()
            #print >>sys.stderr, "**||**" + extraseqfilename
            #print >>sys.stderr,"extraseqfile I/O error({0}): {1}".format(e.errno, e.strerror)
    def getseqnames(self):
        seqnamedict = defaultdict(set)
        for curr in seqlist:
            seqnamedict[currseq] += set(curr.name for curr in readbed(self.seqbed[currseq]))
        return seqnamedict
    def getseqbeds(self):
        return self.seqbed
        
            
class samplefile:
    def __init__(self, samplefilename, bamdir = "./"):
        try:
            samplefile = open(samplefilename)
            samplelist = list()
            samplefiles = dict()
            replicatename = dict()
            
            
            replicatelist = list()
            for i, line in enumerate(samplefile):
                fields = line.split()
                if len(fields) < 2:
                    continue
                samplefiles[fields[0]] = fields[2]
                replicatename[fields[0]] = fields[1]
                
                samplelist.append(fields[0])
                if fields[1] not in set(replicatelist):
                    replicatelist.append(fields[1])
            
            #bamlist = list(curr + "_sort.bam" for curr in samplefiles.iterkeys())
            samplenames = list(curr  for curr in samplefiles.iterkeys())
            if bamdir is None:
                bamdir = "./"
            self.bamdir = bamdir
            self.samplelist = samplelist
            self.samplefiles = samplefiles
            self.replicatename = replicatename
            self.replicatelist = replicatelist
            #self.bamlist = list(curr+ "_sort.bam" for curr in samplelist)
        except IOError:
            print >>sys.stderr, "Cannot read sample file "+samplefilename
            print >>sys.stderr, "exiting..."
            sys.exit(1)
    def getsamples(self):
        return self.samplelist
    def getbamlist(self):
        return list(curr+ ".bam" for curr in self.samplelist)
    def getbam(self, sample):
        return self.bamdir + "/" + sample + ".bam" 
    def getmergebam(self, sample):
        return self.bamdir + "/" + sample + "-merge.bam" 
    def getfastq(self, sample):
        return self.samplefiles[sample]
    def getreplicatename(self, sample):
        return self.replicatename[sample]
    def allreplicates(self):
        return self.replicatelist
    def getrepsamples(self, replicate):
        return list(currsample for currsample in self.samplelist if self.replicatename[currsample] == replicate)
    def getfastqsample(self, fastqname):
        for currsample in self.samplefiles.iterkeys():
            if self.samplefiles[currsample] == fastqname:
                return currsample
        
def getsizefactors( sizefactorfilename):
    sizefactorfile = None
    try:
        sizefactorfile = open(sizefactorfilename)
    except IOError as err:
        print >>sys.stderr, "Cannot read size factor file "+sizefactorfilename
        print >>sys.stderr, "check Rlog.txt"
        sys.exit(1)
    sizefactors = dict()
    bamheaders =list()
    sizes = list()
    for i, line in enumerate(sizefactorfile):
        if i == 0:
            bamheaders = list(curr.strip("\"\n") for curr in line.split())
        elif i == 1:
            sizes = list(float(curr.strip("\"\n")) for curr in line.split())
    for i in range(0, len(bamheaders)):
        sizefactors[bamheaders[i]] = sizes[i]
        #print >>sys.stderr, bamheaders[i]+":"+ str(sizes[i])
    return sizefactors

#special class that uses the read indentifier for hashing in sets

class GenomeRange:
    __slots__ = "dbname", "chrom", "strand","name", "fastafile", "start", "end", "data"
    def __eq__(self, other):
        return self.strand == other.strand and self.chrom == other.chrom and self.start == other.start and self.end == other.end
    def __hash__(self):
        return  self.start + self.end + hash(self.chrom) + hash(self.strand)
    def __init__(self, dbname, chrom, start, end, strand = None,name = None, orderstrand = False, data = None, fastafile = None):
        self.dbname =dbname
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.fastafile = fastafile
        if orderstrand and self.start > self.end:
            temp = self.start
            self.start = self.end
            self.end = temp
            self.strand = invertstrand(strand)
        self.data = data
        self.name = name        
    def coverage(self, other):
        if self.strand == other.strand and self.chrom == other.chrom:
            start = max([self.start,other.start])
            end = min([self.end,other.end])
            
            if end - start < 0:
                return 0
            else:
                return end - start
        else:
            return 0
            
    def bedstring(self, name = None,score = 1000):
        if self.strand == None:
            self.strand = "+"
        if name is None and self.name is None:
            name = "FEAT"
        elif name is None:
            name = self.name
        return "\t".join([self.chrom,str(self.start),str(self.end),name,str(score), self.strand])
    def length(self):
        return self.end - self.start
    def addmargin(self, dist = 50, name = None):
        newname = name
        if name is None:
            newname = self.name
        return GenomeRange(self.dbname, self.chrom, self.start - dist,self.end + dist,self.strand, name = newname, fastafile = self.fastafile)
    def getupstream(self, dist = 50, name = None):
        newname = name
        if name is None:
            newname = self.name
        return GenomeRange(self.dbname, self.chrom, self.start - dist,self.start,self.strand, name = newname, fastafile = self.fastafile)
    def getdownstream(self, dist = 50, name = None):
        newname = name
        if name is None:
            newname = self.name
        if self.strand != '-':
            return GenomeRange(self.dbname, self.chrom, self.end,self.end + dist,self.strand, name = newname, fastafile = self.fastafile)
        else:
            return GenomeRange(self.dbname, self.chrom, self.start - dist,self.start,self.strand, name = newname, fastafile = self.fastafile)
    def getfirst(self, dist = 50, name = None):
        newname = name
        if name is None:
            newname = self.name
        if self.strand == "-":
            return GenomeRange(self.dbname, self.chrom, self.end - dist,self.end ,self.strand, name = newname, fastafile = self.fastafile)
        else:
            return GenomeRange(self.dbname, self.chrom, self.start,self.start + dist,self.strand, name = newname, fastafile = self.fastafile)
    def getlast(self, dist = 50, name = None):
        newname = name
        if name is None:
            newname = self.name
        if self.strand == "-":
            return GenomeRange(self.dbname, self.chrom, self.start,self.start + dist,self.strand, name = newname, fastafile = self.fastafile)
        else:
            return GenomeRange(self.dbname, self.chrom, self.end - dist,self.end ,self.strand, name = newname, fastafile = self.fastafile)
    def getbase(self, basenum, name = None):
        newname = name
        if name is None:
            newname = self.name
        if self.strand == "-":
            return GenomeRange(self.dbname, self.chrom, self.end - basenum - 1 ,self.end -basenum,self.strand, name = newname, fastafile = self.fastafile)
        else:
            return GenomeRange(self.dbname, self.chrom, self.start +basenum ,self.start + basenum + 1,self.strand, name = newname, fastafile = self.fastafile)
    def antisense(self):
        if self.strand == "+":
            newstrand = "-"
        else:
            newstrand = "+"
        return GenomeRange(self.dbname, self.chrom, self.start,self.end,newstrand, name = self.name, fastafile = self.fastafile)
    def bamseq(self):
        if self.strand == "+":
            return self.data["seq"]
        else:
            return revcom(self.data["seq"])
    def getgc(self):
        seq = self.bamseq()
        
        return sum(1 if currbase in set(["G","C"]) else 0 for currbase in seq) 
class GenomeRead(GenomeRange):
    def __eq__(self, other):
        return self.name == other.name
    def __hash__(self):
        return  hash(self.name)
    def __init__(*args, **nargs):
        GenomeRange.__init__(*args, **nargs)

'''
Still need to add trailer fragment code, both here and elsewhere
'''
def getfragtype(currfeat, currread, maxoffset = 10):
    if currread.start < currfeat.start + maxoffset and currread.end > currfeat.end - maxoffset:
        return "Whole"
    elif currread.start < currfeat.start + maxoffset:
        if currfeat.strand == "+":
            return "Fiveprime"
        else:
            return "Threeprime"
    elif currread.end > currfeat.end - maxoffset:
        if currfeat.strand == "+":
            return "Threeprime"
        else:
            return "Fiveprime"
def getendtype(currfeat, currread, maxoffset = 10):
    endtype = None
    if currread.end == currfeat.end:
        endtype = "CCA"
    elif currread.end == currfeat.end -1:
        endtype = "CC"
    elif currread.end == currfeat.end -2:
        endtype = "C"
    elif currread.end == currfeat.end -3:
        endtype = ""
    return endtype
smallrnatypes = set([])            
def readfeatures(filename, orgdb="genome", seqfile= None, removepseudo = False):
    if filename.endswith(".bed") or filename.endswith(".bed.gz"):
        return readbed(filename, orgdb, seqfile)
    elif filename.endswith(".gtf") or filename.endswith(".gtf.gz") or filename.endswith(".gff") or filename.endswith(".gff.gz"):
        #print >>sys.stderr, removepseudo
        return (curr for curr in readgtf(filename, orgdb, seqfile, filterpsuedo = removepseudo, filtertypes =set(['retained_intron','antisense','lincRNA']) ))
    else:
        print >>sys.stderr, filename+" not valid feature file"
        sys.exit()


def readgtf(filename, orgdb="genome", seqfile= None, filterpsuedo = False, replacename = False, filtertypes = set(['retained_intron','antisense','lincRNA'])):
    bedfile = None
    #print >>sys.stderr, "****"
    if filename == "stdin":
        bedfile = sys.stdin
    elif filename.endswith(".gz"):
        bedfile = gzip.open(filename, 'rb')
    else:
        bedfile = open(filename, "r")
    
    for currline in bedfile:
        #print currline
        if currline.startswith('track') or currline.startswith('#'):
            continue
        fields = currline.rstrip().split("\t")
        if len(fields) > 2:
            biotype = "Unknown"
            featname = None
            genename = None
            #print >>sys.stderr, len(fields)
            genesource = fields[1]  
            #retained introns are often other things as well, so I skip em
            if fields[2] != "transcript" or genesource in filtertypes:
                continue

              
            for currattr in fields[8].rstrip(";").split(";"):
                #print >>sys.stderr,  currattr
                currname = currattr.strip().split()[0]
                currvalue = currattr.strip().split()[1]
                if currname == "transcript_name":
                    featname = currvalue.strip('"')
                elif currname == "name" or currname == "gene_id" and featname is None:
                    featname = currvalue.strip('"')

                elif currname == "gene_biotype":
                    biotype = currvalue.strip('"')
                elif currname == "gene_name":
                    genename = currvalue.strip('"')
                #print >>sys.stderr, "***||"
            if genename is None:
                #print >>sys.stderr, "No name for gtf entry "+featname
                genename = featname
                pass
            if filterpsuedo and biotype == "pseudogene":
                #print >>sys.stderr, "*******"
                continue
                
            if genesource == 'ensembl':
                #print >>sys.stderr, biotype
                genesource = biotype
            if not (fields[6] == "+" or fields[6] == "-"):
                print >>sys.stderr, "strand error in "+filename
                skippedlines += 1
            elif not (fields[3].isdigit() and fields[4].isdigit()):
                print >>sys.stderr, "non-number coordinates in "+filename
                print >>sys.stderr, currline
                
                skippedlines += 1
            else:
                if replacename:
                    featname = genename
                yield GenomeRange( orgdb, fields[0],fields[3],fields[4],fields[6], name = featname, fastafile = seqfile, data = {"biotype":biotype, "source":genesource, "genename":genename,"feature":fields[2]})
            
def readbed(filename, orgdb="genome", seqfile= None, includeintrons = False):
    bedfile = None
    if filename == "stdin":
        bedfile = sys.stdin
    elif filename.endswith(".gz"):
        bedfile = gzip.open(filename, 'rb')
    else:
        bedfile = open(filename, "r")
    skippedlines = 0
    
    for currline in bedfile:
        #print currline
        data = dict()
        if currline.startswith('track') or currline.startswith('#'):
            continue
        fields = currline.rstrip().split()
        if len(fields) > 2:
            if len(fields) < 5:
                strand = "+"
            else:
                strand = fields[5]
            if not (strand == "+" or strand == "-"):
                print >>sys.stderr, "strand error in "+filename
                skippedlines += 1
            elif not (fields[1].isdigit() and fields[2].isdigit()):
                print >>sys.stderr, "non-number coordinates in "+filename
                print >>sys.stderr, currline
                skippedlines += 1
            else:
                if includeintrons and len(fields) > 7:
                    data["blockcount"] = int(fields[9])
                    data["blocksizes"] = tuple(int(curr) for curr in fields[10].rstrip(",").split(","))
                    data["blockstarts"] = tuple(int(curr) for curr in fields[11].rstrip(",").split(",")) 
                yield GenomeRange( orgdb, fields[0],fields[1],fields[2],strand, name = fields[3], fastafile = seqfile, data = data)
    
    if skippedlines > 0:
        print >>sys.stderr, "skipped "+str(skippedlines)+" in "+filename
def ifelse(arg, trueres,falseres):
    if arg:
        return trueres
    else:
        return falseres
        
        
def isprimarymapping(mapping):
    return not (mapping.flag & 0x0100 > 0)        
def issinglemapping(mapping):
    return mapping.mapq > 2
    
def mismatchnum(mapping):
    return int(mapping.get_tag("XM"))
    
def getbamrangeshortseq(bamfile, chromrange = None, primaryonly = False, singleonly = False, maxmismatches = None, allowindels=True, skiptags = False):

    try:
        if chromrange is not None:
            bamiter = bamfile.fetch(chromrange.chrom, chromrange.start, chromrange.end)
        else:
            bamiter = bamfile.fetch()
   
        for currline in bamiter: 
            if primaryonly and not isprimarymapping(currline):
                continue

            if singleonly and not issinglemapping(currline):
                continue
            if not allowindels and len(currline.cigar) > 1:
                continue
            if maxmismatches is not None and mismatchnum(currline) > maxmismatches:
                continue
            rname = bamfile.getrname(currline.rname)
            strand = "+"
            strand = ifelse(currline.is_reverse, '-','+')
            yield GenomeRead( "genome",rname,currline.pos,currline.aend,strand, name = currline.qname, data = {"seq":currline.seq} )

    except ValueError as err:
        #print>>sys.stderr, err
        #print>>sys.stderr, bamfile.name
        if chromrange is not None:
            #print >>sys.stderr, chromrange.chrom+":"+ str(chromrange.start)+"-"+str(chromrange.end) +" failed"
            pass
def getbamrangeshort(bamfile, chromrange = None, primaryonly = False, singleonly = False, maxmismatches = None, allowindels=True, skiptags = False):
    bamiter = None
    bamiter = None

    try:
        if chromrange is not None:
            bamiter = bamfile.fetch(chromrange.chrom, chromrange.start, chromrange.end)
        else:
            bamiter = bamfile.fetch()
   
        for currline in bamiter:
            #if mismatchnum(currline) > maxmismatches:
            #    print >>sys.stderr, "*|*"
            #    print >>sys.stderr, maxmismatches
            #    print >>sys.stderr, mismatchnum(currline)
            if primaryonly and not isprimarymapping(currline):
                continue

            if singleonly and not issinglemapping(currline):
                continue
            if not allowindels and len(currline.cigar) > 1:
                continue
            if maxmismatches is not None and mismatchnum(currline) > maxmismatches:
                #print >>sys.stderr, "skipped"
                continue
            rname = bamfile.getrname(currline.rname)
            strand = "+"
            strand = ifelse(currline.is_reverse, '-','+')
            
            yield GenomeRead( "genome",rname,currline.pos,currline.aend,strand, name = currline.qname)

    except ValueError as err:
        #print>>sys.stderr, err
        #print>>sys.stderr, bamfile.name
        if chromrange is not None:
            #print >>sys.stderr, chromrange.chrom+":"+ str(chromrange.start)+"-"+str(chromrange.end) +" failed"
            pass
        

        
def getbamrange(bamfile, chromrange = None, primaryonly = False, singleonly = False, maxmismatches = None, allowindels = True, skiptags = False):
    bamiter = None
    try:
        if chromrange is not None:
            bamiter = bamfile.fetch(chromrange.chrom, chromrange.start, chromrange.end)
        else:
            bamiter = bamfile.fetch()
   
        for currline in bamiter:
            
            #if mismatchnum(currline) > maxmismatches:
            #    print >>sys.stderr, "**"
            #    print >>sys.stderr, maxmismatches
            #    print >>sys.stderr, mismatchnum(currline)
            if primaryonly and not isprimarymapping(currline):
                continue

            if singleonly and not issinglemapping(currline):
                continue
            if maxmismatches is not None and mismatchnum(currline) > maxmismatches:
                #print >>sys.stderr, "skipped"
                continue
            rname = bamfile.getrname(currline.rname)
            strand = "+"
            strand = ifelse(currline.is_reverse, '-','+')
            #yield GenomeRead( "genome",rname,currline.pos,currline.aend,strand, name = currline.qname)
            #continue
            #print rname
            #need to fix this with cigar stuff
            #len(currline.pos)

            #not giving the reverse complement for now
            seq = currline.seq
            #print currline.cigar
            #print >>sys.stderr, currline.qname
            
            #[("YA",len(anticodons))] + [("YM",len(aminos))]  + [("YR",len(trnamappings))]
            uniqueac = True
            uniqueamino = True
            uniquetrna = True
            #print >>sys.stderr, dir(currline)
            mismatches = None
            alignscore = None
            secondbestscore = None
            uniquemapping = False

            if not skiptags:
                if currline.has_tag("YA") and currline.get_tag("YA") > 1:
                    uniqueac = False
                if currline.has_tag("YM") and currline.get_tag("YM") > 1:
                    uniqueamino = False
                if currline.has_tag("YR") and currline.get_tag("YR") > 1:
                    uniquetrna = False
                if currline.has_tag("XM"):
                    mismatches = currline.get_tag("XM")
                if currline.has_tag("XS"):
                    secondbestscore = float(currline.get_tag("XS"))
                if currline.has_tag("AS"):
                    alignscore = float(currline.get_tag("AS"))

                

                    
            if secondbestscore is None or alignscore > secondbestscore:
                uniquemapping = True
            #continue
            if not allowindels and len(currline.cigar) > 1:
                continue
                
                #"flags": currline.flag, 
            yield GenomeRead( "genome",rname,currline.pos,currline.aend,strand, name = currline.qname , data = {"score":currline.mapq, "CIGAR":currline.cigar,"CIGARstring":currline.cigarstring, "seq":seq, "flags": currline.flag, "qual":currline.qual,"bamline":currline,'uniqueac':uniqueac,"uniqueamino":uniqueamino,"uniquetrna":uniquetrna,"uniquemapping":uniquemapping})
    except ValueError as err:
        #print>>sys.stderr, err
        #print>>sys.stderr, bamfile.name
        if chromrange is not None:
            #print >>sys.stderr, chromrange.chrom+":"+ str(chromrange.start)+"-"+str(chromrange.end) +" failed"
            pass
        
class BamRead(GenomeRange):
    __slots__ = "bamline"
    def __init__(self,bamline, dbname, chrom, start, end, strand = None,name = None, fastafile = None):
        GenomeRange.__init__(self, dbname, chrom, start, end, strand = strand,name = name, fastafile = None)
        self.bamline = bamline
    def getmismatches(self):
        if self.bamline.has_tag("XM"):
            return int(self.bamline.get_tag("XM"))
        else:
            return None
    def isuniqueaminomapping(self):
        if self.bamline.has_tag("YM"):
            return not self.bamline.get_tag("YM") > 1
        else:
            return None
    def isuniqueacmapping(self):
        if self.bamline.has_tag("YA"):
            return not self.bamline.get_tag("YA") > 1
        else:
            return None
    def isuniquetrnamapping(self):
        if self.bamline.has_tag("YR"):
            return not self.bamline.get_tag("YR") > 1
        else:
            return None

    def hasindel(self):
        return len(self.bamline.cigar) > 1
    def getlength(self):
        return len(self.bamline.seq)
    def getseq(self):
        if self.strand == '-':
            return revcom(self.bamline.seq)
        else:
            return self.bamline.seq
            
    def getgc(self):
        seq = self.getseq()
        
        return sum(1 if currbase in set(["G","C"]) else 0 for currbase in seq) 
    def issinglemapped(self):
        return self.bamline.mapq >= 2
    def getcigar(self):
        return self.bamline.cigar

def getbam(bamfile, chromrange = None, primaryonly = False, singleonly = False, allowindels=True):
    bamiter = None
    try:
        if chromrange is not None:
            bamiter = bamfile.fetch(chromrange.chrom, chromrange.start, chromrange.end)
        else:
            bamiter = bamfile.fetch()
   
        for currline in bamiter:
            if primaryonly and not isprimarymapping(currline):
                continue

            if singleonly and not issinglemapping(currline):
                continue
           
            rname = bamfile.getrname(currline.rname)
            strand = "+"
            strand = ifelse(currline.is_reverse, '-','+')

            
            #[("YA",len(anticodons))] + [("YM",len(aminos))]  + [("YR",len(trnamappings))]
            
            #continue
            if currline.aend is None:
                continue
                print >>sys.stderr, currline
            if not allowindels and len(currline.cigar) > 1:
                continue
            yield BamRead( currline,"genome",rname,currline.pos,currline.aend,strand, name = currline.qname)
    except ValueError as err:
        #print>>sys.stderr, err
        #print>>sys.stderr, bamfile.name
        if chromrange is not None:
            #print >>sys.stderr, chromrange.chrom+":"+ str(chromrange.start)+"-"+str(chromrange.end) +" failed"
            pass     
#'uniqueac':uniqueac,"uniqueamino":uniqueamino,"uniquetrna":uniquetrna})
def isuniquetrnamapping(read):
    return read.data["uniquetrna"]
def isuniqueaminomapping(read):
    return read.data["uniqueamino"]
def isuniqueacmapping(read):
    return read.data["uniqueac"]
def issinglemapped(read):
    return (read.data["score"] >= 2)  
    
def getpileuprange(bamfile, chromrange = None):
    bamiter = None
    if chromrange is not None:
        bamiter = bamfile.pileup(chromrange.chrom, chromrange.start, chromrange.end)
    else:
        bamiter = bamfile.pileup()


    posdict = defaultdict(lambda: defaultdict(int))
    refdict = defaultdict(lambda:'N') 

    for currpos in bamiter:

        
        readcounts = defaultdict(int)
        reference = 'N'
        for read in currpos.pileups:
            aln = read.alignment
            
            
            if not read.indel and not read.is_del:
                readcounts[aln.seq[read.query_position]] += 1
        #print >>sys.stderr,currpos
        #print >>sys.stderr, currpos
        #print >>sys.stderr, "****"
        #yield currpos.pos, read.query_position, readcounts
        #yield currpos.pos, readcounts
        posdict[currpos.pos] = readcounts
        refdict[currpos.pos] = read.query_position
        
    for i in range(chromrange.start, chromrange.end):
        yield i, refdict[i], posdict[i]
        
def getseqdict(genelist, faifiles = None):
    namedict = getnamedict(genelist)
    allorgs = set(currgene.dbname for currgene in genelist)
    dbdict = dict()
    fastafiles = dict()
    for currorg in allorgs:
        dbdict[currorg] = dict()
        
    
    for currgene in genelist:
        dbdict[currgene.dbname][currgene.name] = currgene
        if currgene.fastafile is not None:
            fastafiles[currgene.dbname] =  currgene.fastafile
            
        else:
            #print >>sys.stderr, currgene.dbname+":"+genomefile(currorg)
            fastafiles[currgene.dbname] =  genomefile(currgene.dbname)        
    seqdict = dict()
    
    for currorg in allorgs:
        
        if faifiles is not None:
            
            currseqs = getseqs(fastafiles[currorg], dbdict[currorg], faindex = faifiles[currorg])
        else:
            currseqs = getseqs(fastafiles[currorg], dbdict[currorg])
        seqdict.update(currseqs)
    return seqdict
    
def getnamedict(genelist):
    namedict = dict()
    for currgene in genelist:
        namedict[currgene.name] = currgene
    return namedict
        
class FastqSeqException(Exception):
    pass
def getseqs(fafile,rangedict, faindex = None):
    if faindex is not None:
        try:
            faifile = fastaindex(fafile, faindex)
        except IOError as e:
            print >>sys.stderr, "Cannot read fasta file "+fafile
            print >>sys.stderr, "Ensure that file "+fafile +" exits and generate fastaindex "+faindex+" with samtools faidx"
            sys.exit(1)
        return faifile.getseqs(rangedict)
    genomefile = open(fafile, "r")
    reheader = re.compile( r"\>([^\s\,]+)")
    allseqs = defaultdict(str)
    currloc = 0
    for line in genomefile:
        line = line.rstrip("\n")
        currheader = reheader.match(line)
        if currheader: #sequence += line[chromstart - currloc:chromend - currloc]
            currseq = currheader.groups(1)[0]
            #print >>sys.stderr, currseq
            currloc = 0
        else:
            for currname, location in rangedict.iteritems():
                if currseq == location.chrom:
                    
                    chromstart = location.start
                    chromend = location.end
                        #print >>sys.stderr, genomefile
                        #print >>sys.stderr, currseq+":"+str(chromstart)+"-"+str(chromend)+":"+str(currloc)
                    if currloc <= chromstart <= currloc + len(line) and currloc <= chromend <= currloc + len(line):
                        allseqs[currname] += line[chromstart - currloc:chromend - currloc]
                    elif currloc <= chromstart <= currloc + len(line):
                        allseqs[currname] += line[chromstart - currloc:]
                    elif currloc <= chromend <= currloc + len(line):
                        allseqs[currname] += line[:chromend - currloc]
                    elif  currloc < chromstart <chromend < currloc + len(line):
                        pass
                    elif chromstart <= currloc < currloc + len(line) < chromend:
                        allseqs[currname] += line
                
            currloc += len(line)
    genomefile.close()
    finalseqs = dict()
    for currname in allseqs.iterkeys():
        #allseqs[currname] = allseqs[currname].upper()
        if (rangedict[currname].strand == "-"):
            seq = list(allseqs[currname].upper())
            seq.reverse()
            comp = {"A":"T","T":"A", "C":"G","G":"C","N":"N","R":"Y","Y":"R","S":"W","W":"S", "K":"M", "M":"K"}
            finalseqs[currname]  = ''.join(comp[base] for base in seq)
        else:
            finalseqs[currname] = allseqs[currname].upper()
    for currseq in rangedict.iterkeys():
        if currseq not in finalseqs:
            #print >>sys.stderr, "**"
            #print >>sys.stderr, fafile
            
            print >>sys.stderr, "No sequence extracted for "+rangedict[currseq].dbname+"."+rangedict[currseq].chrom+":"+str(rangedict[currseq].start)+"-"+str(rangedict[currseq].end)
    return finalseqs        
    
class fastaindex:
    def __init__(self, fafile, faifile):
        self.fafile = fafile
        fai = open(faifile)
        self.chromsize = dict()
        self.chromoffset = dict()
        self.seqlinesize = dict()
        self.seqlinebytes = dict()
        
        for line in fai:
            fields = line.split("\t")
            self.chromsize[fields[0]] = int(fields[1])
            self.chromoffset[fields[0]] = int(fields[2])
            self.seqlinesize[fields[0]] = int(fields[3])
            self.seqlinebytes[fields[0]] = int(fields[4])
    def getchrombed(self, dbname = 'genome'):
        for curr in self.chromsize.iterkeys():
            yield GenomeRange(dbname,curr,0,self.chromsize[curr],name=curr, strand = "+")
    def getseek(self, currchrom,loc):
        #print >>sys.stderr, (self.seqlinebytes[currchrom] - self.seqlinesize[currchrom])
        if currchrom not in self.chromsize:
            raise FastqSeqException("sequence "+currchrom+" not found in index for "+self.fafile)
            #print >>sys.stderr, 
            #sys.exit(1)
        return self.chromoffset[currchrom] + loc + int(loc/(self.seqlinesize[currchrom]))*(self.seqlinebytes[currchrom] - self.seqlinesize[currchrom])
    def getfullseqs(self, names):
        genomefile = open(self.fafile, "r")
        for currchrom in names:
            #print >>sys.stderr, currchrom+":"+str(self.chromsize[currchrom])
            genomefile.seek(self.getseek(currchrom,0))
            #seq = genomefile.read(self.getseek(currchrom,self.chromsize[currchrom]))
            seq = genomefile.read(self.getseek(currchrom,self.chromsize[currchrom]) - self.getseek(currchrom,0))
            #seq = seq.replace("\n","")
            yield currchrom, seq 
    def getseqs(self,  rangedict):
        genomefile = open(self.fafile, "r")
        allseqs = dict()
        for currname, currregion in rangedict.iteritems():
            try:
                currchrom = currregion.chrom
                #faskip = 
                #print >>sys.stderr, int(currregion.start/(self.seqlinebytes[currchrom] - self.seqlinesize[currchrom]))
                genomefile.seek(self.getseek(currchrom,currregion.start))
                seq = genomefile.read(self.getseek(currchrom,currregion.end) - self.getseek(currchrom,currregion.start))
                seq = seq.replace("\n","")
                allseqs[currname] = seq
                #print >>sys.stderr, len(seq)
                #print >>sys.stderr, str(currregion.end - currregion.start)
            except FastqSeqException as e:
                allseqs[currname] = None
                pass
        genomefile.close()
        finalseqs = dict()
        for currname in allseqs.iterkeys():
            if allseqs[currname] is None:
                continue
            #print >>sys.stderr, currname
            #allseqs[currname] = allseqs[currname].upper()
            if (rangedict[currname].strand == "-"):
                
                seq = list(allseqs[currname].upper())
                #print >>sys.stderr, "".join(seq)
                seq.reverse()
                comp = {"A":"T","T":"A", "C":"G","G":"C","N":"N","R":"Y","Y":"R","S":"W","W":"S", "K":"M", "M":"K"}
                finalseqs[currname]  = ''.join(comp[base] for base in seq)
            else:
                finalseqs[currname] = allseqs[currname].upper()
        return finalseqs    
        
        
        
#the object I use for storing groups of genomeranges.  Useful for finding overlaps
class RangeBin:     
    def __init__(self,rangelist, binfactor = 10000):
        self.binfactor = binfactor
        self.bins = []
        self.length = 0
        for curr in rangelist:
            self.additem(curr)
        
    def __len__(self):
        return self.length
    def __iter__(self):
        for currbin in self.bins:
            #print currbin
            for currgene in currbin:
                yield currgene
    def additem(self, item):
        binstart = int(item.start / self.binfactor)
        #print "**"+str(binstart)
        binend = int(item.end / self.binfactor) + 1
        while (binend + 2 >= len(self.bins)):
            self.bins.append(set())
        self.bins[binstart].add(item)
        for i in range(binstart, binend):
            self.bins[i].add(item)
        
        self.length += 1
        #print self.bins[binstart]
    def getrange(self, item):
        #print >>sys.stderr, item.start / self.binfactor

        for i in range(int(item.start / self.binfactor)-1,int(item.end / self.binfactor)+1):
            if i < len(self.bins):
                for currrange in self.bins[i]:
                    if currrange.start >= item.start and currrange.end <= item.end:
                            yield currrange
    #always a pain dealing with long features without losing the optomization                        
    def getbin(self, item):
        #print >>sys.stderr, item.start / self.binfactor
        #print >>sys.stderr, range(int(item.start / self.binfactor)-1,int(item.end / self.binfactor)+1)
        outputregions = set()
        for i in range(int(item.start / self.binfactor)-1,int(item.end / self.binfactor)+1):
            if i < len(self.bins) and i >= 0:
                
                for currrange in self.bins[i]:
                    regionid = tuple([currrange.name, currrange.start, currrange.end])
                    if regionid not in outputregions:
                        outputregions.add(regionid)
                        yield currrange
                    
                        
    def getbinpos(self, item):
        for i in range(int(item / self.binfactor)-1,int(item / self.binfactor)+1):
            if i < len(self.bins) and i >= 0:
                
                for currrange in self.bins[i]:
                    yield currrange
                    
    def getfeatbin(self, name):
        for i, currbin in enumerate(self.bins):
            #print currbin
            for currgene in currbin:
                if currgene.name == name:
                    yield i
    def getbinnums(self, item):
        for i in range(int(item.start / self.binfactor)-1,int(item.end / self.binfactor)+1):
            if i < len(self.bins) and i >= 0:
                yield i
                    
      
def revcom(sequence):
    seq = list(sequence)
    seq.reverse()
    comp = {"A":"T","T":"A", "C":"G","G":"C","N":"N","R":"Y","Y":"R","S":"W","W":"S", "K":"M", "M":"K","-":"-",".":"."}
    return ''.join(comp[base] for base in seq)            
    
    
def get_location(program, allowfail = False):
    progloc = find_executable(program)
    if find_executable(program) is None and not allowfail:
        print >>sys.stderr, "Could not find "+program+" in path"
        print >>sys.stderr, "Aborting"
        sys.exit(1)
    else:
        return progloc
    
def getgithash(scriptdir):
    gitloc = get_location("git")
    
    if gitloc is None:
        print >>sys.stderr, "Cannot find git in path"
        print >>sys.stderr, "Recording of versioning not possible"
    gitjob = subprocess.Popen([gitloc,"--git-dir="+scriptdir+"/.git","rev-parse","HEAD"],stdout = subprocess.PIPE,stderr = subprocess.STDOUT )
    githash = gitjob.communicate()[0].rstrip()
    if gitjob.returncode != 0:
        githash = "Cannot find git hash"
    gitjob = subprocess.Popen([gitloc,"--git-dir="+scriptdir+"/.git","describe"],stdout = subprocess.PIPE,stderr = subprocess.STDOUT )
    gitversion = gitjob.communicate()[0].rstrip()
    if gitjob.returncode != 0:
        gitversion = "Cannot find git version"
    return gitversion, githash 