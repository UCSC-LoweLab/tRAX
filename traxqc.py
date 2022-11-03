#!/usr/bin/env python

import sys
import subprocess
import argparse
from tempfile import NamedTemporaryFile
import os
import os.path
import re
from trnasequtils import *
from multiprocessing import Pool
import itertools
from time import localtime, strftime

'''

mergable eteps
variation on raw read read counts
min reads per sample
percentage of reads mapped 60%-70%
read sizes 35-75 for dm-tRNAseq 50-60% ??  15-60 for armseq ???

percentage of reads that map to tRNA  
pca comp ???
number of tRNAs with any reads in at least one sample
number of tRNA reads with reads in at least one sample
whole vs fragment types

check if replicates replicate


'''

redrgb = "rgb(255, 0, 0)"

yellowrgb = "rgb(255,165,0)"

greenrgb = "rgb(60,170,113)"


style = '''
<style>

table, td, th { border: 1px solid black; border-spacing: 0; padding: 1px; text-align: left; }

table { width: 70%; }

</style>
'''

def htmlconvert(message):
    message = message.replace("<","&lt;")
    message = message.replace(">","&gt;")
    return message

class errorset:
    def __init__(self,shortname,alllist,faillist,failcriteria,dimension, faildict, percentformat = False, critfaillist = None, checkfile = None):
        self.shortname = shortname
        self.samples = alllist
        self.failsamples = faillist

        self.criteria = failcriteria
        self.dimension = dimension
        self.failnum = faildict
        self.percentformat = percentformat
        self.critfaillist = critfaillist
        if checkfile is None:
            self.checkfile = ""
        else:
            self.checkfile = checkfile
        
        
        self.warning = sum(faillist) > 0
        self.fail = critfaillist is not None and sum(critfaillist) > 0
        
        self.failset = set(curr for i, curr in enumerate(alllist) if faillist[i]) 
        
    def gettestcolor(self):
        if self.fail:
            return "rgb(0,0,255)"
        elif self.warning:
            return"rgb(255,165,0)"
        else:
            return "rgb(60,170,113)"
    def getteststatus(self):
        if self.fail:
            return "Failed"
        elif self.warning:
            return "Warning" 
        else:
            return "Passed"       
    def getsamplecolor(self, sample):
        if self.critfaillist is not None and sample in self.critfaillist:
            return "rgb(0,0,255)"
        elif sample in self.failset:
            return"rgb(255,165,0)"
        else:
            return "rgb(60,170,113)"
    def getsamplestatus(self, sample):
        #print >>sys.stderr, sample
        #print >>sys.stderr, self.failset 
        if self.critfaillist is not None and sample in self.critfaillist:
            return "Failed"
        elif sample in self.failset:
            return "Warning" 
        else:
            return "Passed"
    def getcriteria(self):
        return htmlconvert(self.failcriteria)
    def getsampleresult(self,currsample):
        #print >>sys.stderr, self.failnum
        if currsample not in self.failnum:
            return "??"
        if self.percentformat:
            
            return "{0:.2f}%".format(100 * self.failnum[currsample])
        else:
            return "{0:.2f}".format(self.failnum[currsample])



class seqprepinfo:
    def __init__(self, merged, unmerged, discarded):
        self.merged = merged
        self.unmerged = unmerged
        self.discarded = discarded
    def gettotal(self, sample):
        return self.merged[sample] + self.unmerged[sample] + self.discarded[sample]
    def getmergedpercent(self, sample):
        return self.merged[sample] / (1.*self.gettotal(sample) + .01)
    def getsamples(self):
        return tuple(self.merged.keys())
        
class cutadaptinfo:
    def __init__(self, trimmed, untrimmed, discarded):
        self.trimmed = trimmed
        self.untrimmed = untrimmed
        self.discarded = discarded
    def gettotal(self, sample):
        return self.trimmed[sample] + self.untrimmed[sample] + self.discarded[sample]
    def getpassedpercent(self, sample):
        return (self.trimmed[sample] + self.untrimmed[sample]) / (1.*self.gettotal(sample))
    def getsamples(self):
        return tuple(self.trimmed.keys())
        
        
        
def getreadprep(prepfilename, manifestfilename, sampleinfo):
    manifestfile = open(manifestfilename)
    samplenames = dict()
    for currline in manifestfile:
        fields = currline.rstrip().split("\t")
        samplenames[fields[0]] = sampleinfo.getfastqsample(fields[1])
    prepfile = open(prepfilename)
    samples = dict()
    unmerged = dict()
    discarded = dict()
    merged = dict()
    for i, currline in enumerate(prepfile):
        fields = currline.rstrip().split("\t")
        if i == 0:
            runsamples = list(samplenames[curr] for curr in fields)
            continue
            
        if len(fields) != len(runsamples) + 1:
            continue
        for j in range(0, len(runsamples)):
            if fields[0] == "merged":
                merged[runsamples[j]] = int(fields[j + 1])
            elif fields[0] == "unmerged":
                unmerged[runsamples[j]] = int(fields[j + 1])
            elif fields[0] == "discarded":
                discarded[runsamples[j]] = int(fields[j + 1])
                
    return seqprepinfo(merged, unmerged, discarded)

def getcutadapt(prepfilename, manifestfilename, sampleinfo):
    manifestfile = open(manifestfilename)
    samplenames = dict()
    for currline in manifestfile:
        fields = currline.rstrip().split("\t")
        samplenames[fields[0]] = sampleinfo.getfastqsample(fields[1])
    prepfile = open(prepfilename)
    trimmed = dict()
    untrimmed = dict()
    discarded = dict()
    for i, currline in enumerate(prepfile):
        fields = currline.rstrip().split("\t")
        if i == 0:
            runsamples = list(samplenames[curr] for curr in fields)
            continue
            
        if len(fields) != len(runsamples) + 1:
            continue
        for j in range(0, len(runsamples)):
            if fields[0] == "trimmed":
                trimmed[runsamples[j]] = int(fields[j + 1])
            elif fields[0] == "untrimmed":
                untrimmed[runsamples[j]] = int(fields[j + 1])
            elif fields[0] == "discarded":
                discarded[runsamples[j]] = int(fields[j + 1])
                
    return cutadaptinfo(trimmed, untrimmed, discarded)

minmergepercent = .6
def checkreadprep(allpreps, sampleinfo):
    prepdict = dict()
    for prepinfo in allpreps:   
        if os.path.exists(prepinfo+"_sp.txt"):
            prepresults = getreadprep(prepinfo+"_sp.txt", prepinfo+"_manifest.txt",sampleinfo)
            samples = prepresults.getsamples()
            prepdict.update({currsample : prepresults.getmergedpercent(currsample) for currsample in samples}) 
        if os.path.exists(prepinfo+"_ca.txt"):
            prepresults = getcutadapt(prepinfo+"_ca.txt", prepinfo+"_manifest.txt",sampleinfo)
            samples = prepresults.getsamples()
            prepdict.update({currsample : prepresults.getpassedpercent(currsample) for currsample in samples})
    #print >>sys.stderr, prepdict
    #print >>sys.stderr, len(list(currsample for currsample in samples if currsample in prepdict))
    #print >>sys.stderr,samples
    #print >>sys.stderr, "**||"
    if len(set(sampleinfo.getsamples()) & set(prepdict.keys())) == 0:
        return list()
    

    
    lowmergesamples = list(prepdict[currsample] < minmergepercent for currsample in samples)
    
    mergeerr = errorset("merging_rate",samples, lowmergesamples, "Sequencing read merging rate  > "+str(100*minmergepercent)+"%"+"","Merging Rate", prepdict, percentformat = True, checkfile = prepinfo+"_sp.pdf")
    
    return [mergeerr]
    
def percentform(innum):
    return  "{0:.2f}%".format(100 * innum)
def countform(innum):
    return "{0:.2f}".format(innum)
#print  str(len(missingtrnasamples)) +" samples have low tRNA read counts ( > "+str(100*minactivepercent)+"% between "+str(minsizethreshold)+" and "+str(maxsizethreshold)+") [" +",".join(currsample+":"+str(thresholdreadpercent[currsample]) for currsample in missingtrnasamples)+"]"


def failcolor(failed, colortrue, colorfalse):
    if failed:
        return colorfalse
    else:
        return colortrue

def errorline(alllist,faillist, failmessage, failcriteria, faildict, percentformat = False, critfaillist = None, checkfile = None, outputfile = sys.stdin):
    if percentformat:
        outformat = percentform
    else:
        outformat = countform
    color = greenrgb

    if sum(faillist) > 0:
        color = yellowrgb
    endstring = ""    
    if checkfile is not None:
        endstring = checkfile
    failstring = ""

    failstring = " [" +",".join('<b style="color:'+failcolor(faillist[i],greenrgb,yellowrgb)+';">'+currsample+"</b>"+":"+outformat(faildict[currsample]) for i, currsample in enumerate(alllist))+"]"
        
    print("<p>", file=outputfile)
    print('<text style="color:'+color+';">'+str(sum(faillist)) +" samples</text> " + failmessage +" ( "+failcriteria+") " +endstring+"<br/>"+ failstring, file=outputfile)
    
    if critfaillist is not None and sum(critfaillist) > 0:
        color = redrgb
        print('<text style="color:'+color+';">'+str(len(critfails)) +" samples</text> " + failmessage +" ( "+failcriteria+") [" +",".join(currsample+":"+outformat(faildict[currsample]) for currsample in faillist)+"]", file=outputfile)
    print("</p>", file=outputfile)
    
def errorsingle(fail, failmessage, failcriteria,critfail = False):

    color = greenrgb
    message = "No"
    if critfail:
        color = redrgb
        message = "Present"
    elif fail:
        color = yellowrgb
        message = "Critical"
    

    print("<p>")
    print('<text style="color:'+color+';">'+message +"</text> " +failmessage +" ( "+failcriteria+") ")
    print("</p>")
    
class mappingresults:
    def __init__(self, unmap, single, multi):
        self.unmap = unmap
        self.single = single
        self.multi = multi
        
    def totalreadscount(self, sample):
        return self.unmap[sample] + self.single[sample] + self.multi[sample]
    def getmappercent(self, sample):
        totalreads = self.totalreadscount(sample)
        return (totalreads  - self.unmap[sample]) / (1.*totalreads+.01) 
        
minmapreads = 200000
              
minmappercent = .65
def checkreadsmapping(samplename, sampleinfo, tgirtmode = False):
    mapresults = getreadmapping(samplename, sampleinfo)
    samples = sampleinfo.getsamples()
    totalreads = {currsample : mapresults.totalreadscount(currsample) for currsample in samples}
    mappercent = {currsample : mapresults.getmappercent(currsample) for currsample in samples}
    
    lowcountsamples = list(totalreads[currsample] < minmapreads for currsample in samples)
    #print >>sys.stderr, lowcountsamples
    #print str(len(lowcountsamples)) +" contain fewer mappable reads than recommended minimum ("+str(minmapreads)+") [" +",".join(currsample+":"+str(totalreads[currsample]) for currsample in lowcountsamples)+"]"
    lowcounterr = errorset("mappable_read",samples, lowcountsamples, "Mappable reads > "+str(minmapreads)+"","Read Count", totalreads, checkfile = samplename+"-mapinfo.pdf")

    lowmapsamples = list(mappercent[currsample] < minmappercent  for currsample in samples)
    lowmaperr = errorset("mappable_rate",samples,lowmapsamples, "Mapping rate > "+str(100*minmappercent)+"%","Mapping Rate", mappercent, percentformat = True, checkfile = samplename+"-mapinfo.pdf")
    #print str(len(lowmapsamples)) +" have lower mappable read percentage than recommended minimum ("+str(100*minmappercent)+"%) [" +",".join(currsample+":"+str(mappercent[currsample]) for currsample in lowmapsamples)+"]"
    return [lowcounterr, lowmaperr]
def getmapfile(samplename):
    return samplename+"/"+samplename+"-mapinfo.txt"
def getreadmapping(samplename, sampleinfo):
    mappingcounts = dict()
    mapresults = open(getmapfile(samplename))
    allsamples = sampleinfo.getsamples()
    runsamples = None
    multimaps = dict()
    singlemaps = dict()
    unmaps = dict()
    for i, currline in enumerate(mapresults):
        fields = currline.rstrip().split("\t")
        if i == 0:
            runsamples = list(fields)
            if set(runsamples) != set(allsamples):
                print(runsamples, file=sys.stderr)
                print(allsamples, file=sys.stderr)
                print("QAError", file=sys.stderr)
            continue
            
        if len(fields) != len(allsamples) + 1:
            continue
        for j in range(0, len(runsamples)):
            if fields[0] == "unmap":
                unmaps[runsamples[j]] = int(fields[j + 1])
            elif fields[0] == "single":
                singlemaps[runsamples[j]] = int(fields[j + 1])
            elif fields[0] == "multi":
                multimaps[runsamples[j]] = int(fields[j + 1])
                
    return mappingresults(unmaps, singlemaps, multimaps)

def gettypefile(samplename):
    return samplename+"/"+samplename+"-typerealcounts.txt"
def getreadlengthfile(samplename):
    return samplename+"/"+samplename+"-readlengths.txt"
def getfragtypefile(samplename):
    return samplename+"/"+samplename+"-fragtypes.txt"
def gettrnacountfile(samplename):
    return samplename+"/"+samplename+"-trnacounts.txt"
def getsizefactorfile(samplename):
    return samplename+"/"+samplename+"-SizeFactors.txt"
    
class typecount:
    def __init__(self, typecounts):
        self.typecounts = typecounts
    def gettotal(self, sample):
        return sum(self.typecounts[sample].values()) 
    def gettrnapercent(self, sample):
        return (self.typecounts[sample]["tRNA"] + self.typecounts[sample]["pretRNA"] )/ (1.*self.gettotal(sample)+.01)
    def getrrnapercent(self, sample):
        if "rRNA" in self.typecounts[sample]:
            
            return self.typecounts[sample]["rRNA"] / (1.*self.gettotal(sample)+.01)
        else:
            return None
    def getotherpercent(self, sample):
        return (self.typecounts[sample]["other"]) / (1.*self.gettotal(sample))
def gettypecounts(samplename, sampleinfo):
    typeresults = open(gettypefile(samplename))
    allsamples = sampleinfo.getsamples()
    runsamples = None
    typecounts = defaultdict(dict)
    #print >>sys.stderr, gettypefile(samplename)
    for i, currline in enumerate(typeresults):
        fields = currline.rstrip().split("\t")
        if i == 0:
            runsamples = list(fields)
            
            if set(runsamples) != set(allsamples):
                print(runsamples, file=sys.stderr)
                print(allsamples, file=sys.stderr)
                print("QAError", file=sys.stderr)
            continue
            
           
        if len(fields) != len(allsamples) + 1:
            print(runsamples, file=sys.stderr)
            print(fields, file=sys.stderr)
            print("QAError", file=sys.stderr)
            continue
        for j in range(0, len(runsamples)):
            typecounts[runsamples[j]][fields[0]] = int(fields[j + 1])
    return typecount(typecounts)
    
'''
Length	Sample	other	trnas	pretrnas
0	nuc_rep2	0	0	0
1	nuc_rep2	0	0	0
'''
def getmeanfreq(freqtable):
    return sum(curr * freqtable[curr] for curr in list(freqtable.keys())) / (1.*sum(freqtable.values()) +.01)
class lengthcount:
    def __init__(self, trnalengthcounts, pretrnalengthcounts, otherlengthcounts):
        self.trnalengthcounts = trnalengthcounts
        self.pretrnalengthcounts = pretrnalengthcounts
        self.otherlengthcounts = otherlengthcounts
        self.samples = set(itertools.chain(list(trnalengthcounts.keys()), list(pretrnalengthcounts.keys()), list(otherlengthcounts.keys())))
        
        #print >>sys.stderr, list(itertools.chain(list(trnalengthcounts[currsample].keys() for currsample in self.samples),list(pretrnalengthcounts[currsample].keys() for currsample in self.samples),list(otherlengthcounts[currsample].keys() for currsample in self.samples)))

        
        self.maxlength = max(itertools.chain(itertools.chain.from_iterable(list(trnalengthcounts[currsample].keys()) for currsample in self.samples),itertools.chain.from_iterable(list(pretrnalengthcounts[currsample].keys()) for currsample in self.samples),itertools.chain.from_iterable(list(otherlengthcounts[currsample].keys()) for currsample in self.samples)))
        
    def getalllengths(self, sample):
        return {currlength: self.trnalengthcounts[sample][currlength] + self.trnalengthcounts[sample][currlength] + self.trnalengthcounts[sample][currlength] for currlength in range(self.maxlength)}
    def getsamplemean(self, sample):
        return getmeanfreq(self.getalllengths(sample))
    def getthreshold(self, sample, minsize, maxsize):
        alllengths = self.getalllengths(sample)
        maxsize = min([maxsize, max(alllengths.keys())])
        return sum(alllengths[i] for i in range(minsize, maxsize))
    def getthresholdpercent(self, sample, minsize, maxsize):
        
        return self.getthreshold(sample, minsize, maxsize) / (1.*sum(self.getalllengths(sample).values()))

def getreadlengths(samplename,sampleinfo):
    lengthresults = open(getreadlengthfile(samplename))
    allsamples = sampleinfo.getsamples()
    runsamples = None
    trnalengthcounts = defaultdict(lambda: defaultdict(int))
    otherlengthcounts = defaultdict(lambda: defaultdict(int))
    pretrnalengthcounts = defaultdict(lambda: defaultdict(int))
    #print >>sys.stderr, gettypefile(samplename)
    for i, currline in enumerate(lengthresults):
        fields = currline.rstrip().split("\t")
        if i == 0:
            continue
           
        if len(fields) < 4:
            continue
            pass
        trnalengthcounts[fields[1]][int(fields[0])] = int(fields[3])
        otherlengthcounts[fields[1]][int(fields[0])] = int(fields[2])
        pretrnalengthcounts[fields[1]][int(fields[0])] = int(fields[4])
    return lengthcount(trnalengthcounts, pretrnalengthcounts, otherlengthcounts)
    



def getfragtypes(samplefile, trnainfo):
    pass


def filelink(filename): 
    return '<a href="'+filename+'">'+filename+'</a>' #file:///

def checkreadtypes(samplename, sampleinfo, tgirtmode = False):
    trnapercentcutoff = .05
    ribopercentcutoff = .35
    unmappercentcutoff = .35
    
    highmeanlength = 40
    minsizethreshold = 15
    maxsizethreshold = 50
    percentsizethreshold = .70
    lowmeanlength = None
    if tgirtmode:
        
        trnapercentcutoff = .5
        ribopercentcutoff = .35
        unmappercentcutoff = .35
        
        highmeanlength = 75
        lowmeanlength = 40
        minsizethreshold = 40
        maxsizethreshold = 75
        percentsizethreshold = .70
    typecounts = gettypecounts(samplename, sampleinfo)
    getfragtypes(samplename, sampleinfo)
    
    samples = sampleinfo.getsamples()
    
    trnapercent = {currsample : typecounts.gettrnapercent(currsample) for currsample in samples}
    lowtrnasamples = list(trnapercent[currsample] < trnapercentcutoff for currsample in samples)
    #print  str(len(lowtrnasamples)) +" samples have low tRNA read percentage ( < "+str(100*trnapercentcutoff)+"%) [" +",".join(currsample+":"+str(trnapercent[currsample]) for currsample in lowtrnasamples)+"]"
    lowtrnaerr = errorset("trna_read",samples, lowtrnasamples, "tRNA read share > "+str(100*trnapercentcutoff)+"%", "tRNA Read Percentage",trnapercent, percentformat = True,checkfile = samplename+"-typecounts.pdf")

    rrnapercent = {currsample : typecounts.getrrnapercent(currsample) for currsample in samples}
    highribosamples = list(rrnapercent[currsample] > ribopercentcutoff  for currsample in samples)
    #print str(len(highribosamples)) +" samples have high rRNA read percentage ( > "+str(100*ribopercentcutoff)+"%) [" +",".join(currsample+":"+str(rrnapercent[currsample]) for currsample in highribosamples)+"]"
    highriboerr = errorset("rrna_read",samples, highribosamples, "rRNA read share < "+str(100*ribopercentcutoff)+"%","rRNA Read Percentage", rrnapercent, percentformat = True, checkfile = samplename+"-typecounts.pdf")

    otherpercent = {currsample : typecounts.getotherpercent(currsample) for currsample in samples}
    highothersamples = list(otherpercent[currsample] > unmappercentcutoff for currsample in samples)
    #print  str(len(highothersamples)) +" samples have many reads not mapping to annotated genes ( > "+str(100*unmappercentcutoff)+"%) [" +",".join(currsample+":"+str(otherpercent[currsample]) for currsample in highothersamples)+"]"
    highothererr = errorset("unannotated",samples, highothersamples, "Reads mapping to unannotated regions < "+str(100*unmappercentcutoff)+"%","Unannotated Region Mapping Rate", otherpercent, percentformat = True, checkfile = samplename+"-typecounts.pdf")

    allreadlength = getreadlengths(samplename, sampleinfo)

    meanreadlength = {currsample : allreadlength.getsamplemean(currsample) for currsample in samples}
    meanreadsamples = list(allreadlength.getsamplemean(currsample) > highmeanlength for currsample in samples)
    #print  str(len(meanreadsamples)) +" samples have high read length average ( > "+str(highmeanlength)+") [" +",".join(currsample+":"+str(rrnapercent[currsample]) for currsample in meanreadsamples)+"]"
    highlengtherr = errorset("high_read_len",samples, meanreadsamples, "Read length average < "+str(highmeanlength)+ " bases", "Average Read Length",meanreadlength, percentformat = False, checkfile = samplename+"-readlengths.pdf")
    
    lowmeanlenerr = None
    if lowmeanlength is not None:
        lowmeanreadsamples = list(allreadlength.getsamplemean(currsample) < lowmeanlength for currsample in samples)
        lowmeanlenerr = errorset("low_read_len",samples, lowmeanreadsamples, "Read length average > "+str(lowmeanlength)+ " bases", "Average Read Length",meanreadlength, percentformat = False, checkfile = samplename+"-readlengths.pdf")
    
    
    
    thresholdreadpercent = {currsample : allreadlength.getthresholdpercent(currsample, minsizethreshold, maxsizethreshold) for currsample in samples}
    badsizesamples = list(thresholdreadpercent[currsample] <  percentsizethreshold for currsample in samples )
    
    #print  str(len(badsizesamples)) +" samples have high rRNA read percentage ( > "+str(100*percentsizethreshold)+"% between "+str(minsizethreshold)+" and "+str(maxsizethreshold)+") [" +",".join(currsample+":"+str(rrnapercent[currsample]) for currsample in badsizesamples)+"]"
    badsizeerr = errorset("trna_sizes",samples, badsizesamples, " >= "+str(100*percentsizethreshold)+"% of reads between "+str(minsizethreshold)+" and "+str(maxsizethreshold) + " bases", "Read Percentage",thresholdreadpercent, percentformat = True, checkfile = samplename+"-readlengths.pdf")


    return [lowtrnaerr,highriboerr,highothererr,highlengtherr, badsizeerr]

class trnacount:
    def __init__(self, trnacounts):
        self.trnacounts = trnacounts
    def gettrnaactive(self, currsample, cutoff = 20):
        return sum(1 for curr in self.trnacounts[currsample].keys() if int(self.trnacounts[currsample][curr]) > cutoff)
    def gettrnaactivepercent(self, currsample, trnainfo, cutoff = 20):
        #print >>sys.stderr, self.gettrnaactive(currsample, cutoff)
        #print >>sys.stderr, (1.*len(trnainfo.gettranscripts()))
        #print >>sys.stderr, self.trnacounts[currsample].keys()
        return self.gettrnaactive(currsample, cutoff)/ (1.*len(trnainfo.gettranscripts())+.01)
        
        
def gettrnacounts(samplename, sampleinfo, trnainfo):
    typeresults = open(gettrnacountfile(samplename))
    allsamples = sampleinfo.getsamples()
    trnatranscripts = set(trnainfo.gettranscripts())
    runsamples = None
    trnacounts = defaultdict(dict)
    #print >>sys.stderr, gettypefile(samplename)
    for i, currline in enumerate(typeresults):
        fields = currline.rstrip().split("\t")
        if i == 0:
            runsamples = list(fields)
            if set(runsamples) != set(allsamples):
                print(runsamples, file=sys.stderr)
                print(allsamples, file=sys.stderr)
                print("QAError", file=sys.stderr)
            continue
           
        if len(fields) != len(allsamples) + 1:
            print(runsamples, file=sys.stderr)
            print(fields, file=sys.stderr)
            print("QAError", file=sys.stderr)
            continue
        if fields[0] in trnatranscripts:    
            for j in range(0, len(runsamples)):
                
                
                trnacounts[runsamples[j]][fields[0]] = int(fields[j + 1])
    return trnacount(trnacounts)
class sizefactor:
    def __init__(self, sizefactors):
        self.sizefactors = sizefactors

def getsizefactor(samplename, sampleinfo):
    typeresults = open(getsizefactorfile(samplename))
    allsamples = sampleinfo.getsamples()
    runsamples = None
    sizefactors = dict()
    for i, currline in enumerate(typeresults):
        fields = currline.rstrip().split()
        fields = list(curr.strip('"') for curr in fields)
        if i == 0:
            runsamples = list(fields)
            if set(runsamples) != set(allsamples):
                print(list(runsamples), file=sys.stderr)
                print(list(allsamples), file=sys.stderr)
                print("QAError", file=sys.stderr)
            continue
           
        if len(fields) != len(allsamples):
            print(len(runsamples), file=sys.stderr)
            print(len(fields), file=sys.stderr)
            print("QAError", file=sys.stderr)
            continue
        #print >>sys.stderr, "**"
        for j in range(0, len(runsamples)):
            sizefactors[runsamples[j]] = float(fields[j])
    return sizefactor(sizefactors)
    
sizefactordiff = 3.

minactivepercent = .5

minreadcount = 20
def checkgenecounts(samplename, sampleinfo, trnainfo, tgirtmode = False):
    #print >>sys.stderr, "**"+samplename
    readcounts = gettrnacounts(samplename, sampleinfo, trnainfo)
    sizefactors = getsizefactor(samplename, sampleinfo)
    samples = sampleinfo.getsamples()
    
    thresholdreadpercent = {currsample : readcounts.gettrnaactivepercent(currsample, trnainfo, cutoff = minreadcount) for currsample in samples}
    
    
    missingtrnasamples = list(thresholdreadpercent[currsample] < minactivepercent for currsample in samples)
    #print  str(len(missingtrnasamples)) +" samples have low tRNA read counts ( > "+str(100*minactivepercent)+"% between "+str(minsizethreshold)+" and "+str(maxsizethreshold)+") [" +",".join(currsample+":"+str(thresholdreadpercent[currsample]) for currsample in missingtrnasamples)+"]"
    lowcounterr = errorset("trna_read_count",samples, missingtrnasamples, ">= "+str(100*minactivepercent)+"% of tRNAs with more than "+str(minreadcount) + " reads","Percentage of tRNAs" ,thresholdreadpercent, percentformat = True, checkfile = samplename+"-tRNAcounts.txt")
    
    #gotta fix ***
    samplesizefactors = {currsample : sizefactors.sizefactors[currsample] for currsample in samples}

    badsizefactors = list(samplesizefactors[currsample] > sizefactordiff or samplesizefactors[currsample] < (1./(sizefactordiff+.0001)) for currsample in samples)
    sizefactorerr = errorset("size_factors",samples, badsizefactors, "DESeq2 size factor differences < "+str(sizefactordiff)+"x","Size Factor" ,samplesizefactors,checkfile = samplename+"-SizeFactors.txt")

    #errorsingle(min(sizefactors.sizefactors.values())*sizefactordiff <  min(sizefactors.sizefactors.values()), "large DESeq2 sizefactor differences", " > "+str(sizefactordiff)+"x")
    
    #(failmessage, failcriteria, critfail = False)
    #print  "Large DESeq2 sizefactor differences ( >"+sizefactordiff+" ) [" +",".join(currsample+":"+str(sizefactors.sizefactors[currsample]) for currsample in badsizesamples)+"]"

    return [lowcounterr, sizefactorerr]

def checktrnamappings(sampleinfo):
    pass



def checkfragmenttypes(samplefile, trnainfo):
    pass


def readtrimindex(trimindex):
    filelocs = dict()
    indexfile = open(trimindex)
    for currline in indexfile:
        fields = currline.split()
        if len(fields) > 1:
            filelocs[fields[0]] = fields[1]
    indexfile.close()
    return list(filelocs.keys())


mode = '''<h2>TRAX Data Quality Report</h2>

<h4>

Date: {date}<br/>

Run mode: {mode}

</h4>
'''

'''
<a name="merging_rate"><h4>Sequencing read merging rate >= 60% (nuclearsamples_sp.pdf)</h4></a>

<table>

<thead><tr><th width="15%">Status</th><th width="50%">Sample</tH><th>Merging Rate</th></tr></thead>

<tbody>

<tr><td><b style="color:rgb(255,165,0);">Warning</b></td><td>SRR5757128</td><td>55.10%</td></tr>

<tr><td><b style="color:rgb(255,165,0);">Warning</b></td><td>SRR5757129</td><td>30.06%</td></tr>

<tr><td><b style="color:rgb(255,165,0);">Warning</b></td><td>SRR5757133</td><td>33.90%</td></tr>

<tr><td><b style="color:rgb(255,165,0);">Warning</b></td><td>SRR5757132</td><td>20.16%</td></tr>

<tr><td><b style="color:rgb(255,0,0);">Fail</b></td><td>SRR5757137</td><td>2.99%</td></tr>

<tr><td><b style="color:rgb(255,0,0);">Fail</b></td><td>SRR5757136</td><td>4.15%</td></tr>

</tbody>

</table>

<p><a href="#summary">back to Summary</a></p>

<hr>
'''


def main(**args):
    samplename = args["experimentname"]
    runname = None
    if "runname" in args:
        runname = args["runname"]
    if args["output"] is not None:
        outputfile = open(args["output"], "w")
    else:
        outputfile = sys.stdout
    sampleinfo = samplefile(os.path.expanduser(args["samplefile"]))
    trnainfo = transcriptfile(os.path.expanduser(args["databasename"] + "-trnatable.txt"))
    
    tgirtmode = args["tgirt"]
    
    allsamples = sampleinfo.getsamples()
    print("<html>", file=outputfile)
    print("<head>"+style+"</head>", file=outputfile)
    if tgirtmode: 
        modestring = "Full-length tRNAs"
    else:
        modestring = "tRNA fragments"
    date = strftime("%A %B %d, %Y", localtime())
    print(mode.format(date= date,mode =modestring), file=outputfile)
    

    print("<body>", file=outputfile)
    #if tgirtmode:
    #    print >>outputfile, "<p>In TGIRT mode for tRNA transcript analysis</p>"
    #else:
    #    print >>outputfile, "<p>In ARMSeq mode for tRNA fragment analysis</p>"
        
    print('<a name="summary"><h4>Summary</h4></a>', file=outputfile)
    prepresults = list()
    if runname is None and os.path.exists("trimindex.txt"):
        runnames = readtrimindex("trimindex.txt")
        prepresults = checkreadprep(runnames,sampleinfo)
    if runname is not None:
        prepresults = checkreadprep([runname],sampleinfo)
    mappingresults = list()
    if os.path.exists(getmapfile(samplename)):
        mappingresults = checkreadsmapping(samplename, sampleinfo, tgirtmode)
    typeresults = checkreadtypes(samplename, sampleinfo, tgirtmode)
    countresults = checkgenecounts(samplename, sampleinfo, trnainfo, tgirtmode)
    
    
    if prepresults:
        allresults = prepresults+mappingresults+typeresults+countresults
    else:
        allresults = prepresults+mappingresults+typeresults+countresults
    
    print("<p>", file=outputfile)
    for currtest in allresults:
        color = "rgb(60,170,113)"
        errlvl = currtest.getteststatus()
        color = currtest.gettestcolor()
        print('<b style="color:{color};">{msg}</b> <a href="#{testname}">{criteria}</a> ({filename})</br>'.format(color = color, msg = errlvl, testname = currtest.shortname, criteria = currtest.criteria, filename = filelink(currtest.checkfile)), file=outputfile)

    print("</p>", file=outputfile)
    
    print("<hr>\n<hr>", file=outputfile)
    
    for currtest in allresults:
        print('<a name="{testname}"><h4>{msg}</h4></a>'.format(testname = currtest.shortname,msg = currtest.criteria), file=outputfile)

        print('<table>', file=outputfile)

        print('<thead><tr><th width="15%">Status</th><th width="50%">Sample</tH><th>{measure}</th></tr></thead>'.format(measure = currtest.dimension), file=outputfile)

        print('<tbody>', file=outputfile)

        for currsample in allsamples:
            color = currtest.getsamplecolor(currsample)
            errlvl = currtest.getsamplestatus(currsample)
            print('<tr><td><b style="color:{color};">{errlvl}</b></td><td>{samplename}</td><td>{sampleresult}</td></tr>'.format(color = color, errlvl = errlvl, samplename=currsample,sampleresult=currtest.getsampleresult(currsample)), file=outputfile)



        print('</tbody>\n\n</table>\n\n<p><a href="#summary">back to Summary</a></p>\n\n<hr>', file=outputfile)

    print("</html>", file=outputfile)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Map reads with bowtie2 and process mappings')
    parser.add_argument('--runname',
                       help='run name of trimadapters.py')
    parser.add_argument('--samplefile',required=True,
                       help='Sample file in format')
    parser.add_argument('--databasename',required=True,
                       help='tRNA file in format')
    parser.add_argument('--experimentname',required=True,
                       help='Sample file in format')
    parser.add_argument('--tgirt', action="store_true", default=False,
                   help='tgirt mode')
    
    parser.add_argument('--output',
                       help='output file if not stdout')
    
    


    
    args = parser.parse_args()
    main(**vars(args))

