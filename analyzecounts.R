

library("DESeq2")
library(ggplot2)
library(gridExtra)
library(scales)
library(plyr)
library(ggrepel)

reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
}


colgetlogname =  function(currtable, newname){

newtable = currtable[,c(2),drop = FALSE]
colnames(newtable)[1] <- newname
return(newtable)
}

colgetavgname =  function(currtable, newname){

newtable = currtable[,c(2),drop = FALSE]
colnames(newtable)[1] <- newname
return(newtable)
}

colrename =  function(currtable, newname){

newtable = currtable[,c(5),drop = FALSE]
colnames(newtable)[1] <- newname
return(newtable)
}


args = commandArgs(trailingOnly = TRUE)

experimentname = args[1]
inputtable = args[2]
samplefile = args[3]


readcounts = read.table(inputtable,check.names=FALSE)
sampledata = read.table(samplefile,check.names=FALSE)
#print("*|*")
#print(sampledata)
#print(colnames(readcounts))
#print(gsub("-", ".", sampledata[,1])
#sampleinfo = as.character(sampledata[colnames(readcounts) == sampledata[,1],2])




sampleinfo = as.character(sampledata[colnames(readcounts) == gsub("-", ".", sampledata[,1]) ,2])



#gsub("-", ".", sampledata[,1]) 
#sampledata[,2]
#colnames(readcounts)
#sampledata[,1]

samplenames = unique(sampleinfo)
#combn(unique(sampleinfo),2,simplify = FALSE)
#read.table(args[4], stringsAsFactors = FALSE)
#length(args)
#args[4]

if (length(args) > 3){

pairtable = read.table(args[4], stringsAsFactors = FALSE)
pairreduce = pairtable[pairtable[,1] %in% samplenames & pairtable[,2] %in% samplenames,]

#print(pairtable)
#print(samplenames)

comparisons <- apply(pairreduce, 1, list)
comparisons <- lapply(comparisons,unlist)
}else{
comparisons = combn(unique(sampleinfo),2,simplify = FALSE)
}


coldata = data.frame(condition=factor(sampleinfo))

#print(ncol(readcounts))
#print(nrow(coldata))
#print(head(readcounts))


#print(comparisons)
#print("}}***||*")


cds = DESeqDataSetFromMatrix(countData = readcounts,coldata  ,design = ~ condition)
cds = estimateSizeFactors(cds)
normalizedrnas = sweep(readcounts,2,cds$sizeFactor, "/" )
write.table(normalizedrnas,paste(experimentname,"/",experimentname,"-normalizedreadcounts.txt", sep = ""), sep = "\t")
#print(sizefactors)
write.table(rbind(colnames(readcounts),cds$sizeFactor),file=paste(experimentname,"/",experimentname,"-SizeFactors.txt", sep = ""), row.names=FALSE,col.names=FALSE)
cds = DESeq(cds,betaPrior=TRUE)


#print(coldata)



names = lapply(comparisons, function(currcompare){ })

compareresults = lapply(comparisons, function(currcompare){ list(paste(currcompare[[1]],currcompare[[2]] ,sep= "_"),results( cds, contrast=c("condition", currcompare[[1]] ,currcompare[[2]]),cooksCutoff  =TRUE))})


reslist = lapply(compareresults, function(currresult){colrename(currresult[[2]],currresult[[1]])})

resloglist = lapply(compareresults, function(currresult){colgetlogname(currresult[[2]],currresult[[1]])})

resavglist = lapply(compareresults, function(currresult){colgetlogname(currresult[[2]],currresult[[1]])})                                


#print(compareresults)
#print(length(dispersions(cds)))
#print(nrow(readcounts))

#print(cds)
write.table(cbind(rownames(readcounts),dispersions(cds)),file=paste(experimentname,"/",experimentname,"-dispersions.txt", sep = ""), quote=FALSE,row.names=FALSE,col.names=FALSE)



dds = cds

dashinterc = 1.5


#print adjusted p-values
allprobs = Reduce(function(x,y) cbind(x,y), reslist)
#print("***||*")
#print(reslist)

write.table(allprobs,paste(experimentname,"/",experimentname,"-padjs.txt", sep = ""),sep="	")
#stop("Message")                                                              
#Print log values


alllogvals = Reduce(function(x,y) cbind(x,y), resloglist)
write.table(alllogvals,paste(experimentname,"/",experimentname,"-logvals.txt", sep = ""),sep="	")




outputformat = ".pdf"
if(length(args) > 3){
for (currpair in colnames(alllogvals)){

currlogval = alllogvals[,c(currpair)]
currprob = allprobs[,c(currpair)]
genename = rownames(allprobs)

currsampledata = data.frame(genename, currlogval, currprob)

#print(head(currsampledata))


displaygenes = c("Snora35","Snord116l17", "Mirlet7a-2" ,"Mirlet7b","Mirlet7c-2","Mir138-1", "Mir122", "Mir133a-1")
livergenes =  c("Mir33-201","Mir223","Mir30c-1","Mir144","Mir148a","Mir24-1","Mir29a","Mir122")
musclegenes = c("Mir1a-1","Mir133a-1","Mir208a","Mir208b","Mir499")

displaygenes = c(displaygenes, livergenes, musclegenes)

displaygenes = c()
currsampledata$name = rownames(currsampledata)
displayfeats = ifelse(currsampledata$genename %in% displaygenes, as.character(currsampledata$genename), "")

pvalcutoff = sort(currsampledata$currprob)[10]

displayfeats = ifelse(abs(currsampledata$currlogval) > 1.5 & currsampledata$currprob < pvalcutoff, as.character(currsampledata$genename), "")

#print("**") 
#print(rownames(currsampledata))
#print(head(displayfeats))

pairname = sub( ".", "_",currpair,fixed=TRUE)
pairname = sub( ":", "_",currpair,fixed=TRUE)

#currsampledata = cbind(currlogval,currprob) # 
#print(head(currsampledata))
#currsampledata = currsampledata[currsampledata$currprob > .005,]
#print(head(currsampledata))
currplot <- ggplot(currsampledata, aes_string(x="currlogval", y="currprob")) + geom_point() +scale_x_continuous() +  geom_text_repel(label = displayfeats,min.segment.length = unit(0, 'lines'), segment.color="red")+ scale_y_continuous(trans=reverselog_trans(10))+geom_hline(yintercept = .05, linetype = 2)+geom_hline(yintercept = .005, linetype = 2)+geom_vline(xintercept = dashinterc, linetype = 2) + geom_vline(xintercept = -dashinterc, linetype = 2)+theme_bw() + xlab("Log2-Fold Change")+ylab("Adjusted P-value")+ggtitle(currpair)+theme(legend.box="horizontal",aspect.ratio=1,axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
ggsave(paste(experimentname,"/",pairname ,"-volcano",outputformat,sep= ""), currplot) 

}
}


#print(alllogvals)
#stop("Message")
#Print log values
#print(head(alllogvals))
#head(pairtable)
#head(samplenames)
#head(pairreduce)
colnames(alllogvals) <- paste("log2", colnames(alllogvals), sep = "_")
#print("***||")

colnames(allprobs) <- paste("pval", colnames(allprobs), sep = "_")
allcombinevals = cbind(alllogvals,allprobs)

#write.table(allcombinevals,paste(experimentname,"/",experimentname,"-combine.txt", sep = ""),sep="	", col.names=NA,quote=FALSE) 

sortcombinevals = allcombinevals[order(apply(alllogvals,1,max)),]

#apply(allprobs,1,min) < .05 

#write.table(sortcombinevals,paste(experimentname,"/",experimentname,"-combinesort.txt", sep = ""),sep="	", col.names=NA,quote=FALSE) 

#stop("Message")
#Print out the size factors
#write.table(rbind(colnames(readcounts),dds$sizeFactor),file=paste(experimentname,"/",experimentname,"-SizeFactors.txt", sep = ""), row.names=FALSE,col.names=FALSE)
#stop("Message")
#get deseq normalized  raw counts

#

#allcombined = cbind(allcombinevals,normalizedrnas)
#write.table(allcombined,paste(experimentname,"/",experimentname,"-combineall.txt", sep = ""),sep="	", col.names=NA,quote=FALSE)


#write.table(allcombined[apply(normalizedrnas,1,max) > 30 & apply(alllogvals,1,min) < .05 ,],paste(experimentname,"/",experimentname,"-relevnormalized.txt", sep = ""), col.names=NA )

#write.table(data[abs( log2(data$AlkB_control)- log2(data$AlkB_HCC_1)) > 1.5,],paste(experimentname,"/",experimentname,"-significant.txt", sep = ""))

medcounts = list()

samplenames <- as.character(unique(sampledata[,2]))

#print(samplenames)
for (i in 1:length(samplenames)){
cols <- as.character(sampledata[sampledata[,2] == samplenames[i],1])
#print(samplenames[i])
#print(cols)
#print("**")
if (length(cols) > 1){
#print(samplenames[i])
medcounts[[samplenames[i]]] <- apply(normalizedrnas[,cols], 1, median)

}else{
medcounts[[samplenames[i]]] <- normalizedrnas[,cols]
}
}

#print(medcounts[['Liver_PlusAlkB']]['Snora35'])

#print(medcounts) 
medcountmat <- do.call("cbind",medcounts)

colnames(medcountmat) <- names(medcounts)


write.table(medcountmat,paste(experimentname,"/",experimentname,"-medians.txt", sep = ""),quote=FALSE)


#print(head(medcountmat))
medcountmat = as.matrix(medcountmat)
#print(head(medcountmat["Snora35",]))


allcombinevals = as.matrix(allcombinevals)

#medcounts
#typeof(allcombinevals)
#typeof(medcountmat)
#typeof(medcounts)
#typeof(normalizedrnas)
allcombinevals = cbind(allcombinevals,medcountmat)

#write.table(allcombinevals[apply(readcounts,1,max) > 30 & apply(alllogvals,1,min) < .05 ,],paste(experimentname,"/",experimentname,"-relevnormalizedsamples.txt", sep = ""), col.names=NA )

write.table(allcombinevals,paste(experimentname,"/",experimentname,"-combine.txt", sep = ""), col.names=NA )
