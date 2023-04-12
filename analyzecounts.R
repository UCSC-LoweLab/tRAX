

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



cds = DESeqDataSetFromMatrix(countData = readcounts,coldata  ,design = ~ condition)
cds = estimateSizeFactors(cds)
normalizedrnas = sweep(readcounts,2,cds$sizeFactor, "/" )
write.table(normalizedrnas,paste(experimentname,"/",experimentname,"-normalizedreadcounts.txt", sep = ""), sep = "\t")
#print(sizefactors)
write.table(rbind(colnames(readcounts),cds$sizeFactor),file=paste(experimentname,"/",experimentname,"-SizeFactors.txt", sep = ""), row.names=FALSE,col.names=FALSE)



if(length(unique(coldata$condition)) == ncol(readcounts)) {
q()
}
cds = DESeq(cds,betaPrior=TRUE)

deseq2Data <- cds

#not sure why this fails sometimes

heatmaps = FALSE
if(heatmaps){



#Can fail sometimes
deseq2VST <- vst(deseq2Data,fitType = "local")
write.table(as.data.frame(assay(deseq2VST)),paste(experimentname,"/",experimentname,"-vst.txt", sep = ""), col.names=NA )

deseq2rlog <- rlog(deseq2Data)  
write.table(as.data.frame(assay(deseq2rlog)),paste(experimentname,"/",experimentname,"-rlog.txt", sep = ""), col.names=NA )


zscores = as.data.frame(t(scale(t(normalizedrnas))))
write.table(zscores,paste(experimentname,"/",experimentname,"-zscore.txt", sep = ""), col.names=NA )



}
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

#currpair in colnames(alllogvals)
for (currcomp in comparisons){

currpair = paste(currcomp[1], currcomp[2], sep ="_") 
#heatmapstuff


currlogval = alllogvals[,c(currpair)]
currprob = allprobs[,c(currpair)]
genename = rownames(allprobs)

#pairname = sub( ".", "_",currpair,fixed=TRUE)
pairname = sub( ":", "_",currpair,fixed=TRUE)
currsampledata = data.frame(genename, currlogval, currprob)

#print(head(allprobs))
#print("$$**")
#print(head(alllogvals))

if(heatmaps){
library(viridis)

#deseq2Results <- results(cds, contrasts = list(currpair[1], currpair[2]))
#print(head(deseq2Results))
#deseq2ResDF <- as.data.frame(deseq2Results) 
#print("**")
deseq2trans = deseq2rlog
# Convert the DESeq transformed object to a data frame
deseq2trans <- assay(deseq2trans)
deseq2trans <- as.data.frame(deseq2trans)

#center on mean
deseq2trans = sweep(deseq2trans, MARGIN=1, STATS= rowMeans(deseq2trans))

deseq2trans$Gene <- rownames(deseq2trans)

sigGenes <- currsampledata$genename[abs(currsampledata$currlogval) > 1.5 & currsampledata$currprob < .05]


deseq2trans <- deseq2trans[deseq2trans$Gene %in% sigGenes,]

library(reshape2)

# First compare wide vs long version
deseq2trans_wide <- deseq2trans
deseq2trans_long <- melt(deseq2trans, id.vars=c("Gene"))

#head(deseq2trans)
#head(deseq2trans)

# Now overwrite our original data frame with the long format
deseq2trans <- melt(deseq2trans, id.vars=c("Gene"))
#print("**")
#print(currcomp[[1]])
#print(currcomp[[2]])
#print(coldata)
samplenames = as.character(sampledata[,1])
#print(samplenames)
#print(samplenames[coldata$condition == currcomp[[1]] | coldata$condition == currcomp[[2]]])
currsamples = samplenames[coldata$condition == currcomp[[1]] | coldata$condition == currcomp[[2]]]
#print(currsamples)
#print(head(deseq2trans))
deseq2trans = deseq2trans[deseq2trans$variable %in% currsamples,]
print("()")

maxscore = max(abs(deseq2trans$value))
# Make a heatmap
#print(head(deseq2trans))
heatmap <- ggplot(deseq2trans, aes(x=variable, y=Gene, fill=value)) + geom_raster() +  theme_bw() +  scale_fill_gradient2(high="darkred",low="darkblue", limits = c(-maxscore, maxscore))  + ggtitle(paste(currcomp[[1]], currcomp[[2]], sep = " vs "))+theme(axis.text.x=element_text(angle=65, hjust=1),  axis.ticks.y=element_blank()) + labs(fill = "log-fold change from median") #axis.text.y=element_blank(),  scale_fill_viridis(discrete=FALSE) scale_fill_distiller(palette = "RdBu")
#print(paste(currcomp[[1]], currcomp[[2]], sep = " vs "))
#print(1*length(unique(deseq2trans$Gene)))
ggsave(paste(experimentname,"/",pairname,"-deseqlogheatmap",".pdf",sep= ""),height=4+ .2*length(unique(deseq2trans$Gene)), limitsize=FALSE, heatmap )


deseq2trans = zscores

deseq2trans <- as.data.frame(deseq2trans)


#print(head(deseq2trans))

deseq2trans$Gene <- rownames(deseq2trans)
#print(head(deseq2trans))

#sigGenes <- currsampledata$genename[abs(currsampledata$currlogval) > 1.5 & currsampledata$currprob < .05]


#print(head(deseq2trans))

deseq2trans <- deseq2trans[deseq2trans$Gene %in% sigGenes,]
                                                                                                       
deseq2trans <- melt(deseq2trans, id.vars=c("Gene"))

samplenames = as.character(sampledata[,1])
#print(samplenames)
#print(samplenames[coldata$condition == currcomp[[1]] | coldata$condition == currcomp[[2]]])
currsamples = samplenames[coldata$condition == currcomp[[1]] | coldata$condition == currcomp[[2]]]
#print(currsamples)
#print(head(deseq2trans))
deseq2trans = deseq2trans[deseq2trans$variable %in% currsamples,]
#print("()")

maxscore = max(abs(deseq2trans$value))
# Make a heatmap
#print(head(deseq2trans))
heatmap <- ggplot(deseq2trans, aes(x=variable, y=Gene, fill=value)) + geom_raster() +  theme_bw() +  scale_fill_gradient2(high="darkred",low="darkblue", limits = c(-maxscore, maxscore))  + ggtitle(paste(currcomp[[1]], currcomp[[2]], sep = " vs "))+theme(axis.text.x=element_text(angle=65, hjust=1),  axis.ticks.y=element_blank()) + labs(fill = "Z score") #axis.text.y=element_blank(),  scale_fill_viridis(discrete=FALSE) scale_fill_distiller(palette = "RdBu")
print(paste(currcomp[[1]], currcomp[[2]], sep = " vs "))
print(1*length(unique(deseq2trans$Gene)))
ggsave(paste(experimentname,"/",pairname,"-deseqzscoreheatmap",".pdf",sep= ""),height=4+ .2*length(unique(deseq2trans$Gene)), limitsize=FALSE, heatmap )


}



displaygenes = c()
currsampledata$name = rownames(currsampledata)
displayfeats = ifelse(currsampledata$genename %in% displaygenes, as.character(currsampledata$genename), "")

pvalcutoff = sort(currsampledata$currprob)[10]

displayfeats = ifelse(abs(currsampledata$currlogval) > 1.5 & currsampledata$currprob < pvalcutoff, as.character(currsampledata$genename), "")

#print("**") 
#print(rownames(currsampledata))
#print(head(displayfeats))


#currsampledata = cbind(currlogval,currprob) # 
#print(head(currsampledata))
#currsampledata = currsampledata[currsampledata$currprob > .005,]
#print(head(currsampledata))
#print(currcomp[[1]])
#print(currcomp[[2]])
currplot <- ggplot(currsampledata, aes_string(x="currlogval", y="currprob")) + geom_point() +scale_x_continuous() +geom_text_repel(label = displayfeats,min.segment.length = unit(0, 'lines'), segment.color="red")+ scale_y_continuous(trans=reverselog_trans(10))+geom_hline(yintercept = .05, linetype = 2)+geom_hline(yintercept = .005, linetype = 2)+geom_vline(xintercept = dashinterc, linetype = 2) + geom_vline(xintercept = -dashinterc, linetype = 2)+theme_bw() + xlab("Log2-Fold Change")+ylab("Adjusted P-value")+ggtitle(currpair)+theme(legend.box="horizontal",aspect.ratio=1,axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+  labs(caption = c(currcomp[[1]], currcomp[[2]])) +  theme(plot.caption = element_text(size = 16,hjust=c(1, 0))) 


ggsave(paste(experimentname,"/",pairname ,"-volcano",outputformat,sep= ""), currplot) 


trnasampledata = currsampledata[(grepl( "tRNA", as.character(currsampledata$genename), fixed = TRUE) | grepl( "tRX", as.character(currsampledata$genename), fixed = TRUE)) & ( grepl( "fiveprime", as.character(currsampledata$genename)) | grepl( "threeprime", as.character(currsampledata$genename)) | grepl( "other", as.character(currsampledata$genename)) | grepl( "wholecounts", as.character(currsampledata$genename))) ,]

#print(trnasampledata)
trnapvalcutoff = sort(trnasampledata$currprob)[10]

trnadisplayfeats = ifelse(abs(trnasampledata$currlogval) > 1.5 & trnasampledata$currprob < trnapvalcutoff, as.character(trnasampledata$genename), "")


currplot <- ggplot(trnasampledata, aes_string(x="currlogval", y="currprob")) + geom_point() +scale_x_continuous() +  geom_text_repel(label = trnadisplayfeats,min.segment.length = unit(0, 'lines'), segment.color="red")+ scale_y_continuous(trans=reverselog_trans(10))+geom_hline(yintercept = .05, linetype = 2)+geom_hline(yintercept = .005, linetype = 2)+geom_vline(xintercept = dashinterc, linetype = 2) + geom_vline(xintercept = -dashinterc, linetype = 2)+theme_bw() + xlab("Log2-Fold Change")+ylab("Adjusted P-value")+ggtitle(paste(currpair,"_tRNAs",sep = ""))+theme(legend.box="horizontal",aspect.ratio=1,axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+  labs(caption = c(currcomp[[1]], currcomp[[2]])) +  theme(plot.caption = element_text(size = 16,hjust=c(1, 0))) 
ggsave(paste(experimentname,"/",pairname ,"-volcano_tRNA",outputformat,sep= ""), currplot) 

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


#print(medcounts) 
medcountmat <- do.call("cbind",medcounts)

colnames(medcountmat) <- names(medcounts)


write.table(medcountmat,paste(experimentname,"/",experimentname,"-medians.txt", sep = ""),quote=FALSE)


#print(head(medcountmat))
medcountmat = as.matrix(medcountmat)


allcombinevals = as.matrix(allcombinevals)

#medcounts
#typeof(allcombinevals)
#typeof(medcountmat)
#typeof(medcounts)
#typeof(normalizedrnas)
allcombinevals = cbind(allcombinevals,medcountmat)

#write.table(allcombinevals[apply(readcounts,1,max) > 30 & apply(alllogvals,1,min) < .05 ,],paste(experimentname,"/",experimentname,"-relevnormalizedsamples.txt", sep = ""), col.names=NA )

write.table(allcombinevals,paste(experimentname,"/",experimentname,"-combine.txt", sep = ""), col.names=NA )
