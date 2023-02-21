library(ggplot2)
library(gridExtra)
library(scales)
library(plyr)
library(RColorBrewer)
library(ggrepel)


#Rscript /projects/lowelab/users/holmes/pythonsource/TRAX//analyzeunique.R ottrctrlall ottrctrlall/unique/ottrctrlall-trnauniquecounts.txt ottrctrlall/ottrctrlall-anticodoncounts.txt ottrctrlall/ottrctrlall-aminocounts.txt ottrctrlall/ottrctrlall-SizeFactors.txt /soe/holmes/pythonsource/trnatest/trnadbs/mm10/mm10-trnatable.txt  ottrctrlsamples.txt ottrctrlpairs.txt


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

args <- commandArgs(trailingOnly = TRUE)

#args[3]
#args[4]


log2_minor_break = function (...){
  function(x) {
    minx         = floor(min(log2(x), na.rm=T))-1;
    maxx         = ceiling(max(log2(x), na.rm=T))+1;
    n_major      = maxx-minx+1;
    major_breaks = seq(minx, maxx, by=1)
    minor_breaks = 
      rep(log2(seq(1, 9, by=1)), times = n_major)+
      rep(major_breaks, each = 9)
    return(2^(minor_breaks))
  }
}

outputformat <- ".pdf"
#uniqueInitials <- c("a", "A", "c", "C", "j", "J", "R")


trnasize = 6

experimentname <- args[1]
trnauniquecounts <- read.table(args[2], stringsAsFactors = FALSE)
anticodoncounts <- read.table(args[3], stringsAsFactors = FALSE)
aminocounts <- read.table(args[4], stringsAsFactors = FALSE)
sizefactors <- read.table(args[5], stringsAsFactors = FALSE,header= TRUE)

trnatable <- read.table(args[6], stringsAsFactors = FALSE)
sampletable <- read.table(args[7], stringsAsFactors = FALSE)

haspairs = !is.na(args[8])
if(haspairs){

comparisons <- read.table(args[8], stringsAsFactors = FALSE)
}
sizefactorvec = sizefactors[1,]





#trnauniquecounts = sweep(trnauniquecounts, 2, sizefactorvec, FUN = '/')
#anticodoncounts = sweep(anticodoncounts, 2, sizefactorvec, FUN = '/')
#aminocounts = sweep(aminocounts, 2,sizefactorvec , FUN = '/')

trnauniquemat     =   mapply('/',trnauniquecounts,  sizefactorvec)
anticodoncountmat =   mapply('/',anticodoncounts,  sizefactorvec)
aminocountmat     =   mapply('/',aminocounts,  sizefactorvec)

rownames(trnauniquemat)     = rownames(trnauniquecounts)
rownames(anticodoncountmat) = rownames(anticodoncounts)
rownames(aminocountmat)     = rownames(aminocounts)



origtrnauniquecounts = data.frame(trnauniquemat)     
origanticodoncounts  = data.frame(anticodoncountmat) 
origaminocounts      = data.frame(aminocountmat)     

write.table(origtrnauniquecounts, file =  paste(experimentname,"/unique/",experimentname ,"-uniquetrnanormcounts.txt",sep= ""))
write.table(origanticodoncounts,  file = paste(experimentname,"/unique/",experimentname ,"-uniqueacnormcounts.txt",sep= ""))
write.table(origaminocounts,      file = paste(experimentname,"/unique/",experimentname ,"-uniqueaminonormcounts.txt",sep= ""))



samplenames = unique(sampletable[,2])

trnauniquecounts = origtrnauniquecounts + 1
anticodoncounts = origanticodoncounts + 1
aminocounts = origaminocounts + 1
#print(head(trnauniquecounts))
#print(head(aminocounts))

#print("**||")
#q()

#counts <- merge(types,counts, by.x = 1, by.y = 0)
#colnames(counts)[2] <- "type"
#colnames(counts)[3] <- "chrom"
#colnames(counts)[4] <- "readsize"




genetypes = c("trna_fiveprime","trna_threeprime","trna_other","trna_wholecounts","tRNA","snoRNA","miRNA","Mt_tRNA","rRNA","Mt_rRNA","snRNA")
#othertypes = !(counts[,"type"] %in% genetypes)


#counts[,"type"] <- as.character(counts[,"type"])
#counts[,"type"] <- as.factor(counts[,"type"])



#counts[,"type"] <- factor(counts[,"type"], levels = c(genetypes,"other",othertypenames))
trnatypes <- c("trna_fiveprime","trna_threeprime","trna_other","trna_wholecounts","tRNA")
fragtypes <- c("trna_fiveprime","trna_threeprime","trna_other","trna_wholecounts")
                                                                                                                                                       

#colnames(counts)[1] <- "name"
colnames(trnatable) <- c("trnaname", "loci", "amino", "anticodon")

#counts$trnaname = ifelse( counts[,"type"] %in%  fragtypes,sub("^(.*)_[^_]*$", "\\1", as.character(counts[,"name"])),as.character(counts[,"name"]))


trnauniquecounts$name = rownames(trnauniquecounts)
anticodoncounts$name = rownames(anticodoncounts)
aminocounts$name = rownames(aminocounts)

trnauniquecounts$trnaname = sub("^(.*)_[^_]*$", "\\1", as.character(trnauniquecounts[,"name"]))
anticodoncounts$acname =  sub("^(.*)_[^_]*$", "\\1", as.character(anticodoncounts[,"name"]))
aminocounts$aminoname =      sub("^(.*)_[^_]*$", "\\1", as.character(aminocounts[,"name"]))


trnauniquecounts$fragtype = sub("^.*_([^_]*)$", "\\1", as.character(trnauniquecounts[,"name"]))
anticodoncounts$fragtype =  sub("^.*_([^_]*)$", "\\1", as.character(anticodoncounts[,"name"]))
aminocounts$fragtype =      sub("^.*_([^_]*)$", "\\1", as.character(aminocounts[,"name"]))


trnauniquecounts <- merge(trnatable,trnauniquecounts, by.x = "trnaname", by.y = "trnaname", all.y=TRUE)

#print(head(trnauniquecounts))






i = 1


aminos = c(" Glycine","Proline","Alanine","Valine","Leucine","Isoleucine","Methionine","Cysteine","Phenylalanine","Tyrosine","Tryptophan","Histidine","Lysine","Arginine","Glutamine","Asparagine","Glutamic_Acid","Aspartic_Acid","Serine","Threonine","iMethionine","Selenocysteine")
threecodes = c("Gly","Pro","Ala","Val","Leu","Ile","Met","Cys","Phe","Tyr","Trp","His","Lys","Arg","Gln","Asn","Glu","Asp","Ser","Thr","iMet", "SeC")
onecodes = c("G","P","A","V","L","I","M","C","F","Y","W","H","K","R","Q","N","E","D","S","T","X", "Z")





aminoinfo = data.frame(aminos,threecodes,onecodes, stringsAsFactors = FALSE)

#print(trnauniquecounts$amino)


aminoletters <- aminoinfo[match(unique(trnauniquecounts$amino),aminoinfo[,2]),3]
#print(aminoletters)
aminoletters[is.na(aminoletters)] <- "?"


aminoletters <-  unlist(lapply(aminoletters, utf8ToInt))

#q()

amino = c("Ala","Ala","Ala","Arg","Arg","Arg","Arg","Arg","Asn","Asp","Cys","Gln","Gln","Glu","Glu","Gly","Gly","Gly","Gly","His","Ile","Ile","Leu","Leu","Leu","Leu","Leu","Lys","Lys","Met","Phe","Pro","Pro","Pro","SeC","Ser","Ser","Ser","Ser","Ser","Sup","Sup","Thr","Thr","Thr","Trp","Tyr","Val","Val","Val","iMet")
anticodon = c("AGC","CGC","TGC","ACG","CCG","CCT","TCG","TCT","GTT","GTC","GCA","CTG","TTG","CTC","TTC","ACC","CCC","GCC","TCC","GTG","AAT","TAT","AAG","CAA","CAG","TAA","TAG","CTT","TTT","CAT","GAA","AGG","CGG","TGG","TCA","AGA","CGA","GCT","GGA","TGA","TCA","TTA","AGT","CGT","TGT","CCA","GTA","AAC","CAC","TAC","CAT")
anticodoninfo = data.frame(amino,anticodon, stringsAsFactors = FALSE)




anticodoncounts <- merge(anticodoninfo,anticodoncounts, by.x = "anticodon", by.y = "acname", all.y=TRUE)


#print(head(anticodoncounts))
#q()

trnauniquecounts$fragtype <- mapvalues(trnauniquecounts$fragtype, from = c("threeprime", "fiveprime", "other", "whole"), to = c("Three-prime fragments","Five-prime fragments","Other fragments","Whole tRNAs"))
anticodoncounts$fragtype <- mapvalues(anticodoncounts$fragtype, from = c("Threeprime", "Fiveprime", "Other", "Whole"), to = c("Three-prime fragments","Five-prime fragments","Other fragments","Whole tRNAs"))
aminocounts$fragtype <- mapvalues(aminocounts$fragtype, from = c("Threeprime", "Fiveprime", "Other", "Whole"), to = c("Three-prime fragments","Five-prime fragments","Other fragments","Whole tRNAs"))







#remove non-tRNAs


for (i in 1:length(samplenames)){
cols <- sampletable[sampletable[,2] == samplenames[i],1]
if (length(cols) > 1){
#trnacounts[,samplenames[i]] <- apply(trnacounts[,cols], 1, median)

trnauniquecounts[,samplenames[i]] <- apply(trnauniquecounts[,cols], 1, median)
anticodoncounts[,samplenames[i]] <- apply(anticodoncounts[,cols], 1, median)
aminocounts[,samplenames[i]] <- apply(aminocounts[,cols], 1, median)

}else{
#trnacounts[,samplenames[i]] <- trnacounts[,cols]
trnauniquecounts[,samplenames[i]] <- trnauniquecounts[,cols]
anticodoncounts[,samplenames[i]] <- anticodoncounts[,cols]
aminocounts[,samplenames[i]] <- aminocounts[,cols]


}
}

aminocounts$amino = aminocounts$aminoname



data <- origanticodoncounts[rowSums(origanticodoncounts) > 20,]
datapca <- prcomp(t(data),center = TRUE,scale = TRUE) 


percentlabels <- round(datapca$sdev / sum(datapca$sdev) * 100, 2)
percentlabels <- paste( colnames(datapca$x), "(", paste( as.character(percentlabels), "%", ")", sep="") )
#


scores = as.data.frame(datapca$x)

print(head(scores))
# plot of observations

#print(sampledata[match(rownames(scores),sampledata[,1]),2])
#print(rownames(scores))
#print(sampledata)
#aes(colour = factor(cyl))

samplename = sampletable[match(rownames(scores),sampletable[,1]),2]

ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores), color = samplename)) + theme_bw()+labs(color="Sample Name")+ geom_point()+geom_hline(yintercept = 0, colour = "gray65") +geom_vline(xintercept = 0, colour = "gray65") + geom_text(alpha = 0.8, size = 2,vjust="inward",hjust="inward") + ggtitle("Anticodon Read Principle Component Analysis")    + xlab(percentlabels[1]) +    ylab(percentlabels[2]) 
ggsave(filename=paste(experimentname,"/unique/",experimentname,"-anticodonpca",outputformat,sep= ""))


if(haspairs){
for (i in 1:length(rownames(comparisons))){

yaxis = sampletable[sampletable[,2] == comparisons[i,1],1][1]
xaxis = sampletable[sampletable[,2] == comparisons[i,2],1][1]

yname = comparisons[i,1]
xname = comparisons[i,2]
xaxis = xname
yaxis = yname



maxtrnalim = max(c(max(trnauniquecounts[,xaxis]), max(trnauniquecounts[,yaxis])))
maxaminolim = max(c(max(anticodoncounts[,xaxis]), max(anticodoncounts[,yaxis])))
maxaclim = max(c(max(aminocounts[,xaxis]), max(aminocounts[,yaxis]))) 


#corr = cor.test(log(trnacounts[,xaxis]+1),log(trnacounts[,yaxis]+1))
 

dashinterc = 1.5
print(head(aminocounts))
#print(head(anticodoncounts))

#print(head(maxaminolim))
#print(head(aminoletters))

currplot <- ggplot(trnauniquecounts, aes_string(x=xaxis, y=yaxis))+
    ggtitle("Uniquely tRNA mapping  reads")+
    xlab(gsub("_", " ", xname))+ylab(gsub("_", " ", yname))+
    geom_point(aes(shape=amino, color=amino, size = .5))+guides(size=FALSE, ncol = 1)+
    
    scale_shape_manual(values=aminoletters) +facet_wrap( ~fragtype , ncol = 2) + 
    
    scale_size_continuous(range = c(.25,4))+geom_abline(intercept = 0, slope = 1) + 
    geom_abline(intercept = dashinterc, slope = 1,linetype = 2)+
    geom_abline(intercept = 0- dashinterc, slope = 1,linetype = 2)+
    scale_x_continuous(trans=log2_trans(),limits = c(1, maxaminolim),breaks = trans_breaks('log2', function(x) 2^x, n = 10)) + 
    scale_y_continuous(trans=log2_trans(),limits = c(1, maxaminolim),breaks= trans_breaks('log2', function(x) 2^x, n = 10)) + 
    theme_bw() + 
    theme(legend.box="horizontal",aspect.ratio=1,axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
    labs(shape="Acceptor\nType", color="Acceptor\nType") 
ggsave(paste(experimentname,"/unique/",comparisons[i,1],"_",comparisons[i,2] ,"-uniquetrnascatter",outputformat,sep= ""), currplot, height = 10, width = 12)

currplot <- ggplot(anticodoncounts, aes_string(x=xaxis, y=yaxis))+
    ggtitle("anticodon mapping  reads")+
    xlab(gsub("_", " ", xname))+ylab(gsub("_", " ", yname))+
    geom_point(aes(shape=amino, color=amino, size = .5))+guides(size=FALSE, ncol = 1)+
    
    scale_shape_manual(values=aminoletters) +facet_wrap( ~fragtype , ncol = 2) + 
    
    scale_size_continuous(range = c(.25,4))+geom_abline(intercept = 0, slope = 1) + 
    geom_abline(intercept = dashinterc, slope = 1,linetype = 2)+
    geom_abline(intercept = 0- dashinterc, slope = 1,linetype = 2)+
    scale_x_continuous(trans=log2_trans(),limits = c(1, maxaminolim),breaks = trans_breaks('log2', function(x) 2^x, n = 10)) + 
    scale_y_continuous(trans=log2_trans(),limits = c(1, maxaminolim),breaks= trans_breaks('log2', function(x) 2^x, n = 10)) + 
    theme_bw() + 
    theme(legend.box="horizontal",aspect.ratio=1,axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
    labs(shape="Acceptor\nType", color="Acceptor\nType") 
ggsave(paste(experimentname,"/unique/",comparisons[i,1],"_",comparisons[i,2] ,"-uniqueanticodonscatter",outputformat,sep= ""), currplot, height = 10, width = 12)


currplot <- ggplot(aminocounts, aes_string(x=xaxis, y=yaxis))+
    ggtitle("amino mapping reads")+
    xlab(gsub("_", " ", xname))+ylab(gsub("_", " ", yname))+
    geom_point(aes(shape=amino, color=amino, size = .5))+guides(size=FALSE, ncol = 1)+
    
    scale_shape_manual(values=aminoletters) +facet_wrap( ~fragtype , ncol = 2) + 
    
    scale_size_continuous(range = c(.25,4))+geom_abline(intercept = 0, slope = 1) + 
    geom_abline(intercept = dashinterc, slope = 1,linetype = 2)+
    geom_abline(intercept = 0- dashinterc, slope = 1,linetype = 2)+
    scale_x_continuous(trans=log2_trans(),limits = c(1, maxaminolim),breaks = trans_breaks('log2', function(x) 2^x, n = 10)) + 
    scale_y_continuous(trans=log2_trans(),limits = c(1, maxaminolim),breaks= trans_breaks('log2', function(x) 2^x, n = 10)) + 
    theme_bw() + 
    theme(legend.box="horizontal",aspect.ratio=1,axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
    labs(shape="Acceptor\nType", color="Acceptor\nType") 
ggsave(paste(experimentname,"/unique/",comparisons[i,1],"_",comparisons[i,2] ,"-uniqueaminoscatter",outputformat,sep= ""), currplot, height = 10, width = 12)




}
}