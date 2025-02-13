library(ggplot2)
library(gridExtra)
library(scales)
library(plyr)
library(RColorBrewer)
library(ggrepel)




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
counts <- read.table(args[2], stringsAsFactors = FALSE)
trnatable <- read.table(args[3], stringsAsFactors = FALSE)

types <- read.table(args[4])

sampletable <- read.table(args[5], stringsAsFactors = FALSE)

samplenames = unique(sampletable[,2])
comparisons <- read.table(args[6], stringsAsFactors = FALSE)

counts = counts + 1



counts <- merge(types,counts, by.x = 1, by.y = 0)
colnames(counts)[2] <- "type"
colnames(counts)[3] <- "chrom"
colnames(counts)[4] <- "readsize"


genetypes = c("trna_fiveprime","trna_threeprime","trna_other","trna_wholecounts","tRNA","snoRNA","miRNA","Mt_tRNA","rRNA","Mt_rRNA","snRNA")
othertypes = !(counts[,"type"] %in% genetypes)


counts[,"type"] <- as.character(counts[,"type"])
#counts[othertypes,"type"] <- "other"
counts[,"type"] <- as.factor(counts[,"type"])


othertypenames = unique(as.character(counts[othertypes,"type"]))
#print(othertypenames)
                                                                                                                                                       
#print(unique(as.character(counts[,"type"])))
#print(unique(as.character(othertypenames)))
genetypes
#print(unique(as.character(counts[,"type"])))
counts[,"type"] <- factor(counts[,"type"], levels = c(genetypes,"other",othertypenames))
#print(unique(as.character(counts[,"type"])))
trnatypes <- c("trna_fiveprime","trna_threeprime","trna_other","trna_wholecounts","tRNA")
fragtypes <- c("trna_fiveprime","trna_threeprime","trna_other","trna_wholecounts")
                                                                                                                                                       
#print(unique(as.character(counts[,"type"])))
trnagenes =  counts[,"type"] %in% trnatypes


#print("**||**")

#counts <- rbind(counts[counts[,"type"] == "other",],counts[counts[,"type"] == "snoRNA",],counts[counts[,"type"] == "miRNA",], counts[counts[,"type"] == "Mt_tRNA",],counts[counts[,"type"] == "rRNA",],counts[counts[,"type"] == "Mt_rRNA",],counts[counts[,"type"] == "snRNA",],counts[trnagenes,])
                                                                                                                                                       
#print(unique(as.character(counts[,"type"])))

colnames(counts)[1] <- "name"
colnames(trnatable) <- c("trnaname", "loci", "amino", "anticodon")

counts$trnaname = ifelse( counts[,"type"] %in%  fragtypes,sub("^(.*)_[^_]*$", "\\1", as.character(counts[,"name"])),as.character(counts[,"name"]))
#head(counts[counts[,"type"] %in% trnatypes,])
trnacounts <- merge(trnatable,counts, by.x = "trnaname", by.y = "trnaname", all.y=TRUE)


i = 1
#trnacounts[is.na(trnacounts$amino),"amino"] = "non-tRNA"
#trnacounts <- rbind(trnacounts[trnacounts[,"type"] == "other",],trnacounts[trnacounts[,"type"] == "snoRNA",],trnacounts[trnacounts[,"type"] == "miRNA",], trnacounts[trnacounts[,"type"] == "Mt_tRNA",],trnacounts[trnacounts[,"type"] == "rRNA",],trnacounts[trnacounts[,"type"] == "Mt_rRNA",],trnacounts[trnacounts[,"type"] == "snRNA",],trnacounts[trnagenes,])

#trnacounts <- rbind(trnacounts[trnacounts[,"type"] == "other",],trnacounts[trnacounts[,"type"] == "snoRNA",],trnacounts[trnacounts[,"type"] == "miRNA",], trnacounts[trnacounts[,"type"] %in% trnatypes,])
trnacounts$dotsize = ifelse(trnacounts[,"type"] %in% trnatypes, trnasize, .1)

#cor.test(log(counts[,xaxis]+1),log(counts[,yaxis]+1))

onlytrnas <- trnacounts[trnacounts[,"type"] %in% trnatypes,]

#trnacounts[,"type"]

trnacounts$fragtype <- factor(ifelse(as.character(trnacounts[,"type"]) %in% trnatypes, as.character(trnacounts[,"type"]), "nontRNA" ), levels = c("nontRNA",trnatypes))


trnacounts$amino <- as.factor(trnacounts$amino)
#levels(trnacounts$amino) <- sort(trnacounts$amino)

#script.dir <- dirname(sys.frame(1)$ofile)

#frame_files <- lapply(sys.frames(), function(x) x$ofile)
#frame_files <- Filter(Negate(is.null), frame_files)
#script.dir <- dirname(frame_files[[length(frame_files)]])


#aminoinfo <- read.table(paste(script.dir, "aminotable.txt", sep=''), stringsAsFactors = FALSE)

aminos = c(" Glycine","Proline","Alanine","Valine","Leucine","Isoleucine","Methionine","Cysteine","Phenylalanine","Tyrosine","Tryptophan","Histidine","Lysine","Arginine","Glutamine","Asparagine","Glutamic_Acid","Aspartic_Acid","Serine","Threonine","iMethionine","Selenocysteine")
threecodes = c("Gly","Pro","Ala","Val","Leu","Ile","Met","Cys","Phe","Tyr","Trp","His","Lys","Arg","Gln","Asn","Glu","Asp","Ser","Thr","iMet", "SeC")
onecodes = c("G","P","A","V","L","I","M","C","F","Y","W","H","K","R","Q","N","E","D","S","T","X", "Z")

aminoinfo = data.frame(aminos,threecodes,onecodes, stringsAsFactors = FALSE)

aminoletters <- aminoinfo[match(aminoinfo[,2], levels(trnacounts$amino)),3]
aminoletters <- aminoinfo[match(levels(trnacounts$amino),aminoinfo[,2]),3]

aminoletters[is.na(aminoletters)] <- "?"


aminoletters <-  unlist(lapply(aminoletters, utf8ToInt))


trnacounts$fragtype <- mapvalues(trnacounts$fragtype, from = c("trna_threeprime", "trna_fiveprime", "trna_other", "trna_wholecounts"), to = c("Three-prime fragments","Five-prime fragments","Other fragments","Whole tRNAs"))
#remove non-tRNAs
typepal <- c(
  "other" = "#f2dab2",
  "tRNA" = "#6666cc", 
  "pretRNA" = "#a1ade5", 
  "miRNA" = "#b4008d",
  "snRNA" = "#74531d",
  "Mt_tRNA" = "#87e275",
  "Mt_rRNA" = "#00754a",
  "rRNA" = "#a5e5d9",
  "snoRNA" = "#ff9e18",
  "misc_RNA" = "#b2b2b2",
  "trna_fiveprime" = "#f3b2db",
  "trna_threeprime" = "#fd495c",
  "trna_other" = "#000f9f",
  "trna_wholecounts" = "#6666cc"
  
)


if(FALSE)
{
trnacounts = trnacounts[trnacounts$fragtype != "nontRNA",]
}
#samplenames
#samplenames <- unique(sampletable[,2])
for (i in 1:length(samplenames)){
cols <- sampletable[sampletable[,2] == samplenames[i],1]
if (length(cols) > 1){
trnacounts[,samplenames[i]] <- apply(trnacounts[,cols], 1, median)

}else{
trnacounts[,samplenames[i]] <- trnacounts[,cols]
}
}



for (i in 1:length(rownames(comparisons))){

yaxis = sampletable[sampletable[,2] == comparisons[i,1],1][1]
xaxis = sampletable[sampletable[,2] == comparisons[i,2],1][1]
#sampletable[sampletable[,2] == comparisons[i,1],1]
#print(sampletable[sampletable[,2] == comparisons[i,2],1])


yname = comparisons[i,1]
xname = comparisons[i,2]
xaxis = xname
yaxis = yname



#print(xaxis)
#print(yaxis)

#maxlim = max(max(log(trnacounts[,xaxis])), max(log(trnacounts[,yaxis]))) #, xlim = c(0,maxlim), ylim = c(0,maxlim)

#here is error
#"Untreated"
#print(trnacounts)
#print (xaxis)
maxlim = max(c(max(trnacounts[,xaxis]), max(trnacounts[,yaxis]))) #, xlim = c(0,maxlim), ylim = c(0,maxlim) + 1


corr = cor.test(log(trnacounts[,xaxis]+1),log(trnacounts[,yaxis]+1))

sublabel = paste("Pearson Correlation: ",corr$estimate, sep = "")  

dashinterc = 1.5

colourCount = length(unique(counts$type))+1
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
#scale_colour_manual(values = getPalette(colourCount)) 
#print("*||**")
#Snora35, and SNORD116 (there are a few from a big gene cluster) let7a,b,c, mir-138, mir-122, mir-133a
displaygenes = c()
#displaygenes = c("Snora35","SNORD116", "let7a" ,"let7b","let7c","mir-138", "mir-122", "mir-133a")

#displaygenes = c("Snora35","Snord116l17", "Mirlet7a-2" ,"Mirlet7b","Mirlet7c-2","Mir138-1", "Mir122", "Mir133a-1")
#livergenes =  c("Mir33-201","Mir223","Mir30c-1","Mir144","Mir148a","Mir24-1","Mir29a","Mir122")
#musclegenes = c("Mir1a-1","Mir133a-1","Mir208a","Mir208b","Mir499")
#displaygenes = c(displaygenes, livergenes, musclegenes)

 # miR-33, -33*, miR-223, -30c, -144, -148a, -24, -29, and -122
 #miR-1, miR133a, miR-208a/b, and miR-499
 


#print(trnacounts[trnacounts$name %in% displaygenes,])
displayfeats = ifelse(trnacounts$name %in% displaygenes, as.character(trnacounts$name), "")
#print(trnacounts[trnacounts$name %in% displaygenes,"name"])
#print(unique(displayfeats))
#print(unique(trnacounts$type))


extratypes = setdiff(unique(trnacounts$type), names(typepal))
extratypes = sort(extratypes)
#print(extratypes)
#print(unique(countsmelt$seq))
#print(gg_color_hue(length(extratypes)))
extracolors = setNames(gg_color_hue(length(extratypes)), extratypes)
typepal = c(typepal, extracolors)
print("***")

xaxisgg = enquo(xaxis)
yaxisgg = enquo(yaxis)
currplot <- ggplot(trnacounts, aes(x=!!ensym(xaxis), y=!!ensym(yaxis))) + geom_point(aes(color=type), size = 1) + 
    #scale_colour_manual(values = getPalette(colourCount)) + 
    scale_colour_manual(values = typepal)+
    geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = dashinterc, slope = 1,linetype = 2)+geom_abline(intercept = 0- dashinterc, slope = 1,linetype = 2)+
    scale_x_continuous(trans=log2_trans(),limits = c(1, maxlim)) + scale_y_continuous(trans=log2_trans(),limits = c(1, maxlim)) + 
    theme_bw() + 
    #geom_text_repel(label = displayfeats,min.segment.length = unit(0, 'lines'))+
    theme(legend.box="horizontal",aspect.ratio=1,axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 
    
#currplot <- ggplot(data=counts,x=trnacounts[,xaxis],y=trnacounts[,yaxis],xlab = xaxis,ylab = yaxis,color=type, asp=1) + scale_colour_manual(values = getPalette(colourCount)) + geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = dashinterc, slope = 1,linetype = 2)+geom_abline(intercept = 0- dashinterc, slope = 1,linetype = 2)+scale_x_continuous(trans=log2_trans(),limits = c(1, maxlim)) + scale_y_continuous(trans=log2_trans(),limits = c(1, maxlim)) + theme_bw() + theme(legend.box="horizontal",aspect.ratio=1,axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))



#currplot <- arrangeGrob(currplot, sub = textGrob(sublabel, x = 0, hjust = -0.1, vjust=0.1, gp = gpar(fontsize = 14)))
ggsave(paste(experimentname,"/",comparisons[i,1],"_",comparisons[i,2] ,"-typescatter",outputformat,sep= ""), currplot, width = 8, height = 8)
currplot = currplot + facet_wrap(~type, ncol = 3)
ggsave(paste(experimentname,"/",comparisons[i,1],"_",comparisons[i,2] ,"-typescatterfacet",outputformat,sep= ""), currplot, width = 14, height = length(unique(trnacounts$type)))
#print("***||")

trnacounts$logfc = log2(trnacounts[,xaxis]) -  log2(trnacounts[,yaxis])
#print(head(trnacounts))
#print(max(abs(trnacounts$logfc)))
currplot <- ggplot(trnacounts, aes(x=logfc, y=readsize)) + geom_point(aes(color=type), size = 1) + 
    #scale_colour_manual(values = getPalette(colourCount)) + 
    ggtitle(paste(xaxis, " vs ",yaxis, sep = "")) +
    scale_colour_manual(values = typepal)+
    #geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = dashinterc, slope = 1,linetype = 2)+geom_abline(intercept = 0- dashinterc, slope = 1,linetype = 2)+
    scale_x_continuous(limits = c(-max(abs(trnacounts$logfc)), max(abs(trnacounts$logfc)))) +
    theme_bw() + 
    #geom_text_repel(label = displayfeats,min.segment.length = unit(0, 'lines'))+
    theme(legend.box="horizontal",aspect.ratio=1,axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 
    
#currplot <- ggplot(data=counts,x=trnacounts[,xaxis],y=trnacounts[,yaxis],xlab = xaxis,ylab = yaxis,color=type, asp=1) + scale_colour_manual(values = getPalette(colourCount)) + geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = dashinterc, slope = 1,linetype = 2)+geom_abline(intercept = 0- dashinterc, slope = 1,linetype = 2)+scale_x_continuous(trans=log2_trans(),limits = c(1, maxlim)) + scale_y_continuous(trans=log2_trans(),limits = c(1, maxlim)) + theme_bw() + theme(legend.box="horizontal",aspect.ratio=1,axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))



#currplot <- arrangeGrob(currplot, sub = textGrob(sublabel, x = 0, hjust = -0.1, vjust=0.1, gp = gpar(fontsize = 14)))
#ggsave(paste(experimentname,"/",comparisons[i,1],"_",comparisons[i,2] ,"-lengthcompare",outputformat,sep= ""), currplot)

trnacounts$mincounts = min(log2(trnacounts[,xaxis]),log2(trnacounts[,yaxis]))
currplot <- ggplot(trnacounts, aes(x=logfc, y=mincounts)) + geom_point(aes(color=type), size = 1) + 
    #scale_colour_manual(values = getPalette(colourCount)) + 
    ggtitle(paste(xaxis, " vs ",yaxis, sep = "")) +
    scale_colour_manual(values = typepal)+
    #geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = dashinterc, slope = 1,linetype = 2)+geom_abline(intercept = 0- dashinterc, slope = 1,linetype = 2)+
    scale_x_continuous(limits = c(-max(abs(trnacounts$logfc)), max(abs(trnacounts$logfc)))) +
    theme_bw() + 
    #geom_text_repel(label = displayfeats,min.segment.length = unit(0, 'lines'))+
    theme(legend.box="horizontal",aspect.ratio=1,axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 
    
#currplot <- ggplot(data=counts,x=trnacounts[,xaxis],y=trnacounts[,yaxis],xlab = xaxis,ylab = yaxis,color=type, asp=1) + scale_colour_manual(values = getPalette(colourCount)) + geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = dashinterc, slope = 1,linetype = 2)+geom_abline(intercept = 0- dashinterc, slope = 1,linetype = 2)+scale_x_continuous(trans=log2_trans(),limits = c(1, maxlim)) + scale_y_continuous(trans=log2_trans(),limits = c(1, maxlim)) + theme_bw() + theme(legend.box="horizontal",aspect.ratio=1,axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))



#currplot <- arrangeGrob(currplot, sub = textGrob(sublabel, x = 0, hjust = -0.1, vjust=0.1, gp = gpar(fontsize = 14)))
#ggsave(paste(experimentname,"/",comparisons[i,1],"_",comparisons[i,2] ,"-countcompare",outputformat,sep= ""), currplot)

#+scale_shape_manual(values = c(20,15, 17,18,19)) theme(legend.position = "bottom")
#currplot <- ggplot(trnacounts, aes(x=!!ensym((xaxis), y=!!ensym((yaxis)))+geom_point(aes(color=amino, size = dotsize, shape=fragtype))+guides(size=FALSE, ncol = 1)+scale_size_continuous(range = c(.75,2))+geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = dashinterc, slope = 1,linetype = 2)+geom_abline(intercept = 0- dashinterc, slope = 1,linetype = 2)+scale_x_continuous(trans=log2_trans(),limits = c(1, maxlim)) + scale_y_continuous(trans=log2_trans(),limits = c(1, maxlim)) + theme_bw() + theme(legend.box="horizontal",aspect.ratio=1,axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

#scale_shape_manual(values = aminoletters)
#scale_shape_manual(values=1:nlevels(trnacounts$amino))
#print("***DONESCATTER")

#print(head(trnacounts)) #PNK +geom_point(data = transform(trnacounts[trnacounts$fragtype == "nontRNA",], fragtype=NULL), aes(size = dotsize))
currplot <- ggplot(trnacounts[trnacounts$fragtype != "nontRNA",], aes(x=!!ensym(xaxis), y=!!ensym(yaxis)))+xlab(gsub("_", " ", xname))+ylab(gsub("_", " ", yname))+geom_point(aes(shape=amino, color=amino, size = dotsize))+guides(size=FALSE, ncol = 1)+scale_shape_manual(values=aminoletters) +facet_wrap( ~fragtype , ncol = 2) + scale_size_continuous(range = c(.25,4))+geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = dashinterc, slope = 1,linetype = 2)+geom_abline(intercept = 0- dashinterc, slope = 1,linetype = 2)+scale_x_continuous(trans=log2_trans(),limits = c(1, maxlim),breaks = trans_breaks('log2', function(x) 2^x, n = 10)) + scale_y_continuous(trans=log2_trans(),limits = c(1, maxlim),breaks= trans_breaks('log2', function(x) 2^x, n = 10)) + theme_bw() + theme(legend.box="horizontal",aspect.ratio=1,axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + labs(shape="Acceptor\nType", color="Acceptor\nType") 
ggsave(paste(experimentname,"/",comparisons[i,1],"_",comparisons[i,2] ,"-aminoscatter",outputformat,sep= ""), currplot, height = 10, width = 12)

#currplot <- qplot(data=onlytrnas,x=onlytrnas[,xaxis],y=onlytrnas[,yaxis],xlab = xaxis,ylab = yaxis, asp=1)+geom_abline(intercept = 0, slope = 1) + geom_abline(intercept = dashinterc, slope = 1,linetype = 2)+geom_abline(intercept = 0- dashinterc, slope = 1,linetype = 2)+scale_x_continuous(trans=log2_trans(),limits = c(1, maxlim)) + scale_y_continuous(trans=log2_trans(),limits = c(1, maxlim)) +facet_wrap( ~amino, , ncol = 2)
#ggsave(paste(experimentname,"/",comparisons[i,1],"_",comparisons[i,2] ,"-tRNAscatter.pdf",sep= ""), currplot,height=4*length(unique(trnacounts$amino)),width=20, limitsize=FALSE)




}