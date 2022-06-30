library(ggplot2)
library(RColorBrewer)

library(reshape2)
library(scales)
library(plyr)
library(gridExtra)
library(getopt)


args <- commandArgs(trailingOnly = TRUE)
#mousearmseqpaper3-.txt
#Rscript /projects/lowelab/users/holmes/pythonsource/TRAX/endplot.R --ends=mousepnksamples2/mousepnksamples2-trnaendcounts.txt --trna=/projects/lowelab/users/holmes/pythonsource/trnatest/trnadbs/mm10/mm10-trnatable.txt --samples=mousepnksamples.txt --directory=mousepnksamples2/mismatch/

#hcvm.i <- melt(hcvm.i, id.vars=c(grep("^readC", names(hcvm.i), value=TRUE, invert=TRUE)), variable.name="tRNA.basePosition", value.name="read.density")


spec <- matrix(c(
        'ends'     , 'n', 1, "character", "trna end file (required)",
        'samples'    , 's', 1, "character", "Sample file (required)",
        'trna'     , 't', 1, "character", "trna file (required)",
        'directory'     , 'd', 1, "character", "output directory (required)",
        'runname'     , 'r', 1, "character", "run name",        
        'help'   , 'h', 0, "logical",   "this help"
),ncol=5,byrow=T)

opt = getopt(spec);

dotsize = .4
aminodotsize = .8
ccapseudocount = .1


ends <- read.table(opt$ends, header = TRUE,row.names = NULL, stringsAsFactors=FALSE)
Sampletable <- read.table(opt$samples)
trnatable <- read.table(opt$trna)                                    
runname = opt$runname
directory <- opt$directory



outputformat <- ".pdf"
aminos = c(" Glycine","Proline","Alanine","Valine","Leucine","Isoleucine","Methionine","Cysteine","Phenylalanine","Tyrosine","Tryptophan","Histidine","Lysine","Arginine","Glutamine","Asparagine","Glutamic_Acid","Aspartic_Acid","Serine","Threonine","iMethionine")
threecodes = c("Gly","Pro","Ala","Val","Leu","Ile","Met","Cys","Phe","Tyr","Trp","His","Lys","Arg","Gln","Asn","Glu","Asp","Ser","Thr","iMet")
onecodes = c("G","P","A","V","L","I","M","C","F","Y","W","H","K","R","Q","N","E","D","S","T","M")

#ends$trnaname = rownames(ends)
#print(head(ends))
tailmelt = melt(ends,id = c("end", "row.names"))
colnames(tailmelt) = c("tails","trna","sample","counts")
#print(head(tailmelt))

#tailmelt = tailmelt[tailmelt$sample %in% c("Liv_1",	"Liv_2"	,"Liv_AlkB_1","Liv_AlkB_2","Liv_PNK_AlkB_1","Liv_PNK_AlkB_2","Liv_PNK_AlkB_3"),]


tailmelt[tailmelt$tails == "CCA","counts"] <- tailmelt[tailmelt$tails == "CCA","counts"] + ccapseudocount
totalends <- aggregate(tailmelt$counts, by=list(trna = tailmelt$trna, sample = tailmelt$sample), FUN=sum)
colnames(totalends) = c("trna","sample","total")
#print(head(totalends))

#tailmelt$total = totalends[tailmelt$trna %in% totalends$trna & tailmelt$sample %in% totalends$sample,"x"]
#tailmelt$test = totalends[totalends$trna %in% tailmelt$trna & totalends$sample %in% tailmelt$sample,"trna"]
#match(tailmelt$trna,totalends$trna) & match(tailmelt$sample,totalends$sample)

tailmelt = merge(tailmelt, totalends, by = c("trna","sample"), all.x = TRUE)

tailmelt$percentcount = tailmelt$counts / tailmelt$total
#print(head(tailmelt))
tailmelt$amino = trnatable[match(tailmelt$trna,trnatable[,1]),3]
#print("**")
#print(head(tailmelt[tailmelt$percentcount > 1,]))
#print(max(tailmelt$percentcount))


tailmeltagg <- aggregate(tailmelt$percentcount, by=list(feature = tailmelt$trna,tails = tailmelt$tails, sample = Sampletable[match(tailmelt$sample,Sampletable[,1]),2] ), FUN=mean)
colnames(tailmeltagg) = c("feature","tails","sample","percentage")
tailmeltagg$amino = trnatable[match(tailmeltagg$feature,trnatable[,1]),3]
#tailmeltagg <- aggregate(tailmelt$percentcount, by=list(amino = tailmelt$amino,tails = tailmelt$tails, sample = Sampletable[match(tailmelt$sample,Sampletable[,1]),2] ), FUN=mean)
#colnames(tailmeltagg) = c("amino","tails","sample","percentage")

#print(head(tailmeltagg))


#aes(x = variable, y = value,fill = seq, stat="identity")
#deletemeltSampleagg <- aggregate(deletionmelt$percentdelete, by=list(Feature = deletionmelt$Feature, Sample = Sampletable[match(deletionmelt$Sample,Sampletable[,1]),2]), FUN=mean)
#colnames(deletemeltSampleagg) <- c("Feature","Sample","percentdeletions")
tailmeltagg$amino = trnatable[match(tailmeltagg$feature,trnatable[,1]),3]

tailmeltagg$tails = factor(tailmeltagg$tails, levels = c("Trimmed","C","CC","CCA"))
#tailmeltagg = tailmeltagg[tailmeltagg$amino == "Ala" & tailmeltagg$tails == "CCA" ,]

#print(tailmeltagg[tailmeltagg$sample == "Brain_plusAlkB_plusPNK2x",])


posname = paste(directory,"/",runname,"-trnatails",outputformat, sep = "")  #+ facet_wrap(~amino, scale="free")  #geom_jitter( size = dotsize,width = 0.25) 
ggplot(data = tailmeltagg, aes(x=tails, y=percentage,colour = sample,fill=sample)) + geom_boxplot()+theme_bw() +  scale_y_continuous(labels = percent_format(),limits=c(0,1)) +  ggtitle(paste("tRNA tails", sep = ""))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ xlab("Tail")+ylab("Percent")
#print(head(tailmeltagg))
#print(tailmeltagg)



#ggplot(tailmeltagg,aes(x = sample, y = percentage,fill = tails, stat="identity")) + theme_bw() + theme(panel.border = element_rect(linetype = "blank"), panel.grid = element_line(linetype = "blank")) + facet_wrap(~amino, scale="free")+
#	geom_bar(position = "fill",stat="identity") +
#    geom_bar(position = "fill",stat="identity",color="black",show.legend=FALSE) + 
#    scale_y_continuous(labels = percent_format()) +
#    theme(axis.text.x = element_text(size=5))+
#    xlab("Sample") +
#    ylab("Percentage of Tail Reads") + 
#    labs(fill="Tail\nType")+
#    #scale_fill_ucscgb()+
#    #scale_fill_brewer(palette = "Dark2")+
#    #scale_fill_manual(values = getPalette(colourCount))+
#    theme(axis.title.x = element_text(face="bold", size=15), axis.text.x = element_text(face="bold", size=9,angle = 90, vjust = .5)) #+scale_colour_gradient() #+ scale_fill_brewer( palette="RdPu")

    
ggsave(filename=posname,width=7, height=7)


