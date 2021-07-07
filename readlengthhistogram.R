

library(ggplot2)
library(reshape2)
library(scales)

args <- commandArgs(trailingOnly = TRUE)

#args <- c("hg19-genomictrailertable.txt","hg19-genomicbarplot.png")
#args =c("repfragtypenormcounts.txt","barplot.png")


#Rscript trailerbarplot.R hg19-trailertable.txt hg19-barplot.png
readlengths <- read.table(args[1], header = TRUE)
sampledata <- read.table(args[2], header = TRUE)







#temp = cbind(selectcounts, seq = factor(rownames(selectcounts),rev(rownames(selectcounts)), ordered = TRUE))

#levels(temp$seq) <- rev(rownames(selectcounts))

#countsmelt = melt(temp, id.vars = c('seq'))


if(FALSE){
ggplot(readlengths,aes(x = variable, y = value,fill = seq, stat="identity")) +
	geom_bar(position = "fill",stat="identity") + 
    geom_bar(position = "fill",stat="identity",color="black",show_guide=FALSE) + 
    scale_y_continuous(labels = percent_format()) +
    xlab("Sample") +
    ylab("Percentage of Reads") + 
    theme(axis.title.x = element_text(face="bold", size=15), axis.text.x = element_text(face="bold", size=9,angle = 90, vjust = .5))
    }
    

	


binw <- 2
maxlen <- min(c(200, max(readlengths[,1])))
readlengths$bin <- floor(readlengths$Length / binw)*binw

#head(readlengths)
binnedtrnareadlengths <- aggregate(trnas ~ bin * Sample, readlengths, sum)

#binnedreadlengths <- aggregate(other ~ bin, readlengths, sum)
binnedreadlengths <- aggregate(other ~ bin * Sample ,readlengths, sum)
#head(binnedreadlengths)
temp = merge(binnedtrnareadlengths,binnedreadlengths,by=c("bin", "Sample"))
#head(temp)
#temp
#temp = binnedtrnareadlengths
binnedlengths = melt(temp, id.vars = c('bin', 'Sample'))
#binnedlengths = temp
#head(binnedlengths)
binnedlengths[is.na(binnedlengths)] <- 0
#unique(binnedlengths$variable)
binnedlengths <- binnedlengths[order(binnedlengths$Sample,binnedlengths$bin),]

#binnedlengths <- binnedlengths[binnedlengths$variable =='trnas',]

#binnedlengths
#binnedlengths$variable <- factor(binnedlengths$variable, levels = c("trnas", "other"))
binnedlengths$variable <- factor(binnedlengths$variable, levels = c("other", "trnas"))
binnedlengths$Sample <- factor(binnedlengths$Sample, levels = sampledata[,1])


ggplot(data=binnedlengths, aes(x = binnedlengths$bin, y = binnedlengths$value, fill = binnedlengths$variable))+geom_bar(stat="identity")+facet_wrap( ~ Sample, scales="free",ncol = 3)+ xlim(0,maxlen) + ylab("Count") + xlab("Read Length") + scale_fill_discrete(name="Gene Type")  + scale_x_continuous(labels = comma) + theme_bw()
#length(unique(binnedlengths$Sample))
ggsave(filename=args[3],limitsize=FALSE,width = 8, height = .5 * length(unique(binnedlengths$Sample)))#, width = 3 * length(unique(binnedlengths$Sample))) 
    