
library(ggplot2)
library(reshape2)
library(scales)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)

outtype = ""
if(length(args) > 2){
outtype = args[3]
}


#args <- c("hg19-genomictrailertable.txt","hg19-genomicbarplot.png")
#args =c("repfragtypenormcounts.txt","barplot.png")


#Rscript trailerbarplot.R hg19-trailertable.txt hg19-barplot.png
counts <- read.table(args[1],check.names=FALSE)

selectcounts = counts


temp = cbind(selectcounts, seq = factor(rownames(selectcounts),rev(rownames(selectcounts)), ordered = TRUE))

#levels(temp$seq) <- rev(rownames(selectcounts))

countsmelt = melt(temp, id.vars = c('seq'))


countsmelt = within(countsmelt, seq <- factor(seq, 
    rownames(selectcounts), 
    ordered = FALSE))

#head(countsmelt)
sampletotals = aggregate(countsmelt$value, list(countsmelt$variable), sum)

#countsmelt


#sampletotals$x[countsmelt$variable]
#countsmelt = countsmelt[countsmelt$value > 100,]
countsmelt = countsmelt[countsmelt$value > 0,]

if(outtype != "all"){
countsmelt = countsmelt[countsmelt$value > sampletotals$x[countsmelt$variable] / 100,]
}
#countsmelt = countsmelt
#head(countsmelt)

#unique(countsmelt$seq)
colourCount = length(unique(countsmelt$seq))+1
#getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

#unique(head(countsmelt))

#   other #f2dab2
#   tRNA #6666cc
#   pretRNA #a1ade5
#   miRNA #b4008d
#   snRNA #74531d
#   Mt_tRNA #87e275
#   Mt_rRNA #00754a
#   rRNA    #a5e5d9
#   snoRNA  #ff9e18
#   misc_RNA    #b2b2b2


#print(unique(countsmelt$seq))




ggplot(countsmelt,aes(x = variable, y = value,fill = seq, stat="identity")) + theme_bw() + theme(panel.border = element_rect(linetype = "blank"), panel.grid = element_line(linetype = "blank")) + 
	geom_bar(position = "fill",stat="identity") +
    geom_bar(position = "fill",stat="identity",color="black",show.legend=FALSE) + 
    scale_y_continuous(labels = percent_format()) +
    theme(axis.text.x = element_text(size=5))+
    xlab("Sample") +
    ylab("Percentage of Total Reads") + 
    labs(fill="Read\nType")+
    #scale_fill_ucscgb()+
    #scale_fill_brewer(palette = "Dark2")+
    #scale_fill_manual(values = getPalette(colourCount))+
    #scale_fill_manual(values = typepal)+
    theme(axis.title.x = element_text(face="bold", size=15), axis.text.x = element_text(face="bold", size=9,angle = 90, vjust = .5)) #+scale_colour_gradient() #+ scale_fill_brewer( palette="RdPu")

ggsave(filename=args[2])
    