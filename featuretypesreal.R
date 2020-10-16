

library(ggplot2)
library(reshape2)
library(scales)

args <- commandArgs(trailingOnly = TRUE)

#args <- c("hg19-genomictrailertable.txt","hg19-genomicbarplot.png")
#args =c("repfragtypenormcounts.txt","barplot.png")


#Rscript trailerbarplot.R hg19-trailertable.txt hg19-barplot.png
counts <- read.table(args[1],check.names=FALSE)

selectcounts = counts


temp = cbind(selectcounts, seq = factor(rownames(selectcounts),rev(rownames(selectcounts))))

#levels(temp$seq) <- rev(rownames(selectcounts))

countsmelt = melt(temp, id.vars = c('seq'))


countsmelt = within(countsmelt, seq <- factor(seq, 
    rownames(selectcounts)))

#head(countsmelt)
sampletotals = aggregate(countsmelt$value, list(countsmelt$variable), sum)

#countsmelt


#sampletotals$x[countsmelt$variable]
#countsmelt = countsmelt[countsmelt$value > 100,]
countsmelt = countsmelt[countsmelt$value > sampletotals$x[countsmelt$variable] / 100,] #filters out those types below a certain level
#countsmelt = countsmelt
#countsmelt
ggplot(countsmelt,aes(x = variable, y = value,fill = seq, stat="identity")) + theme_bw() + theme(panel.border = element_rect(linetype = "blank"), panel.grid = element_line(linetype = "blank")) + 
	geom_bar(stat="identity") + 
    geom_bar(stat="identity",color="black",show.legend=FALSE) + 
    theme(axis.text.x = element_text(size=5))+
    xlab("Sample") +
    ylab("Total Reads") + 
    labs(fill="Read\nType") + 
    theme(axis.title.x = element_text(face="bold", size=15), axis.text.x = element_text(face="bold", size=9,angle = 90, vjust = .5))

ggsave(filename=args[2])
    