
library(ggplot2)
library(reshape2)
library(scales)

args <- commandArgs(trailingOnly = TRUE)


data <- read.table(args[1])
samplefile = args[2]
sampledata = read.table(samplefile)



#data <- data[,grepl("tRNA",colnames(data),fixed = TRUE)]
#data <- data[grepl("tRNA",rownames(data),fixed = TRUE),]
data <- data[rowSums(data) > 20,]
#length(colnames(data))
datapca <- prcomp(t(data),center = TRUE,scale = TRUE) 


percentlabels <- round(datapca$sdev / sum(datapca$sdev) * 100, 2)
percentlabels <- paste( colnames(datapca$x), "(", paste( as.character(percentlabels), "%", ")", sep="") )
#


scores = as.data.frame(datapca$x)

#print(head(scores))
# plot of observations

#print(sampledata[match(rownames(scores),sampledata[,1]),2])
#print(rownames(scores))
#print(sampledata)
#aes(colour = factor(cyl))

samplename = sampledata[match(rownames(scores),sampledata[,1]),2]

ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores), color = samplename)) + theme_bw()+labs(color="Sample Name")+ geom_point()+geom_hline(yintercept = 0, colour = "gray65") +geom_vline(xintercept = 0, colour = "gray65") + geom_text(alpha = 0.8, size = 2,vjust="inward",hjust="inward") + ggtitle("Principle Component Analysis")    + xlab(percentlabels[1]) +    ylab(percentlabels[2]) 
ggsave(filename=args[3])

