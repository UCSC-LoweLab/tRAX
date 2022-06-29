
library(ggplot2)
library(gridExtra)
library(scales)
library(plyr)
library(RColorBrewer)
library(ggrepel)
library(dplyr)
library(reshape2)




#positions = c( "9", "14", "58","26","37","20","32","19")
distthresh = .05
args <- commandArgs(trailingOnly = TRUE) 



minscore = .02
minreads = 50
#samplecompare file
runname = args[1]
data = read.table(args[2],header = TRUE,row.names = NULL, stringsAsFactors=FALSE)
covdata = read.table(args[3],header = TRUE,row.names = NULL, stringsAsFactors=FALSE)
samplefilename = args[4]

sampletable <- read.table(samplefilename)

#tRNA-Tyr-GTA-2_27pos_G_clust


trimclusts = data[(data$firsamples > 1 | data$firsamples == max(data$firsamples)) & (data$secsamples > 1 | data$secsamples == max(data$secsamples)),]
#trimclusts = data

trimclusts$minpercent = trimclusts$firpercent - trimclusts$secpercent

trimclusts$pair = paste(trimclusts$firname, trimclusts$secname, sep = ":")


#+ scale_y_continuous(trans=reverselog_trans(10))
#unique(trimclusts$firpos)
#print(unique(trimclusts$pair))
allpos = unique(trimclusts$firpos)
for(currposition in allpos) {
#print(currposition)
#print(head(trimclusts))
positionclust = trimclusts[trimclusts$firpos == currposition,]
#print(head(positionclust))


for(currpair in unique(trimclusts$pair)){

pairclust = positionclust[positionclust$pair == currpair,]



#print(head(pairclust))
#print(max(pairclust$bdist,  na.rm = TRUE))
if(max(pairclust$bdist,  na.rm = TRUE) < distthresh){
next
}


#dplyr that should get the rows with max bdist for each trna
#pairclust = pairclust %>% group_by(firfeat) %>% slice(which.max(bdist))
#pairclust = as.data.frame(pairclust)

#print(head(pairclust))

bdistcutoff = -sort(-pairclust$bdist)[10]

displayfeats = ifelse( pairclust$bdist > bdistcutoff, paste(pairclust$firfeat, sep = ""), "")


#print(nrow(pairclust))
#print(nrow(trimclusts))

#maxdiffs = aggregate(positionclust$bdist, )


firsample = unique(pairclust$firname)[1]
secsample = unique(pairclust$secname)[1]
print("**")
ymax = max(pairclust$bdist,  na.rm = TRUE)
ycord = 0 #-(ymax * .01)
#print(ycord)
diffvolc =  ggplot(pairclust,aes(x = minpercent, y = bdist)) + theme_bw() + theme(panel.border = element_rect(linetype = "blank"), panel.grid = element_line(linetype = "blank")) + 
	        geom_point() +
            geom_text_repel(label = displayfeats,min.segment.length = unit(0, 'lines'), segment.color="red") + 
            theme(axis.text.x = element_text(size=5))+
            xlab("Difference in % identity") +
            ylab("Maximum Distance between Base Compositions") + 
            scale_x_continuous(labels=percent, limits = c(-1, 1)) +
            ggtitle(paste("position", currposition,currpair, sep = " "))+         
            coord_cartesian(clip = 'off') + theme(plot.margin = unit(c(1,3,1,1), "lines")) +
            labs(caption = c(firsample, secsample )) +  theme(plot.caption = element_text(size = 16,hjust=c(1, 0))) +
            theme(axis.title.x = element_text(face="bold", size=15), axis.text.x = element_text(face="bold", size=9,angle = 90, vjust = .5))  #+scale_colour_gradient() #+ scale_fill_brewer( palette="RdPu") +
            #annotation_custom(grid::textGrob(firsample, gp=grid::gpar(fontsize=13,fontface="bold")),xmin=.75, xmax=.75, ymin=ycord) + #theme(plot.margin = unit(c(1,3,1,1), "lines"))
            #annotation_custom(grid::textGrob(secsample, gp=grid::gpar(fontsize=13,fontface="bold")),xmin=-.75, xmax=-.75,  ymin=ycord)
            



ggsave(diffvolc, filename=paste(runname, currposition,currpair,"position.pdf", sep = "_"))

alltrnas = unique(pairclust$firfeat[!is.na(pairclust$firfeat) & pairclust$bdist > minscore])
print(head(pairclust[is.na(pairclust$firfeat),]))

for(currtrna in  alltrnas) {


# }
pairclust = pairclust[!is.na(pairclust$fircounts),]
refbase = unique(pairclust$firrefbase)[[1]]
clustname = paste(runname, currtrna, currposition,refbase,currpair, sep = "_")


firset = strsplit(pairclust$fircounts, ",")
#print(head(pairclust$fircounts))

#write.table(pairclust, "pairclust.txt")
#pairclust = head(pairclust)

firset = data.frame(matrix(unlist(strsplit(pairclust$fircounts, ",")), nrow=length(pairclust$fircounts), byrow=T))
#print(head(firset))
colnames(firset) = c("Adenine","Cytosine","Guanine","Thymine","Deletions")
firset$samplename = pairclust$firname

#print(head(firset))

secset = data.frame(matrix(unlist(strsplit(pairclust$seccounts, ",")), nrow=length(pairclust$seccounts), byrow=T))
colnames(secset) = c("Adenine","Cytosine","Guanine","Thymine","Deletions")
secset$samplename = pairclust$secname
#print(head(pairclust))


allset = rbind(firset, secset)

allset = allset[!(duplicated(allset$samplename)),] 

rownames(allset) <- allset$samplename
allset = allset[,c("Adenine","Cytosine","Guanine","Thymine", "Deletions")]
allset = apply(allset,c(1,2),as.numeric)
allset = t(allset)

allsamples = unique(sampletable[,2])

#levels(temp$seq) <- rev(rownames(selectcounts))
countsmelt = melt(allset)
print(head(countsmelt))

ggplot(countsmelt,aes(x = Var2, y = value,fill = Var1, stat="identity")) + theme_bw() + theme(panel.border = element_rect(linetype = "blank"), panel.grid = element_line(linetype = "blank")) + 
	geom_bar(position = "stack",stat="identity") +
    geom_bar(position = "stack",stat="identity",color="black",show.legend=FALSE) + 
    theme(axis.text.x = element_text(size=5))+
    xlab("Sample") +
    ylab("Count Reads at position") + 
    labs(fill="Base")+
    ggtitle(clustname)+                                                              
    #scale_fill_ucscgb()+
    #scale_fill_brewer(palette = "Dark2")+
    #scale_fill_manual(values = getPalette(colourCount))+
    scale_x_discrete(limits =  allsamples)+
    theme(axis.title.x = element_text(face="bold", size=15), axis.text.x = element_text(face="bold", size=9,angle = 90, vjust = .5)) #+scale_colour_gradient() #+ scale_fill_brewer( palette="RdPu")

ggsave(filename=paste(clustname, "baseaggcounts.pdf", sep = "_"))

allreps = unique(sampletable[,1])
#print(head(covdata))
currdata = covdata[covdata$Feature == currtrna & covdata$position == currposition ,c("Sample","adenines","cytosines","guanines","thymines", "deletions")]
allmelt = melt(currdata,  id.vars = "Sample")
#print(head(allmelt))



ggplot(allmelt,aes(x = Sample, y = value,fill = variable, stat="identity")) + theme_bw() + theme(panel.border = element_rect(linetype = "blank"), panel.grid = element_line(linetype = "blank")) + 
	geom_bar(position = "fill",stat="identity") +
    geom_bar(position = "fill",stat="identity",color="black",show.legend=FALSE) + 
    theme(axis.text.x = element_text(size=5))+
    xlab("Sample") +
    ylab("Count Reads at position") + 
    labs(fill="Base")+
    ggtitle(clustname)+                                                              
    #scale_fill_ucscgb()+
    #scale_fill_brewer(palette = "Dark2")+
    #scale_fill_manual(values = getPalette(colourCount))+
    scale_x_discrete(limits =  allreps)+
    theme(axis.title.x = element_text(face="bold", size=15), axis.text.x = element_text(face="bold", size=9,angle = 90, vjust = .5)) #+scale_colour_gradient() #+ scale_fill_brewer( palette="RdPu")

ggsave(filename=paste(clustname, "basecountspercent2.pdf", sep = "_"))


ggplot(allmelt,aes(x = Sample, y = value,fill = variable, stat="identity")) + theme_bw() + theme(panel.border = element_rect(linetype = "blank"), panel.grid = element_line(linetype = "blank")) + 
	geom_bar(position = "stack",stat="identity") +
    geom_bar(position = "stack",stat="identity",color="black",show.legend=FALSE) + 
    theme(axis.text.x = element_text(size=5))+
    xlab("Sample") +
    ylab("Count Reads at position") + 
    labs(fill="Base")+
    ggtitle(clustname)+                                                              
    #scale_fill_ucscgb()+
    #scale_fill_brewer(palette = "Dark2")+
    #scale_fill_manual(values = getPalette(colourCount))+
    scale_x_discrete(limits =  allreps)+
    theme(axis.title.x = element_text(face="bold", size=15), axis.text.x = element_text(face="bold", size=9,angle = 90, vjust = .5)) #+scale_colour_gradient() #+ scale_fill_brewer( palette="RdPu")

ggsave(filename=paste(clustname, "basecounts2.pdf", sep = "_"))

#trimcounts = allmelt[maxdiff[countsmelt$Var2,"x"] > 50,]




}
}
}