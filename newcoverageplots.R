library(ggplot2)
library(reshape2)
library(scales)
library(plyr)
library(gridExtra)
library(getopt)

library(ggplot2)
library(reshape2)
library(scales)
library(getopt)

allaminos = c('Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly','His','Ile','Ile2','Leu','Lys','Met','iMet','fMet','Phe','Pro','Ser','Thr','Trp','Tyr','Val','SeC','Sup','Undet') 

fiveprimename = "5p end"
dstartname = "D-loop start"
dendname = "D-loop end"
acstartname = "AC-loop start"
acendname = "AC-loop end"
tstartname = "T-loop start"
tendname = "T-loop end"
ccatailname= "3p CCA tail"


expand.delimited <- function(x, col1=1, col2=2, sep=",") {
  rnum <- 1
  expand_row <- function(y) { 
    factr <- y[col1]
    strng <- toString(y[col2])
    expand <- strsplit(strng, sep)[[1]]
    num <- length(expand)
    factor <- rep(factr,num)
    return(as.data.frame(cbind(factor,expand),
          row.names=seq(rnum:(rnum+num)-1)))
    rnum <- (rnum+num)-1
  }
  expanded <- apply(x,1,expand_row)
  df <- do.call("rbind", expanded)
  names(df) <- c(names(x)[col1],names(x)[col2])
  return(df)
}

myBreaks <- function(x){
    breaks <- c(min(0),floor(max(x)))
    names(breaks) <- attr(breaks,"labels")
    breaks
}
percentbreaks <- function(x){
    breaks <- c(0,1)
    names(breaks) <- attr(breaks,"labels")
    names(breaks) <- c("0%", "100%")
    breaks
}
makebasiccovplot <-  function(covdata, filename){

#show_col(hue_pal()(4))

#uniquecols <- c("Transcript Specific" = "#C77CFF", "Isodecoder Specific" = "#00BEC4","Isotype Specific" = "#7CAE00", "Not Amino Specific" = "#F8766D")

#uniquecols <- c("#C77CFF","#00BEC4", "#7CAE00","#F8766D")



allfeatures = unique(covdata$Feature)

#firstsample = unique(covdata$Feature)[1]
#scalingpoint <- data.frame(Feature = allfeatures, Sample =  rep(length(allfeatures),firstsample),variable = rep("1",length(allfeatures)), value = rep(50,length(allfeatures)))
# geom_point(data = scalingpoint, aes(x = variable, y = value), alpha = 0)

#uniquecoverage	multitrnacoverage	multianticodoncoverage	multiaminocoverage
coverage   <- ggplot(covdata,aes(x=variable,y=value), size = 2)  +theme_bw()+facet_grid(Feature ~ Sample, scales="free") + expand_limits(y = 50)+ geom_bar(stat="identity") +  geom_vline(aes(xintercept = dist, col = Modification),data = aminomodomicstable,show.legend=TRUE,size=.2,linetype = "longdash")+theme(axis.text.y=element_text(colour="black",size=6),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8)) + ylab("Normalized Read Count") + xlab("tRNA position") + scale_y_continuous(breaks=myBreaks) +scale_x_discrete(breaks=c("1","13","22","31","39","53","61","73"), labels=c(fiveprimename,dstartname,dendname,acstartname,acendname,tstartname,tendname,ccatailname)) + labs(fill="Mappability", vline="RNA\nModification") #+ scale_color_hue("mod")+labs(fill="RNA\nModification")#+scale_fill_manual(name="RNA\nmodifications")#+scale_colour_manual(data = aminomodomicstable, name="RNA\nModification") #+scale_x_discrete(breaks=c("1","37","73"), labels=c(fiveprimename,"anticodon", ccatailname))

coverage <- configurecov(coverage)

ggsave(filename=filename, coverage,height=scalefactor*(2 + 1.5*length(unique(covdata$Feature))),width=scalefactor*(2+5*length(unique(covdata$Sample))), limitsize=FALSE, dpi = 600)

}

makecovplot <-  function(covdata, filename){

#show_col(hue_pal()(4))

uniquecols <- c("Transcript Specific" = "#C77CFF", "Isodecoder Specific" = "#00BEC4","Isotype Specific" = "#7CAE00", "Not Amino Specific" = "#F8766D")



allsamples = unique(covdata$Sample)
allfeatures = unique(covdata$Feature)

#maplist = 

#covdata$maptype = factor(covdata$maptype,levels = c("Transcript Specific", "Isodecoder Specific","Isotype Specific", "Not Amino Specific"))
covdata$maptype = factor(covdata$maptype,levels = c("Not Amino Specific","Isotype Specific","Isodecoder Specific","Transcript Specific"))



#covdata = covdata[order(as.numeric(covdata$maptype)),]

#firstsample = unique(covdata$Feature)[1]
#scalingpoint <- data.frame(Feature = allfeatures, Sample =  rep(firstsample,length(allfeatures)),variable = rep("1",length(allfeatures)), value = rep(50,length(allfeatures)), maptype = rep("Transcript Specific",length(allfeatures)))
#scaling = geom_point(data = scalingpoint, aes(x = variable, y = value), alpha = 0)

#uniquecoverage	multitrnacoverage	multianticodoncoverage	multiaminocoverage
coverage   <- ggplot(covdata,aes(x=variable,y=value, fill = maptype, order=-as.numeric(maptype)), size = 2) +scale_fill_manual(values = uniquecols)+ expand_limits(y = 50)+ theme_bw()+facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity") +  geom_vline(aes(xintercept = dist, col = Modification),data = aminomodomicstable,show.legend=TRUE,size=.2,linetype = "longdash")+theme(axis.text.y=element_text(colour="black",size=6),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8)) + ylab("Normalized Read Count") + xlab("tRNA position") + scale_y_continuous(breaks=myBreaks) +scale_x_discrete(breaks=c("1","13","22","31","39","53","61","73"), labels=c(fiveprimename,dstartname,dendname,acstartname,acendname,tstartname,tendname,ccatailname)) + labs(fill="Mappability", vline="RNA\nModification") #+ scale_color_hue("mod")+labs(fill="RNA\nModification")#+scale_fill_manual(name="RNA\nmodifications")#+scale_colour_manual(data = aminomodomicstable, name="RNA\nModification") 

coverage <- configurecov(coverage)

ggsave(filename=filename, coverage,height=scalefactor*(2 + 1.5*length(unique(covdata$Feature))),width=scalefactor*(2+5*length(unique(covdata$Sample))), limitsize=FALSE, dpi = 600)

}




makeendplot <-  function(covdata, filename){

#show_col(hue_pal()(4))

uniquecols <- c("Read Starts" = "#C77CFF", "Read ends" = "#00BEC4","Coverage" = "#7CAE00")


allsamples = unique(covdata$Sample)
allfeatures = unique(covdata$Feature)
#print("**")

covdata$endtype = factor(covdata$endtype,levels = c("Coverage","Read Starts", "Read ends"))


#print("**||")
 #+scale_fill_manual(values = uniquecols)
coverage   <- ggplot(covdata,aes(x=position,y=value, fill = endtype, order=-as.numeric(endtype)), size = 2) + expand_limits(y = 50)+ theme_bw()+facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity") +  geom_vline(aes(xintercept = dist, col = Modification),data = aminomodomicstable,show.legend=TRUE,size=.2,linetype = "longdash")+theme(axis.text.y=element_text(colour="black",size=6),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8)) + ylab("Normalized Read Count") + xlab("tRNA position") + scale_y_continuous(breaks=myBreaks) +scale_x_discrete(breaks=c("1","13","22","31","39","53","61","73"), labels=c(fiveprimename,dstartname,dendname,acstartname,acendname,tstartname,tendname,ccatailname)) + labs(fill="Mappability", vline="RNA\nModification") #+ scale_color_hue("mod")+labs(fill="RNA\nModification")#+scale_fill_manual(name="RNA\nmodifications")#+scale_colour_manual(data = aminomodomicstable, name="RNA\nModification") #+scale_x_discrete(breaks=c("1","37","73"), labels=c(fiveprimename,"anticodon", ccatailname))

coverage <- configurecov(coverage)

ggsave(filename=filename, coverage,height=scalefactor*(2 + 1.5*length(unique(covdata$Feature))),width=scalefactor*(2+5*length(unique(covdata$Sample))), limitsize=FALSE, dpi = 600)

}

makepercentcovplot <-  function(covdata, filename){


allfeatures = unique(covdata$Feature)

coverage   <- ggplot(covdata,aes(x=variable,y=value), size = 2)   +theme_bw()+facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity") +  geom_vline(aes(xintercept = dist, col = Modification),data = aminomodomicstable,show.legend=TRUE,size=.2,linetype = "longdash")+theme(axis.text.y=element_text(colour="black",size=6),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8)) + ylab("Percentage of Read Coverage") + xlab("tRNA position") + scale_y_continuous(breaks=percentbreaks,limits = c(0, 1.01)) +scale_x_discrete(breaks=c("1","13","22","31","39","53","61","73"), labels=c(fiveprimename,dstartname,dendname,acstartname,acendname,tstartname,tendname,ccatailname)) + labs(fill="Mappability", vline="RNA\nModification") #+ scale_color_hue("mod")+labs(fill="RNA\nModification")#+scale_fill_manual(name="RNA\nmodifications")#+scale_colour_manual(data = aminomodomicstable, name="RNA\nModification") #+scale_x_discrete(breaks=c("X1","X37","X73"), labels=c(fiveprimename,"anticodon", ccatailname))

coverage <- configurecov(coverage)

ggsave(filename=filename, coverage,height=scalefactor*(2 + 1.5*length(unique(covdata$Feature))),width=scalefactor*(2+5*length(unique(covdata$Sample))), limitsize=FALSE, dpi = 600)

}




makecombplot <-  function(covdata, filename, indepscales = FALSE){
if(indepscales){
plotscales = "free"
}else{
plotscales = "free"
}
smallcovsummary <- ggplot(covdata,aes(x=variable,y=value, fill = amino, order = as.numeric(sortacceptor)),width = 2, size = 2) + 
    geom_bar(stat="identity")+ 
    
    xlab("Position")+ 
    ylab("Normalized Read Count") +   
    scale_y_continuous(breaks=myBreaks) +scale_x_discrete(breaks=c("X1","X37","X73"), labels=c(fiveprimename,"anticodon",ccatailname),expand = c(0.05, .01)) +
    scale_fill_discrete(drop=FALSE, name="Acceptor\ntype", breaks = levels(sortacceptor))  +
    facet_grid( ~ Sample, scales="free", space= "free_y") 

smallcovsummary <- configurecombine(smallcovsummary)

combinescale = 3
ggsave(filename=filename, smallcovsummary, width =combinescale*(.5+.75*length(unique(covdata$Sample))), height = combinescale*1,limitsize = FALSE)
}
makelocuscombplot <-  function(covdata, filename){

smallcovsummary <- ggplot(covdata,aes(x=variable,y=value, fill = sortacceptor, order = as.numeric(sortacceptor)),width = 2, size = 2	) + facet_grid( ~ Sample, scales="free") +xlab("Position")+ geom_bar(stat="identity")+ ylab("Normalized Read Count") +   scale_y_continuous(breaks=myBreaks) +scale_x_discrete(breaks=c("X1","X37","X73"), labels=c(fiveprimename,"anticodon",ccatailname),expand = c(0.05, .01)) +scale_fill_discrete(drop=FALSE, name="Acceptor\ntype", breaks = levels(sortacceptor))

smallcovsummary <- configurecombine(smallcovsummary)

combinescale = 3
ggsave(filename=filename, smallcovsummary, width =combinescale*(.5+.75*length(unique(covdata$Sample))), height = combinescale*1,limitsize = FALSE)

}																															#																																	#, expand = c(0.05, 0)




percentbreaks <- function(x){
    breaks <- c(0,1)
    names(breaks) <- attr(breaks,"labels")
    breaks
}
configurecov <- function(covplot){
covplot<- covplot+ theme_bw()
covplot<- covplot+theme(text = element_text(size = 4))
covplot<- covplot+theme(rect = element_rect(size=.01))
              
covplot<- covplot+theme(legend.key.size =  unit(.3, "cm")) # Change key size in the legend
covplot<- covplot+theme(legend.key.height =  unit(.3, "cm")) # Change key size in the legend
              
covplot<- covplot+theme(legend.text = element_text(size=4)) # Change the size labels in the legend
covplot<- covplot+theme(legend.title = element_text(size=4))             
              
covplot<- covplot+theme(axis.ticks = element_line(size = .1)) 
#smallcovall< covplot+theme(axis.ticks.y=element_blank())
              
covplot<- covplot+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=6))

              
covplot<- covplot+theme(axis.line = element_line(size = .1))
covplot<- covplot+theme(panel.grid  = element_line(size = .1))
              
              
covplot<- covplot+theme(axis.line = element_line(size = .1))
covplot<- covplot+theme(panel.grid  = element_line(size = .15))
              
covplot<- covplot+theme(strip.switch.pad.grid = unit(.1, "cm"))

#covplot<- covplot+theme(panel.margin = unit(.25, "lines"))
covplot<- covplot+theme(panel.spacing = unit(.25, "lines"))


#This needs to be different from the combined version for some reason
covplot<- covplot+theme(axis.text.y=element_text(size=6))
            
covplot<- covplot+theme(strip.text.y = element_text(angle=0,size=10))
#as does this  
covplot<- covplot+theme(strip.text.x = element_text(angle=0,size=10))
covplot<- covplot+theme(axis.text=element_text(size=8))
covplot<- covplot+theme(axis.title=element_text(size=8))

covplot
}

configurecombine <- function(covplot){
covplot<- covplot+ theme_bw()
covplot<- covplot+theme(text = element_text(size = 4))
covplot<- covplot+theme(rect = element_rect(size=.01))
              
covplot<- covplot+theme(legend.key.size =  unit(.4, "cm")) # Change key size in the legend
covplot<- covplot+theme(legend.key.height =  unit(.4, "cm")) # Change key size in the legend
              
covplot<- covplot+theme(legend.text = element_text(size=6)) # Change the size labels in the legend
covplot<- covplot+theme(legend.title = element_text(size=8))             
              
covplot<- covplot+theme(axis.ticks = element_line(size = .1)) 
#smallcovall< covplot+theme(axis.ticks.y=element_blank())
              
covplot<- covplot+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=6))


covplot<- covplot+theme(axis.line = element_line(size = .1))
covplot<- covplot+theme(panel.grid  = element_line(size = .15))
              
covplot<- covplot+theme(strip.switch.pad.grid = unit(.1, "cm"))
#covplot<- covplot+theme(panel.margin = unit(.15, "lines"))
#This needs to be different from the combined version for some reason
covplot<- covplot+theme(axis.text.y=element_text(size=6))
            
covplot<- covplot+theme(strip.text.y = element_text(angle=0,size=10))
#as does this  
covplot<- covplot+theme(strip.text.x = element_text(angle=0,size=10))
covplot<- covplot+theme(axis.text=element_text(size=8))
covplot<- covplot+theme(axis.title=element_text(size=8))

covplot
}

coverageprep <- function(coveragemeltagg, samples, trnatable) {

colnames(coveragemeltagg)[colnames(coveragemeltagg) == "x"]  <- "value"
coveragemeltagg$Sample <- factor(coveragemeltagg$Sample,levels = unique(sampletable[,2]), ordered = TRUE)


featnames = unique(as.character(coveragemeltagg$Feature))
#tails = as.numeric(unlist(lapply(strsplit(featnames, "-", fixed = TRUE), tail, 1)))
#anticodonname = as.character(unlist(lapply(strsplit(featnames, "-", fixed = TRUE), function(x) { return( x[length(x) - 1] ) })))
#aminoname = as.character(unlist(lapply(strsplit(featnames, "-", fixed = TRUE), function(x) { return( x[length(x) - 2] ) })))
#coveragemeltagg$Feature = factor(as.character(coveragemeltagg$Feature), levels = featnames)
#coveragemeltagg$amino = factor(as.character(coveragemeltagg$Feature), levels = featnames)


#colnames(trnatable) <- c("trnaname", "loci", "amino", "anticodon")


trnainfo = trnatable[,c("trnaname", "amino", "anticodon")]
coveragemeltagg <- merge(coveragemeltagg,trnainfo, by.x = "Feature", by.y = "trnaname", all.x=TRUE)



    
return (coveragemeltagg)
}

args <- commandArgs(trailingOnly = TRUE)

#source("traxlib.R")
#print(dirname(sys.frame(1)$ofile))
#print(here())
#print(file.path())

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dirname <- dirname(script.name)
#print(paste(script.dirname,"traxlib.R",sep="/"))

#source(paste(script.dirname,"traxlib.R",sep="/"))

#hcvm.i <- melt(hcvm.i, id.vars=c(grep("^readC", names(hcvm.i), value=TRUE, invert=TRUE)), variable.name="tRNA.basePosition", value.name="read.density")


spec <- matrix(c(
        'cov'     , 'c', 1, "character", "coverage file from getcoverage.py (required)",
        'locicov'     , 'l', 1, "character", "coverage file from getcoverage.py (required)",
        'trna'     , 't', 1, "character", "trna file (required)",
        'samples'    , 's', 1, "character", "sample file (required)",
        'runname'    , 'r', 1, "character", "name of TRAX sample",
        'directory'    , 'f', 1, "character", "directory to place amino acid files",
        'allcov'    , 'a', 1, "character", "output coverages for all tRNAs (optional)",
        'multicov'    , 'm', 1, "character", "output coverages for all tRNAs on multiple pages(optional)",
        
        'combinecov'    , 'o', 1, "character", "output coverages for tRNAs combined",
        'modomics'    , 'd', 2, "character", "optional modomics file",
        'help'   , 'h', 0, "logical",   "this help"
),ncol=5,byrow=T)



#trnapositions = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76")

#trnapositions = c('0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','17a','18','19','20','20a','20b','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14','e15','e16','e17','e18','e19','46','47','48','49','50','51','52','53','54','55','56','57','58','59','60','61','62','63','64','65','66','67','68','69','70','71','72','73','74','75','76')

trnapositions = c('-1','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','17a','18','19','20','20a','20b','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','46','47','48','49','50','51','52','53','54','55','56','57','58','59','60','61','62','63','64','65','66','67','68','69','70','71','72','73','74','75','76')

#locuspositions = c("head30","head29","head28","head27","head26","head25","head24","head23","head22","head21","head20","head19","head18","head17","head16","head15","head14","head13","head12","head11","head10","head9","head8","head7","head6","head5","head4","head3","head2","head1","0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","tail1","tail2","tail3","tail4","tail5","tail6","tail7","tail8","tail9","tail10","tail11","tail12","tail13","tail14","tail15","tail16","tail17","tail18","tail19","tail20","tail21","tail22","tail23","tail24","tail25","tail26","tail27","tail28","tail29","tail30")


#locuspositions = c("head30","head29","head28","head27","head26","head25","head24","head23","head22","head21","head20","head19","head18","head17","head16","head15","head14","head13","head12","head11","head10","head9","head8","head7","head6","head5","head4","head3","head2","head1",'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','17a','18','19','20','20a','20b','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14','e15','e16','e17','e18','e19','46','47','48','49','50','51','52','53','54','55','56','57','58','59','60','61','62','63','64','65','66','67','68','69','70','71','72','73','74','75','76',"tail1","tail2","tail3","tail4","tail5","tail6","tail7","tail8","tail9","tail10","tail11","tail12","tail13","tail14","tail15","tail16","tail17","tail18","tail19","tail20","tail21","tail22","tail23","tail24","tail25","tail26","tail27","tail28","tail29","tail30")

locuspositions = c("head30","head29","head28","head27","head26","head25","head24","head23","head22","head21","head20","head19","head18","head17","head16","head15","head14","head13","head12","head11","head10","head9","head8","head7","head6","head5","head4","head3","head2","head1",'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','17a','18','19','20','20a','20b','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','46','47','48','49','50','51','52','53','54','55','56','57','58','59','60','61','62','63','64','65','66','67','68','69','70','71','72','73','74','75','76',"tail1","tail2","tail3","tail4","tail5","tail6","tail7","tail8","tail9","tail10","tail11","tail12","tail13","tail14","tail15","tail16","tail17","tail18","tail19","tail20","tail21","tail22","tail23","tail24","tail25","tail26","tail27","tail28","tail29","tail30")


#locuspositions = c("-30","-29","-28","-27","-26","-25","-24","-23","-22","-21","-20","-19","-18","-17","-16","-15","-14","-13","-12","-11","-10","-9","-8","-7","-6","-5","-4","-3","-2","-1",'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','17a','18','19','20','20a','20b','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','46','47','48','49','50','51','52','53','54','55','56','57','58','59','60','61','62','63','64','65','66','67','68','69','70','71','72','73','74','75','76',"+1","+2","+3","+4","+5","+6","+7","+8","+9","+10","+11","+12","+13","+14","+15","+16","+17","+18","+19","+20","+21","+22","+23","+24","+25","+26","+27","+28","+29","+30")

opt = getopt(spec);

coverages <- read.table(opt$cov, header = TRUE,row.names = NULL, stringsAsFactors=FALSE)
locicoverages <- read.table(opt$locicov, header = TRUE,row.names = NULL, stringsAsFactors=FALSE)

trnatable <- read.table(opt$trna)
colnames(trnatable) <- c("trnaname", "loci", "amino", "anticodon")
sampletable <- read.table(opt$samples)
expname <- opt$name
modomicstable <- data.frame()
scalefactor = .5
modomicstable <- read.table(text ="",col.names = c("trna", "mod", "pos"),colClasses = c("character", "character", "character")) #(trna = character(), mod = character(), pos = character(),stringsAsFactors=FALSE)
if(!is.null(opt$modomics) && file.exists(opt$modomics) ){
    modomicstable <- read.table(opt$modomics, header = TRUE)
    #modomicstable$pos <- paste("X",modomicstable$pos, sep = "")
}

colnames(modomicstable)[colnames(modomicstable) == 'mod'] <- 'Modification'

modomicstable <- modomicstable[as.character(modomicstable$pos) %in% unique(as.character(trnapositions)),]
modomicstable$dist <- match(modomicstable$pos,  levels(trnapositions)) 
#modomicstable$Feature <- factor(modomicstable$trna,levels = levels(coveragemelt$Feature) )

stopmods = c("m1A","m2,2G","m1G","m1I","m3C")
modomicstable <- modomicstable[modomicstable$mod %in% stopmods,]



aminomodomicstable <- modomicstable



runname = opt$runname
uniquename = paste(runname,runname,sep = "/")


outputformat <- ".pdf"
outputindivformat <- ".pdf"





outputfile <- opt$allcov
combinedfile <- opt$combinecov
multipage <- opt$multicov
#colnames(coverages)

myBreaks <- function(x){
    breaks <- c(min(0),floor(max(x)))
    names(breaks) <- attr(breaks,"labels")
    breaks
}
percentbreaks <- function(x){
    breaks <- c(0,1)
    names(breaks) <- attr(breaks,"labels")
    breaks
}
configurecov <- function(covplot){
covplot<- covplot+ theme_bw()
covplot<- covplot+theme(text = element_text(size = 4))
covplot<- covplot+theme(rect = element_rect(size=.01))
              
covplot<- covplot+theme(legend.key.size =  unit(.3, "cm")) # Change key size in the legend
covplot<- covplot+theme(legend.key.height =  unit(.3, "cm")) # Change key size in the legend
              
covplot<- covplot+theme(legend.text = element_text(size=4)) # Change the size labels in the legend
covplot<- covplot+theme(legend.title = element_text(size=4))             
              
covplot<- covplot+theme(axis.ticks = element_line(size = .1)) 
#smallcovall< covplot+theme(axis.ticks.y=element_blank())
              
covplot<- covplot+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=6))

              
covplot<- covplot+theme(axis.line = element_line(size = .1))
covplot<- covplot+theme(panel.grid  = element_line(size = .1))
              
              
covplot<- covplot+theme(axis.line = element_line(size = .1))
covplot<- covplot+theme(panel.grid  = element_line(size = .15))
              
covplot<- covplot+theme(strip.switch.pad.grid = unit(.1, "cm"))
#covplot<- covplot+theme(panel.margin = unit(.25, "lines"))
covplot<- covplot+theme(panel.spacing = unit(.25, "lines"))

#This needs to be different from the combined version for some reason
covplot<- covplot+theme(axis.text.y=element_text(size=6))
            
covplot<- covplot+theme(strip.text.y = element_text(angle=0,size=10))
#as does this  
covplot<- covplot+theme(strip.text.x = element_text(angle=0,size=10))
covplot<- covplot+theme(axis.text=element_text(size=8))
covplot<- covplot+theme(axis.title=element_text(size=8))

covplot
}

configurecombine <- function(covplot){
covplot<- covplot+ theme_bw()
covplot<- covplot+theme(text = element_text(size = 4))
covplot<- covplot+theme(rect = element_rect(size=.01))
              
covplot<- covplot+theme(legend.key.size =  unit(.4, "cm")) # Change key size in the legend
covplot<- covplot+theme(legend.key.height =  unit(.4, "cm")) # Change key size in the legend
              
covplot<- covplot+theme(legend.text = element_text(size=6)) # Change the size labels in the legend
covplot<- covplot+theme(legend.title = element_text(size=8))             
              
covplot<- covplot+theme(axis.ticks = element_line(size = .1)) 
#smallcovall< covplot+theme(axis.ticks.y=element_blank())
              
covplot<- covplot+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=6))

              
covplot<- covplot+theme(axis.line = element_line(size = .1))
covplot<- covplot+theme(panel.grid  = element_line(size = .1))
              
              
covplot<- covplot+theme(axis.line = element_line(size = .1))
covplot<- covplot+theme(panel.grid  = element_line(size = .15))
              
covplot<- covplot+theme(strip.switch.pad.grid = unit(.1, "cm"))
#covplot<- covplot+theme(panel.margin = unit(.15, "lines"))
covplot<- covplot+theme(panel.spacing = unit(.15, "lines"))

#This needs to be different from the combined version for some reason
covplot<- covplot+theme(axis.text.y=element_text(size=6))
            
covplot<- covplot+theme(strip.text.y = element_text(angle=0,size=10))
#as does this  
covplot<- covplot+theme(strip.text.x = element_text(angle=0,size=10))
covplot<- covplot+theme(axis.text=element_text(size=8))
covplot<- covplot+theme(axis.title=element_text(size=8))

covplot
}
coverageall = coverages

coverageall = coverageall[coverageall$position %in% trnapositions,]
coverageall$position = factor(coverageall$position, levels=trnapositions)


coveragemeltagg <- aggregate(coverageall$coverage, by=list(Feature = coverageall$Feature, Sample = sampletable[match(coverageall$Sample,sampletable[,1]),2], variable = coverageall$position), FUN=mean)

coveragemeltagg = coveragemeltagg[coveragemeltagg$variable %in% trnapositions,]

coveragemeltagg$variable = factor(coveragemeltagg$variable, levels=trnapositions)


colnames(coveragemeltagg)[colnames(coveragemeltagg) == "x"]  <- "value"
coveragemeltagg$Sample <- factor(coveragemeltagg$Sample,levels = unique(sampletable[,2]), ordered = TRUE)


#acceptorType = trnatable[match(coveragemeltagg$Feature, trnatable[,1]),3]
#acceptorType <- factor(acceptorType, levels = sort(unique(acceptorType)))
#featnames = unique(as.character(coveragemeltagg$Feature))
##tails = as.numeric(unlist(lapply(strsplit(featnames, "-", fixed = TRUE), tail, 1)))
##anticodonname = as.character(unlist(lapply(strsplit(featnames, "-", fixed = TRUE), function(x) { return( x[length(x) - 1] ) })))
##aminoname = as.character(unlist(lapply(strsplit(featnames, "-", fixed = TRUE), function(x) { return( x[length(x) - 2] ) })))
#coveragemeltagg$Feature = factor(as.character(coveragemeltagg$Feature), levels = featnames)
#coveragemelt <- coveragemeltagg

coveragemelt <- coverageprep(coveragemeltagg, samples, trnatable)

acceptorType = trnatable[match(coveragemelt$Feature, trnatable[,1]),3]
print(unique(acceptorType))


#acceptorType <- factor(acceptorType, levels = sort(unique(acceptorType)))
acceptorType <- factor(acceptorType, levels = allaminos)
sortacceptor <- acceptorType[order(coveragemelt$variable, coveragemelt$Sample,-as.numeric(coveragemelt$Feature))]


endsmeltagg  <- aggregate(coverageall$readstarts, by=list(Feature = coverageall$Feature, Sample = sampletable[match(coverageall$Sample,sampletable[,1]),2], variable = coverageall$position), FUN=mean)
endsmelt <- coverageprep(endsmeltagg, samples, trnatable)

#endsmeltagg  <- aggregate(coverageall$readstarts, by=list(Feature = coverageall$Feature, Sample = sampletable[match(coverageall$Sample,sampletable[,1]),2], variable = coverageall$position), FUN=mean)
#endsmelt <- coverageprep(endsmeltagg, samples, trnatable)


 
#print(coverageall$readstarts)
readends  = cbind(coverageall$readstarts, coverageall$readends, coverageall$coverage - (coverageall$readstarts + coverageall$readends))

colnames(readends) =  c("Read Starts", "Read ends","Coverage")

readsmeltagg  <- aggregate(readends, by=list(Feature = coverageall$Feature, Sample = sampletable[match(coverageall$Sample,sampletable[,1]),2], variable = coverageall$position), FUN=mean)
readsmelt <- coverageprep(readsmeltagg, samples, trnatable)

#print(head(readsmelt))

allreadsmelt = melt(readsmelt, id.vars = c("Feature", "Sample", "variable"))


colnames(allreadsmelt) =  c("Feature", "Sample","position","endtype","value")
#print(head(allreadsmelt))
pcount = 30


deletemeltagg  <- aggregate(coverageall$deletions / (coverageall$coverage + pcount), by=list(Feature = coverageall$Feature, Sample = sampletable[match(coverageall$Sample,sampletable[,1]),2], variable = coverageall$position), FUN=mean)
deletemelt <- coverageprep(deletemeltagg, samples, trnatable)


#deletemelt$variable = factor(deletemelt$variable, levels=trnapositions)


    
    
#write.table( coverageall[coverageall$deletions / (coverageall$coverage + pcount) > .1,], file = paste(opt$directory,"/mismatch/",runname, "-sigdelete.txt",sep= "")) 

#write.table( deletemelt, file = paste(opt$directory,"/mismatch/",runname, "-alldelete.txt",sep= "")) 


mismatchesmeltagg <- aggregate(coverageall$mismatchedbases / (coverageall$coverage + pcount), by=list(Feature = coverageall$Feature, Sample = sampletable[match(coverageall$Sample,sampletable[,1]),2], variable = coverageall$position), FUN=mean)
mismatchesmelt <- coverageprep(mismatchesmeltagg, samples, trnatable)
#mismatchesmelt$variable = factor(mismatchesmelt$variable, levels=trnapositions)

write.table(coverageall[coverageall$mismatchedbases / (coverageall$coverage + pcount) > .1,],file = paste(opt$directory,"/mismatch/",runname, "-sigmismatch.txt",sep= ""))



coverageunique = coverageall[,c("Feature", "Sample","position","uniquecoverage","multitrnacoverage","multianticodoncoverage","multiaminocoverage")]

colnames(coverageunique)[colnames(coverageunique) == "uniquecoverage"]          <- "Transcript Specific"
colnames(coverageunique)[colnames(coverageunique) == "multitrnacoverage"]  <- "Isodecoder Specific"  
colnames(coverageunique)[colnames(coverageunique) == "multianticodoncoverage"]  <- "Isotype Specific"
colnames(coverageunique)[colnames(coverageunique) == "multiaminocoverage"]      <- "Not Amino Specific" 
allmultmelt = melt(coverageunique, id.vars = c("Feature", "Sample", "position"))

allmultmeltagg <- aggregate(allmultmelt$value, by=list(Feature = allmultmelt$Feature, Sample = sampletable[match(allmultmelt$Sample,sampletable[,1]),2], maptype = allmultmelt$variable,variable = allmultmelt$position), FUN=mean)
colnames(allmultmeltagg)[colnames(allmultmeltagg) == "x"]  <- "value"
allmultmeltagg$Sample <- factor(allmultmeltagg$Sample,levels = unique(sampletable[,2]), ordered = TRUE)
allmultmelt <- allmultmeltagg

print("::))")
print(head(coveragemelt))
makecombplot(coveragemelt,filename=paste(combinedfile,sep= ""))
#makecombplot(coveragemelt,filename=paste(uniquename, "-combinedfreecoverages.pdf",sep= ""), indepscales = TRUE)
print("::))")

modomicstable <- data.frame(trna = character(), mod = character(), pos = character(),stringsAsFactors=FALSE)

modomicstable <- read.table(text ="",col.names = c("trna", "mod", "pos"),colClasses = c("character", "character", "character")) #(trna = character(), mod = character(), pos = character(),stringsAsFactors=FALSE)
if(!is.null(opt$modomics) && file.exists(opt$modomics) ){
    modomicstable <- read.table(opt$modomics, header = TRUE)
    modomicstable$pos <- paste("X",modomicstable$pos, sep = "")
}



modomicstable <- modomicstable[as.character(modomicstable$pos) %in% unique(as.character(coveragemelt$variable)),]
modomicstable$dist <- match(modomicstable$pos,  levels(coveragemelt$variable)) 
modomicstable$Feature <- factor(modomicstable$trna,levels = levels(coveragemelt$Feature) )

stopmods = c("m1A","m2,2G","m1G","m1I","m3C")
modomicstable <- modomicstable[modomicstable$mod %in% stopmods,]

colnames(modomicstable)[colnames(modomicstable) == 'mod'] <- 'Modification'
#() # 1 minute


#print("**||2")
canvas = ggplot(coveragemelt,aes(x=variable,y=value), size = 2)
if(max(coveragemelt$value) <= 1.1){
#allcoverages <- canvas + theme_bw() + geom_bar(stat="identity") + facet_grid(Feature ~ Sample, scales="free") + geom_vline(aes(xintercept = dist, col = Modification),size=.4,linetype = "longdash",data = modomicstable,show.legend=TRUE)#+theme(axis.text.y=element_text(colour="black",size=6),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8), strip.text.y = element_text(angle=0),strip.text.x = element_text(size = ,angle=0))+ ylab("Normalized Read Count") + xlab("tRNA position") +   scale_y_continuous(breaks=percentbreaks,limits = c(0, 1.01)) +scale_x_discrete(breaks=c("X1","X9","X26","X37","X44","X58","X65","X73"), labels=c(fiveprimename,"m1g","m22g","anticodon","varloop","m1a","65",ccatailname))   #+scale_x_discrete(breaks=c("X1","X37","X73"), labels=c(fiveprimename,"anticodon", ccatailname))
}else{
allcoverages <- canvas + theme_bw() + geom_bar(stat="identity") + facet_grid(Feature ~ Sample, scales="free", space = "fixed")  +  geom_vline(aes(xintercept = dist, col = Modification),size=.4,linetype = "longdash",data = modomicstable,show.legend=TRUE)#+theme(axis.text.y=element_text(colour="black",size=6),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8), strip.text.y = element_text(angle=0),strip.text.x = element_text(size = ,angle=0))+ ylab("Normalized Read Count") + xlab("tRNA position") +   scale_y_continuous(breaks=myBreaks) +scale_x_discrete(breaks=c("X1","X9","X26","X37","X44","X58","X65","X73"), labels=c(fiveprimename,"m1g","m22g","anticodon","varloop","m1a","65",ccatailname))   #+scale_x_discrete(breaks=c("X1","X37","X73"), labels=c(fiveprimename,"anticodon", ccatailname))
}

#scalefactor*(2 + 1*length(unique(coveragemelt$Feature)))
#scalefactor*(2+4*length(unique(coveragemelt$Sample)))
#q()

ggsave(filename=outputfile, allcoverages,height=scalefactor*(2 + 1*length(unique(coveragemelt$Feature))),width=scalefactor*(2+4*length(unique(coveragemelt$Sample))), limitsize=FALSE, dpi = 600) #real one

#ggsave(filename=outputfile, allcoverages)

#q() #9 mins

#print("**||3")

if(!is.null(opt$directory) && FALSE){
for (curramino in unique(acceptorType)){
aminodata = coveragemelt[acceptorType == curramino,]
aminomodomicstable = modomicstable[modomicstable$Feature %in% unique(aminodata$Feature),]

aminocoverage <- ggplot(aminodata,aes(x=variable,y=value), size = 2) + theme_bw()+facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity") +  geom_vline(aes(xintercept = dist, col = Modification),data = aminomodomicstable,show.legend=TRUE)+theme(axis.text.y=element_text(colour="black",size=6),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8), strip.text.y = element_text(angle=0),strip.text.x = element_text(size = ,angle=0))+ ylab("Normalized Read Count") + xlab("tRNA position") + scale_y_continuous(breaks=myBreaks) +scale_x_discrete(breaks=c("X1","X9","X26","X37","X44","X58","X65","X73"), labels=c(fiveprimename,"m1g","m22g","anticodon","varloop","m1a","65",ccatailname)) + guide_legend(title="RNA modifications")   #+scale_x_discrete(breaks=c("X1","X37","X73"), labels=c(fiveprimename,"anticodon", ccatailname))

aminoname = paste(opt$directory,"/",curramino,"-coverage",outputformat, sep = "")

ggsave(filename=aminoname, aminocoverage,height=2 + scalefactor*1.5*length(unique(aminodata$Feature)),width=scalefactor*5*length(unique(aminodata$Sample)), limitsize=FALSE, dpi = 600)
aminopngname = paste(opt$directory,"/",curramino,"-coverage",".png", sep = "")
scalefactor = 4
ggsave(filename=aminopngname, aminocoverage,height=2 + scalefactor*1.5*length(unique(aminodata$Feature)),width=scalefactor*5*length(unique(aminodata$Sample)), limitsize=FALSE)

}
}


if(!is.null(uniquename)){


#multamino <- read.table(paste(uniquename, "-multaminocoverages.txt",sep= ""), header = TRUE)
#multactable <- read.table(paste(uniquename, "-multaccoverages.txt",sep= ""), header = TRUE)
#multtrnas <- read.table(paste(uniquename, "-multtrnacoverages.txt",sep= ""), header = TRUE)
#uniquetable <- read.table(paste(uniquename, "-uniquecoverages.txt",sep= ""), header = TRUE)

#uniquetable$maptype = "Transcript Specific"
#multactable$maptype = "Isotype Specific"
#multamino$maptype = "Not Amino Specific"
#multtrnas$maptype = "Isodecoder Specific" 
#allmulttables <- rbind(uniquetable,multactable,multamino,multtrnas)

#allmulttables <- rbind(multamino,multtrnas,multactable,uniquetable)
#allmulttables <- rbind(uniquetable, multtrnas,multactable,multamino)



#allmulttables <- allmulttables[ , colSums(is.na(allmulttables)) < nrow(allmulttables)/8]
#allmultmelt = melt(allmulttables, id.vars = c("Feature", "Sample", "maptype"))
#allmultmelt[is.na(allmultmelt)] <- 0



#print(head(coverageall))

#print(head(coverageunique))




#colnames(allmultmelttest)[colnames(allmultmelttest) == "position"]  <- "variable"

#print(head(allmultmelt))

#print(head(allmultmelttest))


#allmultmeltagg <- aggregate(allmultmelt$value, by=list(Feature = allmultmelt$Feature, Sample = sampletable[match(allmultmelt$Sample,sampletable[,1]),2], maptype = allmultmelt$variable,variable = allmultmelt$position), FUN=mean)
#allmultmeltagg = allmultmeltagg[allmultmeltagg$variable %in% trnapositions,]
#colnames(allmultmeltagg)[colnames(allmultmeltagg) == "x"]  <- "value"
#allmultmeltagg$Sample <- factor(allmultmeltagg$Sample,levels = unique(sampletable[,2]), ordered = TRUE)
#allmultmelt <- allmultmeltagg

#allmultmelt$maptype <- factor(allmultmelt$maptype, levels=c("Not Amino Specific","Isotype Specific","Isodecoder Specific","Transcript Specific"))
#allmultmelt$variable = factor(allmultmelt$variable, levels=trnapositions)

#allmultmelt <- allmultmelt[order(-as.numeric(allmultmelt$maptype)),]


allmultmelt <- allmultmelt[order(-as.numeric(allmultmelt$maptype), match(allmultmelt$variable,trnapositions)),]
rownames(allmultmelt) <- c()
#print(unique(allmultmeltagg$variable))

#print(allmultmelt[allmultmelt$Feature == "tRNA-Val-AAC-1" & allmultmelt$Sample == "Truseq_AlkB",])

#print(head(allmultmelt))

featnames = unique(as.character(allmultmelt$Feature))
tails = as.numeric(unlist(lapply(strsplit(featnames, "-", fixed = TRUE), tail, 1)))
anticodonname = as.character(unlist(lapply(strsplit(featnames, "-", fixed = TRUE), function(x) { return( x[length(x) - 1] ) })))
aminoname = as.character(unlist(lapply(strsplit(featnames, "-", fixed = TRUE), function(x) { return( x[length(x) - 2] ) })))
featnames = featnames[order(aminoname, anticodonname,tails)]
allmultmelt$Feature = factor(as.character(allmultmelt$Feature), levels = featnames)

featnames = unique(as.character(endsmelt$Feature))
tails = as.numeric(unlist(lapply(strsplit(featnames, "-", fixed = TRUE), tail, 1)))
anticodonname = as.character(unlist(lapply(strsplit(featnames, "-", fixed = TRUE), function(x) { return( x[length(x) - 1] ) })))
aminoname = as.character(unlist(lapply(strsplit(featnames, "-", fixed = TRUE), function(x) { return( x[length(x) - 2] ) })))
featnames = featnames[order(aminoname, anticodonname,tails)]
endsmelt$Feature = factor(as.character(endsmelt$Feature), levels = featnames)

featnames = unique(as.character(mismatchesmelt$Feature))
tails = as.numeric(unlist(lapply(strsplit(featnames, "-", fixed = TRUE), tail, 1)))
anticodonname = as.character(unlist(lapply(strsplit(featnames, "-", fixed = TRUE), function(x) { return( x[length(x) - 1] ) })))
aminoname = as.character(unlist(lapply(strsplit(featnames, "-", fixed = TRUE), function(x) { return( x[length(x) - 2] ) })))
featnames = featnames[order(aminoname, anticodonname,tails)]
mismatchesmelt$Feature = factor(as.character(mismatchesmelt$Feature), levels = featnames)

featnames = unique(as.character(deletemelt$Feature))
tails = as.numeric(unlist(lapply(strsplit(featnames, "-", fixed = TRUE), tail, 1)))
anticodonname = as.character(unlist(lapply(strsplit(featnames, "-", fixed = TRUE), function(x) { return( x[length(x) - 1] ) })))
aminoname = as.character(unlist(lapply(strsplit(featnames, "-", fixed = TRUE), function(x) { return( x[length(x) - 2] ) })))
featnames = featnames[order(aminoname, anticodonname,tails)]
deletemelt$Feature = factor(as.character(deletemelt$Feature), levels = featnames)


#q()
print("**||2")
#amino specific plots

#print(head(allmultmelt))

for (curramino in unique(acceptorType)){

aminodata = allmultmelt[acceptorType == curramino,]
aminomodomicstable = modomicstable[modomicstable$Feature %in% unique(aminodata$Feature),]



aminonamesec = paste(uniquename, "-",curramino,"_cov",outputformat,sep= "")
makecovplot(aminodata,aminonamesec)


aminoendsdata = endsmelt[acceptorType == curramino,]
aminonamesec = paste(opt$directory,"/mismatch/",runname, "-",curramino,"_fiveprimeends",outputformat,sep= "")
makebasiccovplot(aminoendsdata,aminonamesec)


readsmeltdata = allreadsmelt[acceptorType == curramino,]
#print(head(readsmeltdata))
readsnamesec = paste(opt$directory,"/",runname, "-",curramino,"_endplot",outputformat,sep= "")

#makeendplot(readsmeltdata,readsnamesec)




aminomismatchdata = mismatchesmelt[acceptorType == curramino,]
aminonamemissec = paste(opt$directory,"/mismatch/",runname, "-",curramino,"_mismatch",outputformat,sep= "")
makepercentcovplot(aminomismatchdata,aminonamemissec)

aminodeletedata = deletemelt[acceptorType == curramino,]
aminonamedelsec = paste(opt$directory,"/mismatch/",runname, "-",curramino,"_delete",outputformat,sep= "")
makepercentcovplot(aminodeletedata,aminonamedelsec)
#write.table(aminodeletedata[aminodeletedata$value > .5,], file = paste(opt$directory,"/mismatch/",runname, "-",curramino,"delete.txt",sep= "")) 

#locusplotname = paste(opt$directory,"/pretRNAs/",runname,"-pretRNAcoverage.pdf", sep = "")
#unique only plots
aminonameunique = paste(opt$directory,"/unique/",opt$directory, "-",curramino,"_uniqueonlycov",outputformat,sep= "")
#print(head(aminodata[aminodata$maptype == "Transcript Specific",]))
makecovplot(aminodata[aminodata$maptype == "Transcript Specific",],aminonameunique)

}
#change this back
#tRNA specific plots
if(FALSE){
for (currtranscript in unique(allmultmelt$Feature)){

transcriptdata = allmultmelt[allmultmelt$Feature == currtranscript,]
transcriptmodomicstable = modomicstable[modomicstable$Feature %in% unique(transcriptdata$Feature),]
 


transcriptname = paste(opt$directory,"/indiv/",currtranscript,"-uniqcoverage",outputindivformat, sep = "")
#makecovplot(transcriptdata,transcriptname)

}
}

makecovplot(allmultmelt,paste(uniquename, "-coverages.pdf",sep= ""))




}







#print(head(locicoverages))

locicoveragesagg <- aggregate(locicoverages$coverage, by=list(Feature = locicoverages$tRNA_name, Sample = sampletable[match(locicoverages$sample,sampletable[,1]),2], variable = locicoverages$position), FUN=mean)
colnames(locicoveragesagg)[colnames(locicoveragesagg) == "x"]  <- "value"
locicoveragesagg$Sample <- factor(locicoveragesagg$Sample,levels = unique(sampletable[,2]), ordered = TRUE)





locicoveragesagg = locicoveragesagg[locicoveragesagg$variable %in% locuspositions,]
locicoveragesagg$variable = factor(locicoveragesagg$variable, levels=locuspositions)


#print(unique(locicoveragesagg$variable))
#print("**")
#factor

oldpos = unique(locicoveragesagg$variable)

locicoverages <- locicoveragesagg




#vector of amino acids
locustable = expand.delimited(trnatable,3,2)
acceptorType = locustable[match(locicoverages$Feature, locustable[,2]),1]


acceptorType <- factor(acceptorType, levels = allaminos)
#sortacceptor <- acceptorType[order(locicoverages$variable, locicoverages$Sample,-as.numeric(locicoverages$Feature))]


sortcovmelt <- locicoverages[order(locicoverages$variable, locicoverages$Sample,as.numeric(acceptorType),locicoverages$Feature),]  
sortacceptor <- acceptorType[order(locicoverages$variable, locicoverages$Sample,as.numeric(acceptorType),locicoverages$Feature)]

#print(unique(sortcovmelt$variable))

#print(setdiff(oldpos, levels(sortcovmelt$variable)))

#print(setdiff(levels(sortcovmelt$variable), oldpos))
#print("***")
#print(head(acceptorType))
#print(allaminos)
#print(unique(sortacceptor))
#print(unique(locicoverages$Sample))
#print(unique(sortcovmelt$variable))
print("***||&&")
locuscombinedname = paste(opt$directory,"/pretRNAs/",runname,"-pretRNAcombinedcoverage.pdf", sep = "")
makelocuscombplot(sortcovmelt, locuscombinedname) 
#covsummary <- ggplot(sortcovmelt,aes(x=variable,y=value, fill = sortacceptor, order = as.numeric(sortacceptor))) + facet_grid( ~ Sample, scales="free") +xlab("Position")+ geom_bar(stat="identity")+theme_bw()+theme(axis.text.y=element_text(colour="black",size=8), strip.text.y = element_text(angle=0,size=4),strip.text.x = element_text(angle=0,size=8),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8))+ ylab("Read Share") +   scale_y_continuous(breaks=myBreaks, labels = c("0","1")) +scale_x_discrete(breaks=c("1","13","22","31","39","53","61","73"), labels=c(fiveprimename,dstartname,dendname,acstartname,acendname,tstartname,tendname,ccatailname)) +scale_fill_discrete(drop=FALSE, name="Acceptor\ntype", breaks = levels(sortacceptor))
#ggsave(filename=locuscombinedname,covsummary)
print("***||")
allmultmelt <- sortcovmelt



#print(locicoverages[locicoverages$Feature == "tRNA-Ala-AGC-1-1" & locicoverages$Sample == "MCF7_Truseq_1",])

#print(head(locicoverages))

locusplot <- ggplot(locicoverages,aes(x=variable,y=value)) + facet_grid(Feature ~ Sample, scales="free") + geom_bar(stat="identity")+theme_bw()+theme(axis.title.x=element_blank(), axis.text.y=element_text(colour="black",size=4),axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=8), strip.text.y = element_text(angle=0,size=6),strip.text.x = element_text(angle=0,size=8))+ ylab("read count") +scale_x_discrete(breaks=c("1","13","22","31","39","53","61","73"), labels=c(fiveprimename,dstartname,dendname,acstartname,acendname,tstartname,tendname,ccatailname))
locusplotname = paste(opt$directory,"/pretRNAs/",runname,"-pretRNAcoverage.pdf", sep = "")

ggsave(locusplot, filename=locusplotname,height=.5*length(unique(locicoverages$Feature)),width=2*length(unique(locicoverages$Sample)), limitsize=FALSE)





