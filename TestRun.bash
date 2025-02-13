#!/usr/bin/env bash


REALNAME=$(readlink -f $0)
SCRIPTDIR=$( cd "$( dirname "$REALNAME" )" && pwd )


#Download and remove adapters from small RNA sequencing studies
fastq-dump --gzip  SRR029131 
fastq-dump --gzip  SRR029124 

fastq-dump --gzip  SRR207111
fastq-dump --gzip  SRR207116 

"$SCRIPTDIR/trimadapters.py" --runname adapt1 --runfile $SCRIPTDIR/testadapt1.txt --firadapter TCGTATGCCGTCTTCT --singleend
 
"$SCRIPTDIR/trimadapters.py" --runname adapt2 --runfile $SCRIPTDIR/testadapt2.txt --firadapter TCGTATGCCGTCTTCT  --singleend

#Download and combine hg19 chromosomes 
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar xvf chromFa.tar.gz -O > hg19.fa

#Download ensembl GTF, change chromosome names

wget -q -O - ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz | gzip -cd | grep -v '^#' | awk '{print "chr" $0;}' | grep -e Mt_rRNA -e miRNA -e misc_RNA -e rRNA -e snRNA -e snoRNA -e ribozyme -e sRNA -e scaRNA  >hg19-genes.gtf


#get tRNA information
wget --no-check-certificate https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi19/hg19-tRNAs.tar.gz
tar xvf hg19-tRNAs.tar.gz




#Create the tRNA database
"$SCRIPTDIR/maketrnadb.py" --databasename=hg19 --genomefile=hg19.fa --trnascanfile=hg19-tRNAs-confidence-set.out --namemapfile=hg19-tRNAs_name_map.txt


#Map the tRNAreads
#"$SCRIPTDIR/processsamples.py" --experimentname=TestTrnas --databasename=hg19 --samplefile=${SCRIPTDIR}/TestSamples.txt --ensemblgtf=hg19-genes.gtf --olddeseq
"$SCRIPTDIR/processsamples.py" --experimentname=TestTrnas --databasename=hg19 --samplefile=${SCRIPTDIR}/TestSamples.txt --ensemblgtf=hg19-genes.gtf


