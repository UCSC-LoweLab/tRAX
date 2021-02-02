#!/bin/env bash


#Adapted this from https://gist.github.com/arq5x/b67196a46db5b63bee06
#chrHSCHR15_1_CTG4

# zcat ../hgtrnadb/hg19.gtf.gz | perl -pne "s/^chrGL(\d+)\.\d+/chrUn_gl\$1/" | perl -pne "s/^chrHS(CHR\S+)/lc(\$1)/e;" |  perl -pne "s/^chrHS(CHR\S+)/lc(\$1)/e;"gzip -c >../hgtrnadb/hg19fix.gtf.gz 

bedtools genomecov -bg -ibam $1 -strand +  | awk -v mincov=$2 '{ if ($4 >= mincov) {print $0"\t"mincov} else {print $0"\t"$4}}' | bedtools groupby -g 1,5 -c 1,2,3,4 -o first,first,last,max | cut -f 3-6 | awk -v mincov=$2 '$4 >= mincov'| awk '{print $1 "\t" $2 "\t" $3 "\tFEATURE" NR "_plus_" $4 "\t1000\t+"}'
bedtools genomecov -bg -ibam $1 -strand -  | awk -v mincov=$2 '{ if ($4 >= mincov) {print $0"\t"mincov} else {print $0"\t"$4}}' | bedtools groupby -g 1,5 -c 1,2,3,4 -o first,first,last,max | cut -f 3-6 | awk -v mincov=$2 '$4 >= mincov'| awk '{print $1 "\t" $2 "\t" $3 "\tFEATURE" NR "_minus_" $4 "\t1000\t-"}'

