#!/usr/bin/env bash

# Function to build the database
function db_builder() {
  if test "${1}" = "hg19"
  then
     GTF_URL="ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz"
     gtRNAdb_URL="http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi19/hg19-tRNAs.tar.gz"
     gtRNAdb_OUT="hg19-tRNAs-detailed.out"
     gtRNAdb_NAME="hg19-tRNAs_name_map.txt"
     GENOME_URL="http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
     FASTA=true
  elif test "${1}" = "hg38"
  then
     GTF_URL="ftp://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz"
     gtRNAdb_URL="http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz"
     gtRNAdb_OUT="hg38-tRNAs-detailed.out"
     gtRNAdb_NAME="hg38-tRNAs_name_map.txt"
     GENOME_URL="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
     FASTA=true
 elif test "${1}" = "mm10"
 then
    GTF_URL="ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz"
    gtRNAdb_URL="http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/mm10-tRNAs.tar.gz"
    gtRNAdb_OUT="mm10-tRNAs-confidence-set.out"
    gtRNAdb_NAME="mm10-tRNAs_name_map.txt"
    GENOME_URL="https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit"
    FASTA=false
 elif test "${1}" = "sacCer3"
 then
    GTF_URL="ftp://ftp.ensembl.org/pub/release-97/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.97.gtf.gz"
    gtRNAdb_URL="http://gtrnadb.ucsc.edu/genomes/eukaryota/Scere3/sacCer3-tRNAs.tar.gz"
    gtRNAdb_OUT="sacCer3-tRNAs.out-noChrM"
    gtRNAdb_NAME="sacCer3-tRNAs_name_map.txt"
    GENOME_URL="https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.2bit"
    FASTA=false
  else
    echo "Could not generate RNA database"
    return
  fi

  # GTF File from Ensembl
  echo "Generating GTF"
  wget -q -O - ${GTF_URL} | \
    gzip -cd | \
    grep -v '^#' | \
    awk '{print "chr" $0;}' | \
    sed 's/chrMT/chrM/g' | \
    grep -e Mt_rRNA -e Mt_tRNA -e miRNA -e misc_RNA -e rRNA -e snRNA -e snoRNA -e ribozyme -e sRNA -e scaRNA \
    > /rnadb/genes.gtf
  echo "Generating GTF Done"

  # gtRNAdb Files
  echo "Generating gtRNAdb"
  wget -q -O /rnadb/tse.tar.gz ${gtRNAdb_URL}
  tar zxf /rnadb/tse.tar.gz -C /rnadb
  rm /rnadb/tse.tar.gz
  echo "Generating gtRNAdb Done"

  # Genome Fasta File from UCSC
  echo "Generating Fasta"
  if test ${FASTA} = true
  then
    wget -q -O - ${GENOME_URL} | gzip -cd > /rnadb/genome.fa
  else
    wget -q -O /rnadb/genome.2bit ${GENOME_URL}
    wget -q -O /rnadb/twoBitToFa http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
    chmod +x /rnadb/twoBitToFa
    /rnadb/twoBitToFa /rnadb/genome.2bit /rnadb/genome.fa
  fi
  echo "Generating Fasta Done"

  # TRAX maketrnadb
  echo "Starting TRAX makernadb"
  maketrnadb.py \
    --databasename=/rnadb/db \
    --genomefile=/rnadb/genome.fa \
    --trnascanfile=/rnadb/${gtRNAdb_OUT} \
    --namemapfile=/rnadb/${gtRNAdb_NAME}

  # Change permissions
  #chmod -R 777 /rnadb/
}

db_builder ${1}
