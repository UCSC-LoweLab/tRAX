#!/usr/bin/env bash

# Supported genomes
GENOMES=("hg19" "hg38" "rn6" "mm10" "sacCer3")

# Help function
function print_usage() {
  echo "USAGE: $0 databasename location(optional)" >&2
  echo "  databasename: ${GENOMES[@]}" >&2
  echo "  location: Directory to store the database (default = /rnadb)" >&2
}

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
       GTF_URL="ftp://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.gtf.gz"
       gtRNAdb_URL="http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz"
       gtRNAdb_OUT="hg38-tRNAs-detailed.out"
       gtRNAdb_NAME="hg38-tRNAs_name_map.txt"
       GENOME_URL="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
       FASTA=true
 elif test "${1}" = "mm10"
<<<<<<< HEAD
 then
    GTF_URL="ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz"
    gtRNAdb_URL="http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/mm10-tRNAs.tar.gz"
    gtRNAdb_OUT="mm10-tRNAs-confidence-set.out"
    gtRNAdb_NAME="mm10-tRNAs_name_map.txt"
    GENOME_URL="https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit"
    FASTA=false
=======
   then
      GTF_URL="ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz"
      gtRNAdb_URL="http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/mm10-tRNAs.tar.gz"
      gtRNAdb_OUT="mm10-tRNAs-confidence-set.out"
      gtRNAdb_NAME="mm10-tRNAs_name_map.txt"
      GENOME_URL="https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit"
      FASTA=false
 elif test "${1}" = "rn6"
   then
      GTF_URL="ftp://ftp.ensembl.org/pub/release-102/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.102.gtf.gz"
      gtRNAdb_URL="http://gtrnadb.ucsc.edu/genomes/eukaryota/Rnorv6/rn6-tRNAs.tar.gz"
      gtRNAdb_OUT="rn6-tRNAs-detailed.out"
      gtRNAdb_NAME="rn6-tRNAs_name_map.txt"
      GENOME_URL="https://hgdownload.soe.ucsc.edu/goldenPath/rn6/bigZips/rn6.fa.gz"
      FASTA=true
>>>>>>> 58611287c18ef7395448b36f3f1e32029cbf4fe9
 elif test "${1}" = "sacCer3"
   then
      GTF_URL="ftp://ftp.ensembl.org/pub/release-97/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.97.gtf.gz"
      gtRNAdb_URL="http://gtrnadb.ucsc.edu/genomes/eukaryota/Scere3/sacCer3-tRNAs.tar.gz"
      gtRNAdb_OUT="sacCer3-tRNAs.out-noChrM"
      gtRNAdb_NAME="sacCer3-tRNAs_name_map.txt"
      GENOME_URL="https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.2bit"
      FASTA=false
  else
      echo "Could not generate RNA database, non-supported DB or parameter missing"
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
    > ${DB_LOCATION}/genes.gtf
  echo "Generating GTF Done"

  # gtRNAdb Files
  echo "Generating gtRNAdb"
  wget -q -O ${DB_LOCATION}/tse.tar.gz ${gtRNAdb_URL}
  tar zxf ${DB_LOCATION}/tse.tar.gz -C ${DB_LOCATION}
  rm ${DB_LOCATION}/tse.tar.gz
  echo "Generating gtRNAdb Done"

  # Genome Fasta File from UCSC
  echo "Generating Fasta"
  if test ${FASTA} = true
    then
      wget -q -O - ${GENOME_URL} | gzip -cd > ${DB_LOCATION}/genome.fa
  else
    wget -q -O ${DB_LOCATION}/genome.2bit ${GENOME_URL}
    wget -q -O ${DB_LOCATION}/twoBitToFa http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
    chmod +x ${DB_LOCATION}/twoBitToFa
    ${DB_LOCATION}/twoBitToFa ${DB_LOCATION}/genome.2bit ${DB_LOCATION}/genome.fa
  fi
  echo "Generating Fasta Done"

  # TRAX maketrnadb
  echo "Starting TRAX makernadb"
  maketrnadb.py \
    --databasename=${DB_LOCATION}/db \
    --genomefile=${DB_LOCATION}/genome.fa \
    --trnascanfile=${DB_LOCATION}/${gtRNAdb_OUT} \
    --namemapfile=${DB_LOCATION}/${gtRNAdb_NAME}
}

<<<<<<< HEAD
db_builder ${1}
=======
# Init test
if [ -z "$2" ]
  then
    DB_LOCATION="/rnadb"
else
  DB_LOCATION=${2}
  mkdir ${DB_LOCATION}
fi

if [ -z "$1" ]
  then
    print_usage
else
db_builder ${1} ${DB_LOCATION}
fi
>>>>>>> 58611287c18ef7395448b36f3f1e32029cbf4fe9
