#!/usr/bin/env bash

# Supported genomes
GENOMES=("hg19" "hg38" "rn6" "mm10" "sacCer3" "hg19mito" "hg38mito" "mm10mito")

# Help function
function print_usage() {
  echo "USAGE: $0 databasename location(optional)" >&2
  echo "  databasename: ${GENOMES[@]}" >&2
  echo "  location: Directory to store the database (default = /rnadb)" >&2
}

# Function to build the database
function db_builder() {
  downloaddb=false
  echo "${2}"
  if test "${1}" = "hg19"; then
    GTF_URL="ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz"
    gtRNAdb_URL="http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi19/hg19-tRNAs.tar.gz"
    gtRNAdb_OUT="hg19-tRNAs-detailed.out"
    gtRNAdb_NAME="hg19-tRNAs_name_map.txt"
    GENOME_URL="http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
    FASTA=true
  elif test "${1}" = "hg38"; then
    GTF_URL="ftp://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.gtf.gz"
    gtRNAdb_URL="http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz"
    gtRNAdb_OUT="hg38-tRNAs-detailed.out"
    gtRNAdb_NAME="hg38-tRNAs_name_map.txt"
    GENOME_URL="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    FASTA=true
  elif test "${1}" = "mm10"; then
    GTF_URL="ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz"
    gtRNAdb_URL="http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/mm10-tRNAs.tar.gz"
    gtRNAdb_OUT="mm10-tRNAs-confidence-set.out"
    gtRNAdb_NAME="mm10-tRNAs_name_map.txt"
    GENOME_URL="https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit"
    FASTA=false
  elif test "${1}" = "mm10mito"; then
    GTF_URL="ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz"
    gtRNAdb_URL="http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/mm10-tRNAs.tar.gz"
    gtRNAdb_OUT="mm10-tRNAs-confidence-set.out"
    gtRNAdb_NAME="mm10-tRNAs_name_map.txt"
    GENOME_URL="https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit"
    FASTA=false
  elif test "${1}" = "rn6"; then
    GTF_URL="ftp://ftp.ensembl.org/pub/release-102/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.102.gtf.gz"
    gtRNAdb_URL="http://gtrnadb.ucsc.edu/genomes/eukaryota/Rnorv6/rn6-tRNAs.tar.gz"
    gtRNAdb_OUT="rn6-tRNAs-detailed.out"
    gtRNAdb_NAME="rn6-tRNAs_name_map.txt"
    GENOME_URL="https://hgdownload.soe.ucsc.edu/goldenPath/rn6/bigZips/rn6.fa.gz"
    FASTA=true
  elif test "${1}" = "sacCer3"; then
    GTF_URL="ftp://ftp.ensembl.org/pub/release-97/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.97.gtf.gz"
    gtRNAdb_URL="http://gtrnadb.ucsc.edu/genomes/eukaryota/Scere3/sacCer3-tRNAs.tar.gz"
    gtRNAdb_OUT="sacCer3-tRNAs.out-noChrM"
    gtRNAdb_NAME="sacCer3-tRNAs_name_map.txt"
    GENOME_URL="https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.2bit"
    FASTA=false
  elif test "${1}" = "hg19mito"; then
    GTF_URL="ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz"
    dburl="http://trna.ucsc.edu/tRAX/data/refdb/hg19mito.tar.gz"
    downloaddb=true
  elif test "${1}" = "hg38mito"; then
    GTF_URL="ftp://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/Homo_sapiens.GRCh38.102.gtf.gz"
    dburl="http://trna.ucsc.edu/tRAX/data/refdb/hg38mito.tar.gz"
    downloaddb=true
  elif test "${1}" = "mm10mito"; then
    GTF_URL="ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz"
    dburl="http://trna.ucsc.edu/tRAX/data/refdb/mm10mito.tar.gz"
    downloaddb=true
  else
    echo "Could not generate RNA database, non-supported DB or parameter missing"
    return
  fi

  # GTF File from Ensembl
  echo "Generating GTF"
  wget -q -O - ${GTF_URL} |
    gzip -cd |
    grep -v '^#' |
    awk '{print "chr" $0;}' |
    sed 's/chrMT/chrM/g' |
    grep -e Mt_rRNA -e Mt_tRNA -e miRNA -e misc_RNA -e rRNA -e snRNA -e snoRNA -e ribozyme -e sRNA -e scaRNA \
      >${2}/genes.gtf
  echo "Generating GTF Done"

  if test ${downloaddb} = true; then
    echo "Downloading TRAX db"
    wget -O ${2}/mitodb.tar.gz ${dburl}
    tar zxf ${2}/mitodb.tar.gz -C ${2}
    exit 0
  fi

  # gtRNAdb Files
  echo "Generating gtRNAdb"
  wget -q -O ${2}/tse.tar.gz ${gtRNAdb_URL}
  tar zxf ${2}/tse.tar.gz -C ${2}
  rm ${2}/tse.tar.gz
  echo "Generating gtRNAdb Done"

  # Genome Fasta File from UCSC
  echo "Generating Fasta"
  if test ${FASTA} = true; then
    wget -q -O - ${GENOME_URL} | gzip -cd >${DB_LOCATION}/genome.fa
  else
    wget -q -O ${2}/genome.2bit ${GENOME_URL}
    twoBitToFa ${DB_LOCATION}/genome.2bit ${DB_LOCATION}/genome.fa
  fi
  echo "Generating Fasta Done"

  # TRAX maketrnadb
  echo "Starting TRAX makernadb"
  maketrnadb.py \
    --databasename=${2}/db \
    --genomefile=${2}/genome.fa \
    --trnascanfile=${2}/${gtRNAdb_OUT} \
    --namemapfile=${2}/${gtRNAdb_NAME}
}

# Init test
if [ -z "$2" ]; then
  DB_LOCATION="/rnadb"
else
  DB_LOCATION=${2}
  mkdir -p ${2}
fi

if [ -z "$1" ]; then
  print_usage
else
  db_builder ${1} ${DB_LOCATION}
fi
