#!/usr/bin/env bash

# Supported genomes
GENOMES=("hg19" "hg38" "rn6" "mm10" "sacCer3")

# Help function
function print_usage() {
  echo "USAGE: $0 tool databasename data" >&2
  echo "  tool: make, build, run, manual" >&2
  echo "    make: Build the docker container" >&2
  echo "    build: Build RNA database into a docker volume" >&2
  echo "    manual: Run container with prebuilt docker volume" >&2
  echo "  databasename: ${GENOMES[@]}" >&2
  echo "  data: Directory to mount with the data (optional)" >&2
}

# Function to build the Docker container
function docker_make() {
  docker build --no-cache -f Dockerfile -t trax .
}

# Function to start the container and build a RNA database
function docker_build_db() {
  if [ -z "$1" ]
    then
      echo "supply a databasename: ${GENOMES[@]}"
  else
    docker volume create rnadb-${1}
    docker run --rm -it --name trax-build-rnadb-${1} \
      -v `pwd`:/opt/trax \
      -v rnadb-${1}:/rnadb \
      trax \
      quickdb_trax.bash ${1}
  fi
}

# Function to start a manual Docker TRAX container
function docker_manual() {
  docker run --rm -it --name trax-${USER}-2 \
  --user=`id -u`:`id -g` \
  -v rnadb-${1}:/rnadb \
  trax
}

# Init function
if [[ ${1} = "make" ]]
then
 docker_make
elif [[ ${1} = "build" ]]
then
  [[ ${GENOMES[*]} =~ ${2} ]] && docker_build_db ${@:2} || print_usage
elif [[ ${1} = "manual" ]]
then
  [[ ${GENOMES[*]} =~ ${2} ]] && docker_manual ${@:2} || print_usage
else
  print_usage
fi
