## tRAX

### [Official docs](http://trna.ucsc.edu/tRAX)

### About
tRNA Analysis of eXpression (tRAX) is a software package built for in-depth analyses of tRNA-derived small RNAs (tDRs), mature tRNAs, and inference of RNA modifications from high-throughput small RNA sequencing data. While the tRAX workflow adopts popular RNA sequencing data analysis methods, which includes adapter trimming of raw sequencing data, read alignment to the reference, transcript abundance estimation, and differential expression analysis among samples, it specifically consists of features designed to support special characteristics of tRNAs and tDRs. To ensure alignment of sequencing reads to tRNA transcripts, tRAX uses a custom-built reference database that not only includes the reference genome but also mature tRNA transcripts with the addition of 3′ CCA tail not encoded genomically. Unlike popular read counting tools that only consider or recommend for uniquely mapped reads for RNA sequencing, tRAX allows reads to be mapped to multiple transcripts and gene loci, which is necessary for conserved tRNA isodecoders (tRNAs with the same anticodon but different sequences in the gene body) and identical tRNA genes that are commonly found in eukayotes. Read coverage for each tRNA transcript is reported in four categories – transcript-specific, isodecoder-specific, isotype-specific, and non-specific – that provide precise results on the level of uniqueness of the aligned reads. Moreover, tRAX computes separate read counts of tRNA fragments that align at the 5′ end, 3′ end, and the middle region of tRNA transcripts to distinguish the abundance of different fragment types. Differential expression comparison across samples is performed using read counts for tRNA transcripts and tRNA fragments to provide better understanding of possible distinction in different tRNA isotypes or isodecoders and fragment types. In addition, tRAX measures the base frequency at each position aligned to the tRNA transcripts for estimating the mis-incorporations that may represent RNA modifications essential for function, stability, and regulation.

 To enable researchers to study the analysis results at multiple levels, tRAX presents over 50 types of visual images and data files such as read distribution summarizations, read coverage per tRNA isodecoder, volcano plots of differential expression comparison between samples, and misincorporation location charts. More details are available at [Understanding Output Results](http://trna.ucsc.edu/tRAX/outputs/) and the [Reference Manual](http://trna.ucsc.edu/tRAX/references/).

### System requirements
tRAX requires to be run on a Linux/Unix system with at least 8 cores and 16 GB memory. Due to the large size of sequencing data, we do not recommend using tRAX on a regular desktop or laptop. The following dependencies have been tested, use newer versions at your own risk.

#### Dependencies
* Python 2.7
* pysam 0.15.3
  * Older versions have a memory leak, make sure you have an updated version
* bowtie2 2.3.5
* samtools 1.9
* R 4.0.2 or higher
* Deseq2 R library 1.30.0 or higher
* getopt R library 1.20.3 or higher
* ggplot2 R library 3.3.0 or higher
* ggrepel R library 0.8.2 or higher
* gridextra R library 2.3 or higher
* reshape2 R library 1.4.4 or higher
* Infernal 1.1 or higher
* UMI-tools 1.1.2 or higher
* The trackhub generation flag requires:
  * ucsc-bedgraphtobigwig v377
  * ucsc-twobittofa v377
  * bedtools 2.29.2
* The TestRun.bash script requires:
  * SRA toolkit(fastq-dump) 2.8.0
  * cutadapt 1.18
  * seqprep 1.3.2


### Using Docker Image
To eliminate the need of installing dependencies, you can download the Docker image from our [DockerHub repository](https://hub.docker.com/r/ucsclowelab/trax) using the command
```
docker pull ucsclowelab/trax
```

### Using Conda Enviroment
In addition to Docker you can alternatively use a Conda environment using the command
```
conda env create -f trax_env.yaml
```

### [Quickstart and tutorial](http://trna.ucsc.edu/tRAX/#tutorial)
