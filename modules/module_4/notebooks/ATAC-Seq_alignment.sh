#!/usr/bin/env bash

#########
# Demo acquisition and alignment of ATAC-Seq reads and sample metadata for study "Layer-specific chromatin accessibility landscapes reveal regulatory networks in adult mouse visual cortex".
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87548
# https://elifesciences.org/articles/21883
# https://doi.org/10.7554/eLife.21883
#########

# install needed components (if needed)
#conda install -c bioconda fastqc bowtie2 igv
#conda install jq
#pip install ffq gget multiqc

# Stop script on any errors.
set -e

# Create acceession variable and directory structure
GEO="GSE87548"
mkdir -p data/raw/${GEO}
cd data/raw/${GEO}

# Get .json file with metadata for all samples
# This operation sometimes fails with errors because of rate limiting at NCBI.  If it fails, just try again.
ffq $GEO > $GEO.json

# Use jq to extract the ftp links for the fastq files from the .json and pass to wget to download (4 at a time)
# Note: this will download all fastq files, and will be approximately 60 Gb (24 * 2.5Gb) of raw data. Please be sure you have sufficient disk space available.
# Note: These samples are paired-end, so there will be 2 files for each sample to download.
jq -r '.[].geo_samples | .[].samples | .[].experiments | .[].runs | .[].files.ftp | .[].url ' $GEO.json | xargs -n 1 -P 8 wget -c

# Use jq to extract the metadata for each sample and write to a .tsv file
jq -r '.[].geo_samples | .[].samples | .[] | [.accession, (.experiments | .[].runs | .[] | .files.ftp | .[] | .[]), (.attributes | .[])] | @tsv' $GEO.json > $GEO.tsv

# Some samples were sequenced more than once. To see which SRX experiments have more than one run, run the following command:
# jq -r '.[].geo_samples | .[].samples | .[].experiments | [.[].accession, (.[].runs | keys)] ' $GEO.json
# To see which SRS samples have more than one run, run the following command:
# jq -r '.[].geo_samples | .[].samples | .[] | [.accession, (.experiments | .[].runs | keys)] ' $GEO.json
# We will aggregate the reads for each sample after we align them below.

# Run FastQC on all fastq files (4 at a time)
# This too will take some time to run on all fastq files
echo $(ls *.fastq.gz) | xargs -n 1 -P 8 fastqc

# Step back to main directory
cd ../../../

# Create directory structure for metadata
mkdir -p metadata/$GEO/bowtie2_idx
cd metadata/$GEO/bowtie2_idx

# Download the mouse genome and build bowtie2 index
#gget ref --ftp -w dna,gtf mus_musculus | xargs -n 1 -P 4 wget -c # Downloads the genome and gtf files for the most recent mouse genome release
#bowtie2-build Mus_musculus.GRCm39.dna.primary_assembly.fa.gz GRCm39 # Builds the bowtie2 index

# OR download a pre-built bowtie2 index
wget -c https://genome-idx.s3.amazonaws.com/bt/GRCm39.zip
unzip GRCm39.zip

# Step back to main directory
cd ../../../
# Build dir structure for results
mkdir -p results/$GEO/bowtie2
cd results/$GEO/bowtie2

# Run bowtie pseudoalignments
for RUN in $(jq -r '.[].geo_samples | .[].samples | .[].experiments | .[].runs | .[].accession' ../../../data/raw/$GEO/$GEO.json);
    do
        echo "Performing Bowtie2 alignment for $RUN"
        bowtie2 -p 4 -x ../../../metadata/$GEO/bowtie2_idx/GRCm39/GRCm39 -1 ../../../data/raw/$GEO/${RUN}_1.fastq.gz -2 ../../../data/raw/$GEO/${RUN}_2.fastq.gz 2>$RUN.log | samtools view -bS - > $RUN.bam;
    done;

# Run MultiQC to generate a project QC report
cd ../../../
multiqc .