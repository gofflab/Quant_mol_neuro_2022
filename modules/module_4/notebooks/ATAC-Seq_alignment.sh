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
jq -r '.[].geo_samples | .[].samples | .[].experiments | .[].runs | .[].files.ftp | .[].url ' $GEO.json | xargs -n 1 -P 4 wget -c

# Use jq to extract the metadata for each sample and write to a .tsv file
jq -r '.[].geo_samples | .[].samples | .[] | [.accession, (.experiments | .[].runs | .[] | .files.ftp | .[] | .[]), (.attributes | .[])] | @tsv' $GEO.json > $GEO.tsv

# Some samples were sequenced more than once. To see which samples have more than one run, run the following command:
# jq -r '.[].geo_samples | .[].samples | .[].experiments | .[].runs | keys ' $GEO.json
# We will aggregate the reads for each sample when we align them below.

# Run FastQC on all fastq files (4 at a time)
# This too will take some time to run on all fastq files
echo $(ls *.fastq.gz) | xargs -n 1 -P 4 fastqc

# Step back to main directory
cd ../../../