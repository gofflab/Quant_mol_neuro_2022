#!/usr/bin/env bash
#########
# Fetch raw RNA-Seq reads and sample metadata for study "Hipposeq: an RNA-seq based atlas of gene expression in excitatory hippocampal neurons", and pseudoalign to the mouse reference genome.
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74985
# https://elifesciences.org/articles/14997
# https://doi.org/10.7554/eLife.14997
#########

# install needed components (if needed)
#conda install -c bioconda MultiQC
#conda install -c bioconda fastqc
#conda install -c bioconda kallisto
#conda install -c bioconda igv
#conda install jq
#pip install ffq gget

# Create acceession variable and directory structure
GEO="GSE74985"
mkdir -p data/raw/${GEO}
cd data/raw/${GEO}

# Get .json file with metadata for all samples
# This operation sometimes fails with errors because of rate limiting at NCBI.  If it fails, just try again.
ffq $GEO > $GEO.json

# Use jq to extract the ftp links for the fastq files from the .json and pass to wget to download
# Note: this will download all fastq files, and will be approximately 60 Gb (24 * 2.5Gb) of raw data. Please be sure you have sufficient disk space available.
jq -r '.[].geo_samples | .[].samples | .[] | .experiments | .[].runs | .[] | .files.ftp | .[] | .url ' $GEO.json | xargs wget

# Use jq to extract the metadata for each sample and write to a .tsv file
jq -r '.[].geo_samples | .[].samples | .[] | [.accession, (.experiments | .[].runs | .[] | .files.ftp | .[] | .[]), (.attributes | .[])] | @tsv' $GEO.json > $GEO.tsv

# Run FastQC on all fastq files
fastqc *.fastq.gz

# Step back to main directory, build dir structure for metadata, and fetch the mouse reference transcriptome and annotation from GENCODE
cd ../../../
vM30_FASTA_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/gencode.vM30.transcripts.fa.gz"
vM30_GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/gencode.vM30.annotation.gtf.gz"
mkdir -p metadata/$GEO/kallisto_index
cd metadata/$GEO/kallisto_index
wget $vM30_FASTA_URL
wget $vM30_GTF_URL

# Build kallisto index from mouse reference transcriptome
# Note: this may take a while to complete and is memory intensive
kallisto index -i gencode.vM30.idx gencode.vM30.transcripts.fa.gz

kallisto inspect gencode.vM30.idx

# Step back to main directory, build dir structure for results, and run kallisto pseudoalign
cd ../../../
mkdir -p results/$GEO/kallisto
cd results/$GEO/kallisto
for fastq in $(ls ../../raw/$GEO/*.fastq.gz); 
    do 
        SAMPLE_NAME=$(basename $fastq | cut -d "_" -f 1)
        kallisto quant -i ../../../metadata/$GEO/kallisto_index/gencode.vM30.idx -o $SAMPLE_NAME --single -l 200 -s 20 $fastq; 
    done

# Run MultiQC to generate a project QC report
cd ../../../
multiQC .

