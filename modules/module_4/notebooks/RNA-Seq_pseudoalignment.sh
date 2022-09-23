#!/usr/bin/env bash
#########
# Fetch raw RNA-Seq reads and sample metadata for study "Hipposeq: an RNA-seq based atlas of gene expression in excitatory hippocampal neurons", and pseudoalign to the mouse reference genome.
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74985
# https://elifesciences.org/articles/14997
# https://doi.org/10.7554/eLife.14997
#
# Please try and work through this script by placing it into a new empty directory and working through the steps to see how it works.
#
#########

# install needed components (uncomment and run if needed).  You will need to have anaconda installed.
#conda install -c bioconda fastqc kallisto igv
#conda install jq
#pip install ffq gget multiqc

# Stop script on any errors.
set -e

# Create acceession variable and directory structure
GEO="GSE74985"
mkdir -p data/raw/${GEO}
cd data/raw/${GEO}

# Get .json file with metadata for all samples
# This operation sometimes fails with errors because of rate limiting at NCBI.  If it fails, just try again.
ffq $GEO > $GEO.json

# Use jq to extract the ftp links for the fastq files from the .json and pass to wget to download (4 at a time)
# Note: this will download all fastq files, and will be approximately 60 Gb (24 * 2.5Gb) of raw data. Please be sure you have sufficient disk space available.
jq -r '.[].geo_samples | .[].samples | .[].experiments | .[].runs | .[].files.ftp | .[].url ' $GEO.json | xargs -n 1 -P 4 wget -c

# Use jq to extract the metadata for each sample and write to a .tsv file
jq -r '.[].geo_samples | .[].samples | .[] | [.accession, (.experiments | .[].runs | .[] | .files.ftp | .[] | .[]), (.attributes | .[])] | @tsv' $GEO.json > $GEO.tsv

# Run FastQC on all fastq files (4 at a time)
# This too will take some time to run on all fastq files
echo $(ls *.fastq.gz) | xargs -n 1 -P 4 fastqc

# Step back to main directory
cd ../../../

# Build dir structure for metadata including kallisto reference & index
mkdir -p metadata/$GEO/kallisto_index
cd metadata/$GEO/kallisto_index

# Fetch the mouse reference transcriptome and annotation from GENCODE
#vM30_FASTA_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/gencode.vM30.transcripts.fa.gz"
#vM30_GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M30/gencode.vM30.annotation.gtf.gz"
#wget $vM30_FASTA_URL
#wget $vM30_GTF_URL

# Build kallisto index from mouse reference transcriptome
# Note: this may take a while to complete and is memory intensive
#kallisto index -i gencode.vM30.idx gencode.vM30.transcripts.fa.gz
#kallisto inspect gencode.vM30.idx

# Download pre-built mouse kallisto index
IDX_URL="https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/ensembl-96/mus_musculus.tar.gz"
wget -c $IDX_URL
tar -xzvf mus_musculus.tar.gz # Expands into a new directory called mus_musculus

# Step back to main directory
cd ../../../
# Build dir structure for results
mkdir -p results/$GEO/kallisto
cd results/$GEO/kallisto

# Run kallisto pseudoalignments
# Note: here we are running kallisto with the --genomebam and --gtf options to generate a alignment (.bam) file from the pseudoaligned reads. This is actually the longest part of the pseudoalignment process, and will take a while to complete.
# We are doing this here to demonstrate the use of these options, and to produce a mapping of the reads that can be visualized in IGV.
# However, this is not necessary for the downstream analysis, and can generally be omitted to save time and disk space.
for fastq in $(ls ../../../data/raw/$GEO/*.fastq.gz); 
    do 
        SAMPLE_NAME=$(basename $fastq .fastq.gz);
        echo "Performing kalliso pseudoalignment for $SAMPLE_NAME";
        time kallisto quant -i ../../../metadata/$GEO/kallisto_index/mus_musculus/transcriptome.idx -o $SAMPLE_NAME --single -l 200 -s 20 $fastq 2>$SAMPLE_NAME/$SAMPLE_NAME.log; 
        #kallisto quant -i ../../../metadata/$GEO/kallisto_index/mus_musculus/transcriptome.idx --genomebam -g ../../../metadata/$GEO/kallisto_index/mus_musculus/Mus_musculus.GRCm38.96.gtf -o $SAMPLE_NAME --single -l 200 -s 20 $fastq 2>$SAMPLE_NAME/$SAMPLE_NAME.log; 
    done

# Run MultiQC to generate a project QC report
cd ../../../
multiqc .

