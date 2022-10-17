#!/bin/bash
# Setup environment
conda activate kallistobus

# Fetch data
wget http://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/3.0.0/neuron_10k_v3/neuron_10k_v3_fastqs.tar /data/neuron_10k_v3_fastqs.tar
cd ./data
tar -xvf ./neuron_10k_v3_fastqs.tar

# Download reference kallisto index for mouse
mkdir index
cd index
kb ref -d mouse -i index.idx -g t2g.txt -f1 transcriptome.fasta

# Run kallisto and bustools (via kb-python)
cd ..
mkdir kallisto
cd kallisto
kb count --h5ad -i ../index/index.idx -g ../index/t2g.txt -x 10xv3 -o output --filter bustools -t 2 \
    ../neuron_10k_v3_fastqs/neuron_10k_v3_S1_L001_R1_001.fastq.gz \
    ../neuron_10k_v3_fastqs/neuron_10k_v3_S1_L001_R2_001.fastq.gz \
    ../neuron_10k_v3_fastqs/neuron_10k_v3_S1_L002_R1_001.fastq.gz \
    ../neuron_10k_v3_fastqs/neuron_10k_v3_S1_L002_R2_001.fastq.gz
