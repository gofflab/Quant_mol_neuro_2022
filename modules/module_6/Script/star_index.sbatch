#!/bin/bash

#SBATCH -J star_index # Job name
#SBATCH -o star_index.out # Standard output
#SBATCH -e star_index.err # Standard error

#SBATCH --partition=defq # rockfish partition
#SBATCH -c 6  # number of cores requested
#SBATCH --mem=50g
#SBATCH -t 10:00:00 # requested time, hh:mm:ss.

module load STAR
cd /data/lgoff2/ME-440/taeyoung/

STAR --runMode genomeGenerate \
--runThreadN 6 \
--genomeFastaFiles ./Genome/GRCm39.primary_assembly.genome.fa \
--sjdbGTFfile ./Genome/gencode.vM30.annotation.gtf \
--sjdbOverhang 99 \
--genomeDir ./Genome/STAR_index

