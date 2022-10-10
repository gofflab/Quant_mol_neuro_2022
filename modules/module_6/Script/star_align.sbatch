#!/bin/bash

#SBATCH -J star # Job name
#SBATCH -o star.out # Standard output
#SBATCH -e star.err # Standard error

#SBATCH --partition=defq # rockfish partition
#SBATCH -c 4  # number of cores requested
#SBATCH --mem=30g
#SBATCH -t 10:00:00 # requested time, hh:mm:ss.

module load STAR

cd /home/thwang12/ME-440/taeyoung

for sample in SRR5715043 SRR5715044 SRR5715045 SRR5715063 SRR5715064 SRR5715065
do
STAR --genomeDir ./Genome/STAR_Index \
     --runThreadN 4 \
     --readFilesCommand zcat \
     --readFilesIn ./1.Fastq/${sample}.trimmed.fastq.gz \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix ./2.Alignment/${sample}. &> ./2.Alignment/${sample}.star.out
done

