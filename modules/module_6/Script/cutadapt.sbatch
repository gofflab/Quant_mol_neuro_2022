#!/bin/bash

#SBATCH -J cutadapt # Job name
#SBATCH -o cutadapt.out # Standard output
#SBATCH -e cutadapt.err # Standard error

#SBATCH --partition=defq # rockfish partition
#SBATCH -c 3  # number of cores requested
#SBATCH --mem=3g
#SBATCH -t 10:00:00 # requested time, hh:mm:ss.

module load cutadapt

cd /home/thwang12/ME-440/taeyoung/ # change this to your working directory

for sample in SRR5715043 SRR5715044 SRR5715045 SRR5715063 SRR5715064 SRR5715065; do
    cutadapt -j 3 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        -m 20 \
        -o ./1.Fastq/${sample}.trimmed.fastq.gz \
        ./1.Fastq/${sample}.fastq.gz &>./1.Fastq/${sample}.cutadapt.log
done
