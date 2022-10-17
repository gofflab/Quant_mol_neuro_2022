#!/bin/bash

#SBATCH -J featureCounts # Job name
#SBATCH -o featureCounts.out # Standard output
#SBATCH -e featureCounts.err # Standard error

#SBATCH --partition=defq # rockfish partition
#SBATCH -c 4  # number of cores requested
#SBATCH --mem=2g
#SBATCH -t 10:00:00 # requested time, hh:mm:ss.

cd /home/thwang12/ME-440/taeyoung # change this to your working directory

for sample in SRR5715043 SRR5715044 SRR5715045 SRR5715063 SRR5715064 SRR5715065; do
    /home/thwang12/ME-440/taeyoung/Software/subread-2.0.3-source/bin/featureCounts ./2.Alignment/${sample}.Aligned.sortedByCoord.out.bam \
        -o ./3.Count/${sample}.fCounts \
        -a ./Genome/gencode.vM30.annotation.gtf \
        -s 2 \
        -T 4 -g gene_id -C &>./3.Count/${sample}.fCounts.out
done
