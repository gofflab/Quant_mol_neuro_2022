
# 0. Set-up
### Directory
I recommend to organize a directory similar to that shown below by creating several sub-folders. It is just the best practice to use my script in the class. Of course, you can have your own structure for your project.  
```
taeyoung (at my home directory)  
+-- Software  
+-- Genome  
+-- Results  
+-- Script  
+-- 1.Fastq  
+-- 2.Alignment  
+-- 3.Count  
```

### Fastq files
Copy 6 fastq files to your FASTQ directory.
```sh
cp /home/thwang12/ME-440/taeyoung/1.Fastq/*[0-9].fastq.gz YOUR_FASTQ_DIRECOTRY
```

# 1. To run Fastqc
To run Fastqc, we will use an interactive session to use compute nodes at rockfish. An interactive session is useful to perform small tasks and test workflows. Let's request 3 cores (-n 4) and 3GB RAM (-m 3G) and 1 hr working time (-t 1hr).
```sh
interact -n 4 -m 3g -t 01:00:00
```

Now you are at the compute nodes. The promot should be seen as `[thwang12@c714_teyoung]$`, not a log-in node like `[thwang12@login01 taeyoung]$`.  

```sh
cd 1.Fastq  
module load fastqc  
fastqc --help  
fastqc -t 3 *.fastq.gz  
exit  
```
**Note that 3 fastq files are processed at a time.**

# (optional) Trimming adapters
You can still use an interactive mode, but this time, we will submit a job to use computational nodes.

First, check the batch script and modify it necessarily. Here, I use an editor, nano, but you can use your preferred editior, or you can edit it at your local machine and upload it.
```sh
nano ./Script/cutadapt.sbatch
```
Then, go to the folder and submit a job.  
```sh
cd 1.Fastq  
sbatch ../Script/cutadapt.sbatch
sacct # check the status of your job.
```

# 2. Alignment
To map RNA-seq reads on the genome, we will use STAR. You first need to create a genome index then run STAR for fastq files.

### (1) Creating a genome index 
**This steps takes time, so we will use a pre-made index in `/home/thwang12/ME-440/taeyoung/Genome/STAR_index` in the class.**  
First, check the batch script and modify it necessarily.
```sh
nano ./Script/star_index.sbatch
```
Then, go to the folder and submit a job.  
```
cd Genome  
sbatch ../Script/star_index.sbatch
sacct # check the status of your job.
```
### (2) Mapping
First, check the batch script and modify it necessarily.
```sh
nano ./Script/star_align.sbatch
```
Then, go to the folder and submit a job.  
```sh
cd 2.Alignment  
sbatch ../Script/star_align.sbatch
sacct # check the status of your job.
```
# 3. Counting

### (1) featureCounts
First, check the batch script and modify it necessarily.
```sh
nano ./Script/featureCounts.sbatch
```
Then, go to the folder and submit a job.  
```sh
cd 3.Count  
sbatch ../Script/featureCounts.sbatch
sacct # check the status of your job.
```

### (2) Summary
In this step, we will use a R script to save count results in a single Rda file.  
Use an interactive mode to run R script line by line or use an Rstudio session through portoal.rockfish.edu  

#### Interactive mode
```sh
interact -n 1 -m 3g -t 01:00:00
module load R  
``` 
Execute each line in `../Script/makeCountTable.R`.

#### Rstudio on portal.rockfish.edu
After open an Rstudio session, execute each line in `../Script/makeCountTable.R`.

# 5. Differential expression analysis
I usually download the rda file that has generated in the above step to my mac/PC and analyze it locally to avoid connection issue. You can also use Rstudio session on portal.rockfish.edu. Follow `DE.Rmd`.


