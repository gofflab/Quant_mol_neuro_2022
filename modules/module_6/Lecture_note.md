# 0. Set-up

### Directory

I recommend to organize a directory similar to that shown below by creating
several sub-folders. It is just the best practice to use my script in the class.
Of course, you can have your own structure for your project.

```
MODUE_6
+-- Software
+-- Genome
+-- Results
+-- Script
+-- 1.Fastq
+-- 2.Alignment
+-- 3.Count
```



# 1. Go to RStudio Server at portal.rockfish.jhu.edu

<img width="1340" alt="image" src="https://user-images.githubusercontent.com/34997334/195887394-088d935c-66a7-4df5-9a16-ec4c91a5dde8.png">

Use these settings

<img width="488" alt="image" src="https://user-images.githubusercontent.com/34997334/195887680-c8b918b7-96f4-4dd8-a917-bd5a9fb05ee6.png">


<img width="599" alt="image" src="https://user-images.githubusercontent.com/34997334/195887484-25775aa1-6446-42d7-a978-75fa4ce55eee.png">

Wait until the box turn green and click Launch RStudio server

<img width="717" alt="image" src="https://user-images.githubusercontent.com/34997334/195887791-2861a1a0-1769-416a-95b7-a18212a0dac5.png">



# 2. Go to the folder you created on Monday
### Fastq files

Copy 6 fastq files to your FASTQ directory.

```sh
cd 1.Fastq
cp /data/lgoff2/ME-440/taeyoung/1.Fastq/*[0-9].fastq.gz . # . means the current directory, you can modify it if necessary.
```

Run fastqc

```sh
cd 1.Fastq
fastqc --help
fastqc -t 3 *.fastq.gz
exit
```

**Note that 3 fastq files are processed at a time.**

# (optional) Trimming adapters

You can still use an interactive mode, but this time, we will submit a job to
use computational nodes.

First, check the batch script and modify it necessarily. Here, I use an editor,
nano, but you can use your preferred editior, or you can edit it at your local
machine and upload it. **Modify the script according to your working
directory.**

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

To map RNA-seq reads on the genome, we will use STAR. You first need to create a
genome index then run STAR for fastq files.

### (1) Creating a genome index

**This steps takes time, so we will use a pre-made index in
`/data/lgoff2/ME-440/taeyoung/Genome/STAR_index` in the class.**

First, check the batch script and modify it necessarily. **Modify the script
according to your working directory.**

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

First, check the batch script and modify it necessarily. **Modify the script
according to your working directory.**

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

First, check the batch script and modify it necessarily. **Modify the script
according to your working directory.** At rockfish, featureCounts is not
installed by default, so you need to install it first under your home directory.
But you can use the installed one that is available from my home directory. This
script will use featureCounts installed under my home directory.

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
Use an interactive mode to run R script line by line or use an Rstudio session
through portal.rockfish.jhu.edu

#### Interactive mode

Execute each line in `../Script/makeCountTable.R`. **Modify the script according
to your working directory.**

#### Rstudio on portal.rockfish.edu

After open an Rstudio session, execute each line in
`../Script/makeCountTable.R`. **Modify the script according to your working
directory.**

# 4. Differential expression analysis

I usually download the rda file that has generated in the above step to my
mac/PC and analyze it locally to avoid connection issue. You can also use
Rstudio session on portal.rockfish.edu. Follow `DE.Rmd`.
