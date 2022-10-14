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

### Fastq files

Copy 6 fastq files to your FASTQ directory.

```sh
cd 1.Fastq  # assuming that you are currently at "Module_6" folder
cp /data/lgoff2/ME-440/taeyoung/1.Fastq/*[0-9].fastq.gz . # . means the current directory, you can modify it if necessary.
```

### Script files

Copy script files to your Script directory. 

```sh
cd Script # assuming that you are currently at "Module_6" folder
cp /data/lgoff2/ME-440/taeyoung/Script/* . # . means the current directory, you can modify it if necessary.
```

### Editor to modify sbatch or text files
Here, I use an editor, nano, but you can use your preferred editior, for example, Rstudio editor on rockfish. In case you use Rstudio server, just click the file name in Rstudio file panel to open it, and modify it. 


### How to use RStudio Server at portal.rockfish.jhu.edu

<img width="1340" alt="image" src="https://user-images.githubusercontent.com/34997334/195887394-088d935c-66a7-4df5-9a16-ec4c91a5dde8.png">

Use these settings

<img width="488" alt="image" src="https://user-images.githubusercontent.com/34997334/195887680-c8b918b7-96f4-4dd8-a917-bd5a9fb05ee6.png">


<img width="599" alt="image" src="https://user-images.githubusercontent.com/34997334/195887484-25775aa1-6446-42d7-a978-75fa4ce55eee.png">

Wait until the box turn green and click Launch RStudio server

<img width="717" alt="image" src="https://user-images.githubusercontent.com/34997334/195887791-2861a1a0-1769-416a-95b7-a18212a0dac5.png">



# 1. QC Fastq files

### (1) FastQC
To run Fastqc, we will use an interactive mode to use compute nodes at rockfish. An interactive mode is useful to perform small tasks and test workflows. Let's request 3 cores (-n 4) and 3GB RAM (-m 3G) and 1 hr working time (-t 1hr). Put this line on your terminal, not Rstudio terminal panel.

```sh
interact -n 4 -m 3g -t 01:00:00
```
Now you are at the compute nodes. The promot should be seen as `[thwang12@c714_teyoung]$`, not a log-in node like `[thwang12@login01 taeyoung]$`.  

**Run fastqc**
Note that 3 fastq files are processed at a time.

```sh
cd 1.Fastq # assuming that you are currently at "Module_6" folder
fastqc --help
fastqc -t 3 *.fastq.gz
exit
```


### (2) Trimming adapters

You can still use an interactive mode, but this time, we will submit a job to
use computational nodes. Make sure that you are at the login node, not terminal node.

First, check the batch script and modify it necessarily. **Modify the script according to your working
directory** (modify the line starting "cd" to change the directory to your working directory.)

```sh
nano ./Script/cutadapt.sbatch # assuming that you are currently at "Module_6" folder
```

Then, go to the folder and submit a job. Submitting a job should be done in login node. Use terminal, not Rstudio session at rockfish.

```sh
cd 1.Fastq # assuming that you are currently at "Module_6" folder
sbatch ../Script/cutadapt.sbatch
sacct # check the status of your job.
```

# 2. Alignment

To map RNA-seq reads on the genome, we will use STAR. You first need to create a
genome index then run STAR for fastq files.

### (1) Creating a genome index

**This steps takes time, so we will use a pre-made index in
`/data/lgoff2/ME-440/taeyoung/Genome/STAR_index` in the class. Skip this step. **

First, check the batch script and modify it necessarily. **Modify the script
according to your working directory** (modify the line starting "cd" to change the directory to your working directory.)

```sh
nano ./Script/star_index.sbatch # assuming that you are currently at "Module_6" folder
```

Then, go to the folder and submit a job. Submitting a job should be done in login node. Use terminal, not Rstudio session at rockfish.

```
cd Genome
sbatch ../Script/star_index.sbatch
sacct # check the status of your job.
```

### (2) Mapping

First, check the batch script and modify it necessarily. **Modify the script
according to your working directory** (modify the line starting "cd" to change the directory to your working directory.)

```sh
nano ./Script/star_align.sbatch # assuming that you are currently at "Module_6" folder
```

Then, go to the folder and submit a job. Submitting a job should be done in login node. Use terminal, not Rstudio session at rockfish.

```sh
cd 2.Alignment # assuming that you are currently at "Module_6" folder
sbatch ../Script/star_align.sbatch
sacct # check the status of your job.
```

# 3. Counting

### (1) featureCounts

First, check the batch script and modify it necessarily. **Modify the script
according to your working directory** (modify the line starting "cd" to change the directory to your working directory).  At rockfish, featureCounts is not
installed by default, so you need to install it first under your home directory.
But you can use the installed one that is available from my home directory. Note that this
script will use featureCounts installed under my home directory.

```sh
nano ./Script/featureCounts.sbatch # assuming that you are currently at "Module_6" folder
```

Then, go to the folder and submit a job. Submitting a job should be done in login node. Use terminal, not Rstudio session at rockfish.

```sh
cd 3.Count # assuming that you are currently at "Module_6" folder
sbatch ../Script/featureCounts.sbatch
sacct # check the status of your job.
```

### (2) Summary

In this step, we will use a R script to save count results in a single Rda file.
Use an interactive mode to run R script line by line or use an Rstudio session
through portal.rockfish.jhu.edu

#### Interactive mode
If you are still at an interactive mode that was initated during the fastqc step above, you don't have to run the following lines on your terminal.
```sh
interact -n 1 -m 3g -t 01:00:00
module load R  
``` 
Execute each line in `../Script/makeCountTable.R`. **Modify the script according
to your working directory** (line starting with "cd").

#### Rstudio on portal.rockfish.edu

After open an Rstudio session, execute each line in
`../Script/makeCountTable.R`. **Modify the script according to your working
directory** (line starting with "cd").

# 4. Differential expression analysis

I usually download the rda file that has generated in the above step to my
mac/PC and analyze it locally to avoid connection issue. You can also use
Rstudio session on portal.rockfish.edu. Follow `DE.Rmd`.
