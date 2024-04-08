# Project Proposal for Practical Computing Skills for Omics Data

### This project will focus on creating an assembly-first metagenomics workflow using samples that were taken from potatoes suspected to be infected with a bacterial pathogen. 

<br>

## Project Description
### In this project, there will be **five** key steps covered. They are: 
1. quality check raw reads
2. trim adapters and low quality reads
3. assembly the reads
4. map the reads back to the original genome
5. bin the reads, creating a metagenome-assembled genome (MAG)
<br>

It is expected that this project will create a workflow from raw metagenomic reads to metagenome-assembled genomes. Raw data will be in the form of fastq files, and after binning, MAGs will be in fasta file format.

Not only will this project produce MAGs, but also a well-defined, organized workflow for myself and others to utilize. As of right now, this project will require user interference at each step, however, with the knowledge gained in this course and further studying, it is a goal to automate this process with the submission of just one batch job to achieve MAGs.  
<br>
## Project Organization
The root directory for this project is located at:
``` bash
/fs/ess/PAS2700/htoth99/project
```
In this directory, you will find the following directories and subdirectories:
* scripts
* raw data
* results
    * binning
    * fastqc
    * meta_spades
    * read_map
    * trim_out
<br>

The results directory will contain the output for each of the five steps listed above in Project Description.

Each step will require external bioinformatic tools (outside bash), and the read mapping step is a series of steps. I plan to incorporate runner scripts, however, this aspect is still needed to be covered in the course. Therefore, I am uncertain this aspect will be integrated correctly.

Since this data is fairly large, slurm batch jobs will be utilized for each step. I would like to structure the code in such a way where the scripts remain unchanged after originally drafted, and inputs via the command line or runner scripts are the only changes being made. However, since this is a fairly new concept to me, I am uncertain that this aspect will make it into the final project.

<br>

## Technical Aspects
This project will utilize data produced from Dr. Jonathan Jacobs Infectious Emerging Ecology Laboratory from the Department of Plant Pathology. Specifically, it will use samples that were gathered from symptomatic potato tubers with an unknown causal agent. Samples were pair-end sequneced at 150 bp with a depth of 25 M reads with the NextSeq2000. 

For each step, the following programs are likely to be used:
- quality assessment -> fastqc
- read trimming -> trimmomatic
- assembly -> meta-spades
- read mapping -> bwa
- binning -> Maxbin2

It is a possibility that a program could be swapped with a similar one if technical issues arise.

<br>

## Reasoning Behind This Project
I chose this project because it relates directly to one of my objectives for my master's thesis. Although I am aware of most of these processes, this organized workflow would significantly benefit my understanding of both the MGS raw data to MAG pipeline and bash scripting practices. I am a user of hard coding when it comes to bash scripting, and I would like to step away from this practice. 

Lastly, as mentioned above, I would love to automate this process for other and future members of my lab. It is a goal of mine to make a raw MGS data -> MAGs pipeline that is user friendly and efficient. Creating an orgnaized workflow is a great first step in this direction.


