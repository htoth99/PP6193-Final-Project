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

It is expected that this project will create a workflow from raw metagenomic reads to metagenome-assembled genomes. 
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

## Technical Aspects
This project will utilize data produced from Dr. Jonathan Jacobs Infectious Emerging Ecology Laboratory from the Department of Plant Pathology. Specifically, it will use samples that were gathered from symptomatic potato tubers with an unknown causal agent. Samples were pair-end sequneced at 150 bp with a depth of 25 M reads with the NextSeq2000. 


