#!/bin/bash
#SBATCH -A PAS2700
#SBATCH --job-name=bwa
#SBATCH --mem=5GB
#SBATCH --time=06:00:00
#SBATCH --output=bwa.out

set -euo pipefail

# Load modules for read mapping
module load bwa
module load samtools/1.16.1

# Set bash variables
contigs=$1
read1=$2
aln_sam=$3
aln_bam=$4
sorted_bam=$5

# Infer R2
read2=${R1/_R1/_R2}

# First, indext the contigs from the assembly
bwa index "$contigs"

# Next, align the pair-end sequences back to the indexed genome
bwa mem -t 16 -a "$contigs" \
"$read1" "$read2" > "$aln_sam"

# Make a bam file from the alignment using samtools
samtools view -b -S "$aln_sam" > "$aln_bam"

# Lastly, sort the bam file since many binners need sorted bam files
samtools sort -o "$aln_bam" "$sorted_bam"
