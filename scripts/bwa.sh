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
R1_in=$2
bwa_outdir=$3

# Infer R2
R2_in=$(echo $R1_in | sed -e "s/_R1/_R2/")

# File basename
N=$(basename "$R1_in" .fastq)

# First, index the contigs from the assembly
bwa index "$contigs"

# Next, align the pair-end sequences back to the indexed genome
bwa mem -t 16 -a "$contigs" \
"$R1_in" "$R2_in" > "$bwa_outdir"/"$N"_aln_sam.sam

# Make a bam file from the alignment using samtools
samtools view -b -S "$bwa_outdir"/"$N"_aln_sam.sam > "$bwa_outdir"/"$N"_aln_bam.bam

# Lastly, sort the bam file since many binners need sorted bam files
samtools sort -o "$bwa_outdir"/"$N"_aln_bam.bam "$bwa_outdir"/"$N"_sorted_bam.bam

# Remove files - aligning creates massive files that aren't needed
rm "$bwa_outdir"/"$N"_aln_sam.sam "$bwa_outdir"/"$N"_aln_bam.bam


