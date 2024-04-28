#!/bin/bash
#SBATCH -A PAS2700
#SBATCH --job-name=bwa
#SBATCH --mem=20GB
#SBATCH --time=06:00:00
#SBATCH --output=bwa-%j.out

set -euo pipefail

# Set start time
start_time=$(date)
echo "Starting time and date: $start_time"

# Load modules for read mapping
module load bwa
module load samtools/1.16.1

# Error if incorrect number of arguments were given
if [[ ! "$#" -eq 3 ]]; then
    echo "Error: You provided $# arguments, while 3 are required."
    echo "Please provide the assembly, trimmed R1.fastq file, and an output directory"
    echo "Example: sbatch metaspades.sh results/trim/sample_R1.fastq -o results/metaspades"
    echo "Your arguments: $*"
    exit 1
fi

# Set bash variables
contigs=$1
R1_in=$2
bwa_outdir=$3

# Infer R2
R2_in=$(echo "$R1_in" | sed -e "s/_R1/_R2/")

# File basename
N=$(basename "$contigs" .fasta)

# First, index the contigs from the assembly
bwa index "$contigs"

# Next, align the pair-end sequences back to the indexed genome
bwa mem -t 16 -a "$contigs" \
"$R1_in" "$R2_in" > "$bwa_outdir"/"$N"_aln_sam.sam

# Make a bam file from the alignment using samtools
samtools view -b -S "$bwa_outdir"/"$N"_aln_sam.sam > "$bwa_outdir"/"$N"_aln_bam.bam

# Lastly, sort the bam file since many binners need sorted bam files
samtools sort -o "$bwa_outdir"/"$N"_sorted_bam.bam "$bwa_outdir"/"$N"_aln_bam.bam

# Remove files - aligning creates massive files that aren't needed
rm "$bwa_outdir"/"$N"_aln_sam.sam "$bwa_outdir"/"$N"_aln_bam.bam

# When was this analysis run and with what sample?
end_time=$(date)
echo "This analysis ended at $end_time with $N"
