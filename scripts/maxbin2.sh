#!/bin/bash
#SBATCH -A PAS2700
#SBATCH --job-name=maxbin2
#SBATCH --cpus-per-task=20
#SBATCH --time=20:00:00
#SBATCH --output=maxbin2-%j.out

set -euo pipefail

# Set start time
start_time=$(date)
echo "Starting time and date: $start_time"

# Activate maxbin2 envrionment
mamba activate maxbin2

# Error if incorrect number of arguments were given
if [[ ! "$#" -eq 2 ]]; then
    echo "Error: You provided $# arguments, while 3 are required."
    echo "Usage: run_MaxBin.pl -contigs <contigs> -reads <R1> -reads2 <R2> -out <outdir>"
    echo "Example: sbatch maxbin2.sh results/metaspades/sample.fasta results/trim/sample_R1.fastq -o results/binning"
    echo "Your arguments: $*"
    exit 1
fi

# Define variables
contigs=$1
R1_in=$2
outdir=$3

# Infer R2 from R1
R2_in=$(echo "$R1_in" | sed -e "s/_R1/_R2/")

# Run Maxbin2
run_MaxBin.pl \
-contig "$contigs" \
-reads "$R1_in" \
-reads2 "$R2_in" \
-out "$outdir"

# When was this analysis run and with what sample?
name=$(basename "$contigs" .fasta)
end_time=$(date)
echo "This analysis ended at $end_time with $name"