#!/bin/bash
#SBATCH --job-name=metabat2
#SBATCH -A PAS2700
#SBATCH --time=6:00:00
#SBATCH --mem=60GB
#SBATCH --out=metabat2-%j.out

set -euo pipefail

# Set start time
start_time=$(date)
echo "Starting time and date: $start_time"

# Activate metabat2 environment
mamba activate metabat2

# Error if incorrect number of arguments were given
if [[ ! "$#" -eq 2 ]]; then
    echo "Error: You provided $# arguments, while 2 are required."
    echo "Usage: python spades.py -1 R1 -2 R2 -o <output_dir>"
    echo "Example: sbatch metaspades.sh results/trim/sample_R1.fastq -o results/metaspades"
    echo "Your arguments: $*"
    exit 1
fi

# Define bash variables
contigs=$1
output=$2

# File name
N=$(basename "$contigs" .fasta)

# Run metabat2
metabat2 \
-i "$contigs" \
-o "$output"

# When was this analysis run and with what sample?
end_time=$(date)
echo "This analysis ended at $end_time with $N"
