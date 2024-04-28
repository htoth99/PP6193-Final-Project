#!/bin/bash
#SBATCH -A PAS2700
#SBATCH --job-name=fastqc
#SBATCH --out=fastqc-%j.out

set -euo pipefail

# Load fastqc
module load fastqc

# Set start time
start_time=$(date)
echo "Starting time and date: $start_time"

# Error if incorrect number of arguments were given
if [[ ! "$#" -eq 2 ]]; then
    echo "Error: You provided $# arguments, while 2 are required."
    echo "Usage: fastqc.sh <FASTQ-file> -o <output-dir>"
    echo "Example: sbatch fastqc.sh rawdata/sample_R1.fastq -o results/fastqc"
    echo "Your arguments: $*"
    exit 1
fi

# Bash variables
fastq_file=$1
output_dir=$2

# Run fastqc: fastqc <input> -o <output>
fastqc "$fastq_file" -o "$output_dir"

# When was this analysis run and with what sample?
name=$(basename "$fastq_file" .fastq)
end_time=$(date)
echo "This analysis ended at $end_time with $name"
