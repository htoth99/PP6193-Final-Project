#!/bin/bash
#SBATCH -A PAS2700
#SBATCH --job-name=fastqc
#SBATCH --out=fastqc-%j.out

set -euo pipefail

# Load fastqc
module load fastqc

# Bash variables
fastq_file=$1
output_dir=$2

# Run fastqc: fastqc <input> -o <output>
fastqc "$fastq_file" -o "$output_dir"

