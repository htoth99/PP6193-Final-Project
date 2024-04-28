#!/bin/bash
#SBATCH --job-name=metaspades
#SBATCH -A PAS2700
#SBATCH --time=20:00:00
#SBATCH --mem=150GB
#SBATCH --out=metaspades-%j.out

set -euo pipefail

# Set start time
start_time=$(date)
echo "Starting time and date: $start_time"

# Load python
module load python/3.6-conda5.2

# Error if incorrect number of arguments were given
if [[ ! "$#" -eq 2 ]]; then
    echo "Error: You provided $# arguments, while 2 are required."
    echo "Usage: python spades.py -1 R1 -2 R2 -o <output_dir>"
    echo "Example: sbatch metaspades.sh results/trim/sample_R1.fastq -o results/metaspades"
    echo "Your arguments: $*"
    exit 1
fi

# Load SPAdes - must place computer insides SPAdes directory
cd /fs/scratch/PAS1568/htoth99/SPAdes-3.15.5-Linux/bin

# Load bash variables
R1=$1
output_dir=$2

# Infer R2
R2=$(echo "$R1" | sed -e "s/_R1/_R2/")

# Run metaspades from home directory
python spades.py \
-1 "$R1" \
-2 "$R2" \
-o "$output_dir" --only-assembler --meta --memory 150

# When was this analysis run and with what sample?
name=$(basename "$R1" _R1.fastq)
end_time=$(date)
echo "This analysis ended at $end_time with $name sample"