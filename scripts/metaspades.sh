#!/bin/bash
#SBATCH --job-name=metaspades
#SBATCH -A PAS2700
#SBATCH --time=20:00:00
#SBATCH --mem=40GB
#SBATCH --out=metaspades-%j.out

set -euo pipefail

# Load python
module load python/3.6-conda5.2

# Run SPAdes from scratch direcotry, for safety/memory purposes
cd /fs/scratch/PAS1568/htoth99/SPAdes-3.15.5-Linux/bin

# Load bash variables
R1=$1
output_dir=$2

# Infer R2
R2=$(echo "$R1" | sed -e "s/_R1/_R2/")

# Run metaspades
python spades.py \
-1 "$R1" \
-2 "$R2" \
-o "$output_dir" --only-assembler --meta