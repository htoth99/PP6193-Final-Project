#!/bin/bash
#SBATCH -A PAS2700
#SBATCH --job-name=trimmomatic
#SBATCH --mem=10GB
#SBATCH --time=05:00:00
#SBATCH --output=trimmomatic-%j.out

# Load trimmomatic module (0.36 for owens, 0.38 for pitzer)
#module load trimmomatic/0.36
module load trimmomatic/0.38

# Load bash variables
R1=$1
output_dir=$2

# Find R2 from R1
R2=${R1/_R1/_R2}

# Define basename
n1=$(basename "$R1" .fastq)
n2=$(basename "$R2" .fastq)

# Run trimmomatic
java -jar $TRIMMOMATIC PE -phred33 \
"$R1" "$R2" \
"$output_dir"/"$n1"_paired1.fastq \
"$output_dir"/"$n1"_unpaired1.fastq \
"$output_dir"/"$n2"_paired2.fastq \
"$output_dir"/"$n2"_unpaired2.fastq \
ILLUMINACLIP:/fs/ess/PAS1568/Taylor/BacterialSpot/NexteraPE-PE.fa:2:30:10 LEADING:10 TRAILING:2 SLIDINGWINDOW:4:15 MINLEN:36
