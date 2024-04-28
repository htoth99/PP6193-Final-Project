#!/bin/bash
#SBATCH -A PAS2700
#SBATCH --job-name=trimmomatic
#SBATCH --mem=10GB
#SBATCH --time=05:00:00
#SBATCH --output=trimmomatic-%j.out

set -euo pipefail

# Set start time
start_time=$(date)
echo "Starting time and date: $start_time"

# Load trimmomatic module (0.36 for owens, 0.38 for pitzer)
module load trimmomatic/0.36
#module load trimmomatic/0.38



# Error if incorrect number of arguments were given
if [[ ! "$#" -eq 2 ]]; then
    echo "Error: You provided $# arguments, while 2 are required."
    echo "Usage: java -jar $TRIMMOMATIC PE -phred33 R1 R2 paired unpaired paired unpaired"
    echo "Example: sbatch trimmomatic.sh rawdata/sample_R1.fastq results/trim"
    echo "Your arguments: $*"
    exit 1
fi

# Load bash variables
R1_in=$1
output_dir=$2

# Find R2 from R1
R2_in=$(echo $R1_in | sed -e "s/_R1/_R2/")

# Define basename
n1=$(basename "$R1_in" .fastq)
n2=$(basename "$R2_in" .fastq)

# Run trimmomatic
java -jar $TRIMMOMATIC PE -phred33 \
"$R1_in" "$R2_in" \
"$output_dir"/"$n1"_paired1.fastq \
"$output_dir"/"$n1"_unpaired1.fastq \
"$output_dir"/"$n2"_paired2.fastq \
"$output_dir"/"$n2"_unpaired2.fastq \
ILLUMINACLIP:/fs/ess/PAS1568/Taylor/BacterialSpot/NexteraPE-PE.fa:2:30:10 LEADING:10 TRAILING:2 SLIDINGWINDOW:4:15 MINLEN:36

# When was this analysis run and with what samples?
end_time=$(date)
echo "This analysis ended at $end_time with $n1 and $n2"
