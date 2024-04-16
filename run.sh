# shellcheck disable=SC2148
# This is a runner script for my final project in PP6193
# It is designed to be a workflow for an assembly-first MGS approach
# This workflow contains 5 steps - please refer to the README.md and USAGE.md for more details
# Copy and paste each section into the terminal to run each step

sample_name=PDMG-52


# First, we need to run fastqc to test quality by running fastqc

# Define fastqc directories
input_fastqc=rawdata
output_fastqc=results/fastqc

# Step 1: Run fastqc in a loop
for file in "$input_fastqc"/*.fastq; do
	sbatch scripts/fastqc.sh "$file" "$output_fastqc"
done

# Next, we trim our reads using trimmomatic
# note - I'm reusing my variable above as my input
output_trim=results/trim

# Step 2: Run trimmommatic in a loop
for R1_t in "$input_fastqc"/*fastq; do
	sbatch scripts/trimmomatic.sh "$R1_t" "$output_trim"
done

# After trimming, we need to rename the files using bash syntax

# Step 3: Assembling the genome using metaspades
output_metaspades=results/metaspades

for R1_a in "output_trim"/*.fastq; do
	sbatch scripts/metaspades.sh "$R1_a" "$output_metaspades"
done

# After assembling genome, rename the "contigs.fasta" to the actual name
mv results/metaspades/contigs.fasta results/metaspades/"$sample_name".fasta

# Step 4: Read mapping - lots of steps to this one!
readmap_dir=results/read_map

sbatch scripts/bwa.sh "$sample_name" "<enter R1 here>" "$readmap_dir"/aln_sam "$readmap_dir"/aln_bam "$readmap_dir"/sorted_bam

# Step 5: Binning

