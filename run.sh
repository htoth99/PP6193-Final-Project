# shellcheck disable=SC2148
# This is a runner script for my final project in PP6193
# It is designed to be a workflow for an assembly-first MGS approach
# This workflow contains 5 steps - please refer to the README.md and USAGE.md for more details
# Copy and paste each section into the terminal to run each step

sample_name=PDMG-52

### Step 1: Run fastqc in a loop
# Define fastqc directories
input_fastqc=rawdata
output_fastqc=results/fastqc

# Run fastqc
for file in "$input_fastqc"/*.fastq; do
	sbatch scripts/fastqc.sh "$file" "$output_fastqc"
done
# Before running next step, move log files to log directory
mv fastqc-* scripts/logs

### Step 2: Run trimmommatic in a loop
# Define trimmomatic directories
output_trim=results/trim

# Run trimmomatic
for R1_t in "$input_fastqc"/*fastq; do
	sbatch scripts/trimmomatic.sh "$R1_t" "$output_trim"
done

# Before running next step, move log files to log directory
mv trimmomatic-* scripts/logs

# After trimming, we need to...
# Delete unpaired files
rm results/trim/*_unpaired1.fastq
rm results/trim/*_unpaired2.fastq
# Rename paired files to only contain R1 and R2 at the end
for f in results/trim/*.fastq;do mv "$f" "${f/_paired1./.}";done
for f in results/trim/*.fastq;do mv "$f" "${f/_paired2./.}";done

### Step 3: Assembling the genome using metaspades
# Metaspades directory
output_metaspades=results/metaspades

for R1_a in "$output_trim"/*.fastq; do
	sbatch scripts/metaspades.sh "$R1_a" "$output_metaspades"
done
# Before running next step, move log files to log directory
mv metaspades-* scripts/logs

# After assembling genome, rename the "contigs.fasta" to the actual name
mv results/metaspades/contigs.fasta results/metaspades/"$sample_name".fasta

### Step 4: Read mapping - lots of steps to this one!
readmap_dir=results/read_map

sbatch scripts/bwa.sh "$sample_name" "<enter R1 here>" "$readmap_dir"/aln_sam "$readmap_dir"/aln_bam "$readmap_dir"/sorted_bam
# Before running next step, move log files to log directory
mv bwa-* scripts/logs

### Step 5: Binning

