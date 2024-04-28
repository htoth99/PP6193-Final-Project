# This is a runner script for my final project in PP6193
# Created by Hannah Toth
# It is designed to be a workflow for an assembly-first MGS approach
# This workflow contains 5 steps - please refer to the README.md and USAGE.md for more details
# Copy and paste each section into the terminal to run each step

### Step 1: Run fastqc in a loop ###
# Define fastqc directories
input_fastqc=rawdata
output_fastqc=results/fastqc

# Run fastqc
for file in "$input_fastqc"/*.fastq; do
	sbatch scripts/fastqc.sh "$file" "$output_fastqc"
done
# Before running next step, move log files to log directory
mv fastqc-* results/logs

### Step 2: Run trimmommatic in a loop ###
# Define trimmomatic directories
input_fastqc=rawdata
output_trim=results/trim

# Run trimmomatic
for R1_t in "$input_fastqc"/*_R1.fastq; do
	sbatch scripts/trimmomatic.sh "$R1_t" "$output_trim"
done

# Before running next step, move log files to log directory
mv trimmomatic-* results/logs

# After trimming, we need to...
# Delete unpaired files
rm results/trim/*_unpaired1.fastq
rm results/trim/*_unpaired2.fastq
# Rename paired files to only contain R1 and R2 at the end
for f in results/trim/*.fastq;do mv "$f" "${f/_paired1./.}";done
for f in results/trim/*.fastq;do mv "$f" "${f/_paired2./.}";done

### Step 3: Assembling the genome using metaspades ###
# Metaspades directories
# Work around coding because SPAdes needs to run in its own directory!
root_dir=$(pwd)
input_metaspades="$root_dir"/results/trim
output_metaspades="$root_dir"/results/metaspades

for R1_a in "$input_metaspades"/*_R1.fastq; do
	file_name=$(basename "$R1_a" .fastq)
	mkdir -p "$output_metaspades"/"$file_name"
	sbatch scripts/metaspades.sh "$R1_a" "$output_metaspades"/"$file_name"
done

# SPAdes will name generically, so we should rename and move assembleed contigs to another directory
cd results/metaspades
mkdir fastas
for f in **/contigs.fasta; do newf=$(dirname "$f"); mv $f fastas/$newf.fasta;done
cd ../..

# Before running next step, move log files to log directory
mv metaspades-* results/logs

### Step 4: Read mapping - lots of steps to this one! ###
# Read mapping is only required for metabat2
# Define read mapping directories
readmap_dir=results/read_map
assemblies_dir=results/metaspades/fastas
output_trim=results/trim

for file in "$assemblies_dir"/*.fasta; do
	name=$(basename "$file" .fasta)
	sbatch scripts/bwa.sh "$file" "$output_trim"/"$name"_R1.fastq "$readmap_dir"
done

# Before running next step, move log files to log directory
mv bwa-* results/logs

### Step 5: Binning ###
# When using metabat2, we need a)contigs and b)sorted bam file
# Define metabat2 and maxbin2 directories
readmap_dir=results/read_map
assemblies_dir=results/metaspades/fastas
metabat_out=results/binning/metabat2
output_trim=results/trim
maxbin_out=results/binning/maxbin2/

#jgi_summarize_bam_contig_depths -outputDepth "$output"/"$N"_depth.txt $bam_file

# Run metabat2
for file in "$assemblies_dir"/*.fasta; do
	name=$(basename "$file" .fasta)
	sbatch scripts/metabat2.sh "$file" "$metabat_out"/"$name"
done

# Or run maxbin2
for file in "$assemblies_dir"/*.fasta; do
	name=$(basename "$file" .fasta)
	sbatch scripts/maxbin2.sh "$file" "$output_trim"/"$name"_R1.fastq "$maxbin_out"/"$name"
done

# Move log file to log directory
mv metabat2-* results/logs

