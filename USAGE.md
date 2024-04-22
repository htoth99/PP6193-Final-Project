# Metagenomic workflow usage guide
### Created by Hannah Toth

## What does this workflow include?
There are **five** steps detailed in this workflow. They are:
- quality assessment (fastqc)
- trimming of raw reads (trimmomatic)
- assembly of metagenome (metaSPAdes)
- read mapping (bwa)
- binning (maxbin2/TBD)

## How should I set up my project directory?
Begin with a root directory. Within this directory, the following files and directories will need to be created prior to running any step, where bold denotes a directory, an italicized denotes a file
* **rawdata**
   * *samplename1_R1.fastq*
   * *samplename1_R2.fastq*
   * *samplename2_R1.fastq*
   * *samplename2_R2.fastq*
* **scripts**
   * **logs**
* **results**
   * **fastqc**
   * **trim**
   * **metaspades**
   * **read_map**
   * **binning**
*   *USAGE.md*
*   *run.sh*
<br>
<br>
It's important to set up the project like this, as the runner script uses relative paths. **The runner script will not work if the project is not set up exactly like this.**

## How do I run each step?
Each step has it's own script in the script directory of this project. However, they are run by copying and pasting the corresponding block from the **runner script**, which is located in the home directory. **Step 1 must come before step 2, and so on.**
<br>

Each step is denoted like this:
```bash
### Step 1: Example of step 1 ###
```
Actions that should be copied and pasted together are grouped together in paragraphs. Each step follows the same formula:
1. Define directories (input, output, etc.)
2. Run the for loop to submit the job(s)
3. If applicable, rename, move, or delete files
4. Move log file to appropriate location

## Do I need to edit any of the scripts?
Yes, you will. First and foremost, you will need to edit each script for your funding source if running on the OSC.
```bash
#SBATCH -A <funding project>
```
Next, you may need to edit the time and memory for each script depending on your data size. The runner script is set up to run one job per sample (with the exception of fastqc and trimmomatic, where a job for each file is run).
```bash
#SBATCH --mem=<memoryGB or MB>
#SBATCH --time=<00:00:00>
```
### Specific cases
Trimmomatic has two different versions on the OSC. 0.36 is only on the Owens cluster, and 0.38 is only on the pitzer cluster. Be sure to change this before running on the appropriate cluster
```bash
module load trimmomatic/0.36
module load trimmomatic/0.38
```
An envrionment is required for metabat2. My envrionment for metabat2 is named "metabat2" and I activate my envrionment through mamba, as seen in the script. It can alternatively be activated through conda. To install, use the following commands:
```bash
conda create -n metabat2
conda activate metabat2
conda install bioconda::metabat2
```
SPAdes is downloaded to a scratch folder of mine. I recommend you put SPAdes in a scratch folder as well, as it can be reinstalled easily, and is prone to core dump if it is having issues.
To install, navigate to the directory where you would like to install SPAdes and do the following:
```bash
wget https://github.com/ablab/spades/releases/download/v3.15.5/SPAdes-3.15.5-Linux.tar.gz
tar -xzf SPAdes-3.15.5-Linux.tar.gz
cd SPAdes-3.15.5-Linux/bin/
```
After installation, change the path in the metaspades.sh script to where your SPAdes bin directory is located.

## What about the runner script?
Nope! It should be fine how it is. Just ensure that your runner script is in the root directory and the rest of the project is set up as decsrcibed above.
