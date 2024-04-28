# Metagenomic workflow usage guide
### Created by Hannah Toth
This workflow is for use on the Ohio Supercomputer. Please read the entire document throughly before attempting any steps, as you need to edit some items, which is described below.
## What does this workflow include?
There are **five** steps detailed in this workflow. They are:
- quality assessment (fastqc)
- trimming of raw reads (trimmomatic)
- assembly of metagenome (metaSPAdes)
- read mapping (bwa)
- binning (maxbin2 or metabat2)

## How should I set up my project directory?
Begin with a root directory. Within this directory, the following files and directories will need to be created prior to running any step, where bold denotes a directory, an italicized denotes a file
* **rawdata**
   * *samplename1_R1.fastq*
   * *samplename1_R2.fastq*
   * *samplename2_R1.fastq*
   * *samplename2_R2.fastq*
* **scripts**
   * *fastqc.sh*
   * *trimmomatic.sh*
   * *metaspades.sh*
   * *bwa.sh*
   * *maxbin2.sh*
   * *metabat2.sh*
* **results**
   * **fastqc**
   * **trim**
   * **metaspades**
   * **read_map**
   * **binning**
   * **logs**
*   *USAGE.md*
*   *run.sh*
<br>
<br>
It's important to set up your project directory exactly as described above, as the runner script uses relative paths. If downloading from github, then the results and rawdata directories and subdirectories must be made manually. This is detailed in the runner script but commented out. Fill the rawdata directory with your raw, paired-end reads. **The runner script will not work if the project is not set up exactly like this.**

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

## Editing scripts and installation of software
This workflow was created so users do not need to create or edit scripts. There are a few exceptions to this concept, and they are described below. *These changes are needed in order for the workflow to work*

## What do I need to install?
An envrionment is required for metabat2 and maxbin2. My envrionment for metabat2 is named "metabat2" and I activate my envrionment through mamba, as seen in the script. My maxbin2 envrionment follows the same format and activation as metabat2. They can alternatively be activated through conda rather than mamba.
<br>
**Metbat2 installation**
```bash
conda create -n metabat2
conda activate metabat2
conda install bioconda::metabat2
```
**Maxbin2 installation**
```bash
conda create -n maxbin2
conda activate maxbin2
conda install bioconda::maxbin2
```
SPAdes is downloaded to a scratch folder of mine. I recommend you put SPAdes in a scratch folder as well, as it can be reinstalled easily, and is prone to core dump if it is having issues.
To install, navigate to the directory where you would like to install SPAdes and do the following:
```bash
wget https://github.com/ablab/spades/releases/download/v3.15.5/SPAdes-3.15.5-Linux.tar.gz
tar -xzf SPAdes-3.15.5-Linux.tar.gz
cd SPAdes-3.15.5-Linux/bin/
```
After installation, change the path in the metaspades.sh script to where your SPAdes bin directory is located.

## What do I need to edit?
### SBATCH lines in each script
You will need to edit the **funding source** and the **time** and **mem** for each script just one time. You may need to adjust the time and memory for smaller or larger projects in the future. Please change the following lines in each script:
```bash
#SBATCH -A <funding project>
#SBATCH --mem=<memoryGB or MB>
#SBATCH --time=<00:00:00>
```
### trimmomatic.sh
Trimmomatic has two different versions on the OSC. 0.36 is only on the Owens cluster, and 0.38 is only on the pitzer cluster. Be sure to change this in the *trimmomatic.sh* script located in scripts/ before running on the appropriate cluster.
```bash
module load trimmomatic/0.36
module load trimmomatic/0.38
```
### metaspades.sh
After spades installation, change the path in the *metaspades.sh* script, located in scripts/ and change line 27 to your spades location. 

# FAQ

## There's 2 binning programs, do I need to run both?
No, you do not need to run both. You just need to run one. I like both and I have not decided which one I will continue using.

## What will the logs tell me?
Logs are named by program-batchjob.out and are located in results/logs. In them, you will find what the verbose from the program. You will also find the start time, end time, and what sample was run as well. 

## What doesn't this workflow include?
This workflow does not (currently) include host removal or checkm. Both of these components are important for plant pathogenic metagenomic samples. I plan to add both of these componets in the near future.

This workflow also does not include extensive error messaging/guiding. This is something I will continue to work on once I've established and utilized this pipeline across many samples. I also want to have others test it this out so I can see the error others may be encountering. Since the set up currently works for me, and is based on copying and pasting from the runner scripts, it was challenging for me to generate errors!

## What if I'm getting an error or can't get something to work?
We all run into errors! Check the log files to find out why your job didn't work. Here are some common steps to take to solving your errors.
* Do you have enough memory and/or time allocated?
* Does your directory and/or file exist? Is it the correct path?
* Do your envrionment exist? Is it activated?
* Did you copy and paste all of the directories into the terminal before submitting the job?

If all else fails, and you continue running into errors, please contact me at toth.302@buckeyemail.osu.edu