# Metagenomic workflow usage guide
### Created by Hannah Toth

## How is this workflow set up?
There are **five** steps detailed in this workflow. They are:
- quality assessment (fastqc)
- trimming of raw reads (trimmomatic)
- assembly of metagenome (metaSPAdes)
- read mapping (bwa)
- binning (maxbin2/TBD)

Each step has it's own script in the script directory of this project. However, they are run by copying and pasting the corresponding block from the **runner script**, which is located in the home directory.