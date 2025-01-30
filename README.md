# star_make
Snakemake pipeline for bulk RNAseq alignment with STAR

## Prerequisites
- This pipeline assumes the genome has been properly built using STAR.
- This pipeline also assumes that the FASTQ files are present in a directory called "samples" that is within the snakemake pipeline directory. The relative paths are setup to deal with input files in this manner.


## Input files 
We need a tab-delimited input file called STAR_input_samples.txt. There should be 3 columns in this file: 
- sample_name
- r1_file
- r2_file

This is read into the snakemake pipeline, then a separate STAR job is started for each input sample. This pipeline handles paired-end reads into STAR appropriately. 

## Testing setup
We can test the setup with: 
snakemake --dry-run 

## Running the pipeline
Run the pipeline with:
snakemake --profile=slurm_htc --jobs={number_of_files_to_align + 1}

Note the +1 is to capture the "all" feature of snakemake.
