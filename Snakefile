# Read in data
import pandas as pd
samples_table = pd.read_csv("STAR_input_files.txt",sep="\t").set_index("sample_name",drop=False)

# Definitions for rule all
SAMPLES_ALIGN = samples_table.loc[:,"sample_name"].values.tolist()

# Add samples prefix to all fastq files to specify directory
samples_table["r1_file"] = samples_table["r1_file"].apply(lambda x: 'samples/' + x)
samples_table["r2_file"] = samples_table["r2_file"].apply(lambda x: 'samples/' + x)

# Define local rule
localrules:
	all

# Rule to define output directories
rule all:
	input:
		expand("star_outs/{final_aligned}",final_aligned=SAMPLES_ALIGN)

# Function to read in gex fastq paths from samples_table
def get_align_input(wildcards):
	return samples_table.loc[wildcards.sample][["r1_file","r2_file"]].to_list()

# STAR align and count rule
rule count:
	input:
		get_align_input
	output:
		directory("star_outs/{sample}/")
	threads:
		8
	resources:
		runtime="1d",
		mem_mb=128000
	shell:
		"""
		module load star/2.7.11b
		"""
		"""
		mkdir -p star_outs
		mkdir -p star_outs/{wildcards.sample}
		STAR --runThreadN 8 \
			--genomeDir /ix1/acillo/arc85/references/star_murine_ref_250129/murine_star_indices \
			--readFilesIn {input} \
			--readFilesCommand zcat \
			--outFileNamePrefix star_outs/{wildcards.sample}/{wildcards.sample} \
			--outSAMtype BAM Unsorted \
			--quantMode GeneCounts
		"""