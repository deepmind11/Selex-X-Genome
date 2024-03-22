#!/usr/bin/env Rscript

# Lenght 1 specfying the full path to the input file
args = commandArgs(trailingOnly=TRUE)

outputdir = get_parent_dir(args[1])
filename = basename(args[1])


# Defining the index (human for now)
dir <- "../data/genomes/GRCh38/"
index_prefix <- file.path(dir, "GCA_000001405.15_GRCh38_no_alt_analysis_set")


if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Rbwa")
library(Rbwa)

# BWA align <70bp
fastq <- args[1]

# strtrim(readLines(file.path(fastq)), 100)

bwa_aln(index_prefix=index_prefix,
        fastq_files=fastq,
        sai_files=file.path(outputdir, "output.sai"))

# Create SAM file
bwa_sam(index_prefix=index_prefix,
        fastq_files=fastq,
        sai_files=file.path(outputdir, paste(filename[1:12], "sai", sep = '.')),
        sam_file=file.path(outputdir, paste(filename[1:12], "sam", sep = '.')))








