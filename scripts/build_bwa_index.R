if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Rbwa")
library(Rbwa)

dir <- "../data/genomes/GRCh38/"

# Building the index for Homo Sapiens
fasta <- "../data/genomes/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

index_prefix <- file.path(dir, "GCA_000001405.15_GRCh38_no_alt_analysis_set")

bwa_build_index(fasta,
                index_prefix=index_prefix)

