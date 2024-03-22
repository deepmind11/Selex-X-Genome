#!/bin/bash




BWT_INDEX=/burg/home/hg2604/Remap/Human/GRCh38_noalt_BWT_index/GRCh38_noalt_as
REFERENCE=/burg/hblab/users/hg2604/ProboundMetrics/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna


# Path to subsampled fastq file
read1=$1
read2=$2

currentdr=$(dirname "$read1")
filename=$(basename "$read1")
stemname=${filename:0:11}



#By default, bowtie2 filters fragments greater than 500bp
bowtie2 -x $BWT_INDEX -1 "$read1" -2 "$read2" --threads 1  --no-unal | samtools view -bS - >  "$currentdr"/"$stemname".bam

#Sorting the bam file
samtools sort "$currentdr"/"$stemname".bam -o "$currentdr"/"$stemname"_sorted.bam
