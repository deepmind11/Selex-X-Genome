#!/bin/bash

BWT_INDEX=/burg/hblab/users/hg2604/Remap/Human/GRCh38_noalt_BWT_index/GRCh38_noalt_as
REFERENCE=/burg/hblab/users/hg2604/Projects/Selex-X-Genome/data/genomes/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

# Full path
fastq_file=$1
currentdr=$(dirname "$fastq_file")
filename=$(basename "$fastq_file")
stemname=${filename:0:11}

#Mapping the processed fastq files
bowtie2 -x $BWT_INDEX -U "$fastq_file" --threads 1  --no-unal | samtools view -bS - >  "$currentdr"/"$stemname".bam

#Sorting the bam file
samtools sort "$currentdr"/"$stemname".bam -o "$currentdr"/"$stemname"_sorted.bam

#Getting the bed file
bedtools bamtobed -i "$currentdr"/"$stemname"_sorted.bam > "$currentdr"/"$stemname".bed

#Processing the bed file
awk 'BEGIN {OFS="\t"} {original2 = $2; original3 = $3; $2 = ($6 == "+") ? original2 - 0 : original3 - 200; $3 = ($6 == "+") ? original2 + 200 : original3 + 0; print}' "$currentdr"/"$stemname".bed > "$currentdr"/temp
awk '$5 > 20 && $2 > 0 && $3 > 0'  "$currentdr"/temp >  "$currentdr"/pr"$stemname".bed
rm "$currentdr"/temp

# Bed to Seq
bedtools getfasta -fi $REFERENCE -bed "$currentdr"/pr"$stemname".bed -s -fo "$currentdr"/"$stemname".fasta

# Deleting all the temporary files
rm "$currentdr"/"$stemname".bed
rm "$currentdr"/pr"$stemname".bed
rm "$currentdr"/"$stemname".bam
rm "$currentdr"/"$stemname"_sorted.bam
rm "$currentdr"/"$stemname".fastq.gz
#rm "$fastq_file"
