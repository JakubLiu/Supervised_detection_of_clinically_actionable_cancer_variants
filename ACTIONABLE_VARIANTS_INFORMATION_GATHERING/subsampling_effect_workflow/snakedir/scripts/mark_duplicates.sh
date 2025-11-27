#!/usr/bin/bash


# This script does not remove duplicates, it just makrs them.

input_bam="$1"
output_bam="$2"
sample_name="$3"

samtools sort -n -o ${sample_name}.namesort.bam $input_bam   # namesort

samtools fixmate -m ${sample_name}.namesort.bam ${sample_name}.fixmate.bam   # fixmate 

samtools sort -o ${sample_name}.positionsort.bam ${sample_name}.fixmate.bam   # (normal) position sort

samtools markdup ${sample_name}.positionsort.bam $output_bam    # mark and remove (-r) duplicates

samtools index $output_bam    # index

# remove intermediate files
rm ${sample_name}.namesort.bam
rm ${sample_name}.fixmate.bam
rm ${sample_name}.positionsort.bam
