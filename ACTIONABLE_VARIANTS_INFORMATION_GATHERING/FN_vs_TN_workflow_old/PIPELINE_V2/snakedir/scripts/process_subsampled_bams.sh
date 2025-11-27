#!/usr/bin/bash


input_bam="$1"
output_sorted_bam="$2"

samtools sort -o $output_sorted_bam $input_bam

samtools index $output_sorted_bam

rm $input_bam  # remove the original unsorted subsampled bamfile