#!/usr/bin/bash

input_bam="$1"
subsampling_proportion_file="$2"
subsampling_random_seed="$3"
output_subsamled_bam="$4"

subsampling_proportion=$(awk -F "," '{print $3}' "$subsampling_proportion_file" | tail -n 1)

samtools view --subsample "$subsampling_proportion" --subsample-seed "$subsampling_random_seed" -b "$input_bam" > "$output_subsamled_bam"