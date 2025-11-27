#!/usr/bin/bash

input_R1="$1"
input_R2="$2"
output_paired_R1="$3"
output_paired_R2="$4"
output_unpaired_R1="$5"
output_unpaired_R2="$6"
threads="$7"
adapter_seq="$8"

# add more memory for Trimmomatic (because it seems to ignore the resource allocation from Snakemake)
export _JAVA_OPTIONS="-Xms4G -Xmx32G"


trimmomatic PE -threads $threads -phred33 \
  $input_R1 $input_R2 \
  $output_paired_R1 $output_unpaired_R1 \
  $output_paired_R2 $output_unpaired_R2 \
  ILLUMINACLIP:${adapter_seq}:2:30:10
