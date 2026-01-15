#!/usr/bin/bash

input_bam="$1"
output_bam="$2"
threads="$3"

samtools sort -@ "$threads" -o "$output_bam" "$input_bam"

samtools index "$output_bam"