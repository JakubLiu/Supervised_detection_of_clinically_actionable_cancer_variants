#!/usr/bin/bash

tumor_bam="$1"
normal_bam="$2"
normal_sample_name="$3"
reference_genome="$4"
output_vcf="$5"
threads="${6:-8}"  # default to 8 threads if not provided


gatk Mutect2 \
    -R "$reference_genome" \
    -I "$tumor_bam" \
    -I "$normal_bam" \
    -normal "$normal_sample_name" \
    -O "$output_vcf"
