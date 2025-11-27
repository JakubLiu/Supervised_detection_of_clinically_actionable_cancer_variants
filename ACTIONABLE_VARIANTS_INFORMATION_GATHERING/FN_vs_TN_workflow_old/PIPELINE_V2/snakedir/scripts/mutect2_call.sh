#!/usr/bin/bash

tumor_bam="$1"
normal_bam="$2"
normal_sample_name="$3"
reference_genome="$4"
output_vcf="$5"

vcf_dir=$(dirname "$output_vcf")

sample=$(basename "$output_vcf" | cut -d'_' -f2)

# Run Mutect2
gatk Mutect2 \
    -R "$reference_genome" \
    -I "$tumor_bam" \
    -I "$normal_bam" \
    -normal "$normal_sample_name" \
    -O "$output_vcf"

