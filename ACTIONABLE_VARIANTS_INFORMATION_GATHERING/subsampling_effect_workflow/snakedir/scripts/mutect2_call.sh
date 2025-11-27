#!/usr/bin/bash

tumor_bam="$1"
normal_bam="$2"
normal_sample_name="$3"
reference_genome="$4"
chromosome="$5"
output_vcf="$6"


gatk Mutect2 \
    -R "$reference_genome" \
    -I "$tumor_bam" \
    -I "$normal_bam" \
    -normal "$normal_sample_name" \
    --intervals "$chromosome" \
    -O "$output_vcf"

