#!/usr/bin/bash

input_vcf="$1"
output_vcf="$2"
bcftools view -f PASS "$input_vcf" -Oz -o "$output_vcf"
tabix -p vcf "$output_vcf"