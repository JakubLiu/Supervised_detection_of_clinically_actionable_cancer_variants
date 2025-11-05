#!/usr/bin/bash

INPUT_VCF="$1"
OUTPUT_VCF="$2"

bcftools view -f PASS "$INPUT_VCF" -Oz -o "$OUTPUT_VCF"
tabix -p vcf "$OUTPUT_VCF"
