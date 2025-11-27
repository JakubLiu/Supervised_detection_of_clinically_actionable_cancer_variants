#!/usr/bin/bash

input="$1"
output="$2"

bcftools sort -Oz $input -o $output
bcftools index $output
tabix -p vcf $output
