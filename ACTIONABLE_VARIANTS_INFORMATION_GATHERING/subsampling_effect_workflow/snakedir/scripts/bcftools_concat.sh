#!/usr/bin/bash

input_vcf_list="$1"
output_merged_vcf="$2"

bcftools concat -f $input_vcf_list -O z -o $output_merged_vcf