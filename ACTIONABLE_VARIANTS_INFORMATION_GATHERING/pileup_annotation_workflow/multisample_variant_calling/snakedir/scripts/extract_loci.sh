#!/usr/bin/bash

multi_vcf="$1"
output="$2"

awk -F ',' '{print $1}' $multi_vcf > $output