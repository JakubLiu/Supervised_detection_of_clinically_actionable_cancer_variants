#!/usr/bin/bash

# note that this script just adds the filtration flags to the vcf, it does not remove any variant calls

reference_genome="$1"
input_vcf="$2"
output_vcf="$3"

gatk FilterMutectCalls -R $reference_genome -V $input_vcf -O $output_vcf