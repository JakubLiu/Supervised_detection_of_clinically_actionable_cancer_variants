#!/usr/bin/bash

vcf_in="$1"
annotation_file="$2"
annotated_vcf_out="$3"
snpEff_dir="$4"

java -Xmx4g -jar "$snpEff_dir/snpEff.jar" $annotation_file $vcf_in > $annotated_vcf_out