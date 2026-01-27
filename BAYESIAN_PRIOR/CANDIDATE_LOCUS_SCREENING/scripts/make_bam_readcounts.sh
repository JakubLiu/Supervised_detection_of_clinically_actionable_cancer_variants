#!/usr/bin/bash


reference_genome="$1"
bedfile="$2"
input_bamfile="$3"
output_file="$4"

bam-readcount -f $reference_genome -l $bedfile $input_bamfile > $output_file

