#!/usr/bin/bash

bamlist_file="$1"
bedfile="$2"
reference_genome="$3"
output="$4"

samtools mpileup -l $bedfile -f $reference_genome -b $bamlist_file -s > $output