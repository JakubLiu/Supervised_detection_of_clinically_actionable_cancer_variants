#!/usr/bin/bash

input_bam="$1"
output="$2"
regions_bed="$3"

samtools depth -b $regions_bed $input_bam -o $output

# usage in the snakefile ./samtools_depth {input} {output} {params.regions} 