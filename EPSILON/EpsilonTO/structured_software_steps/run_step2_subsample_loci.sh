#!/usr/bin/bash

raw_bed="merged.processed.mpileup.txt"
output_subsampled_bed="subsampled.mpileup.txt"
subsampling_strength=0.1

Rscript scripts/subsample_loci.R $raw_bed $output_subsampled_bed $subsampling_strength
