#!/usr/bin/bash

callable_loci_bed="callable.bed"
input_bed="subsampled.mpileup.txt"
output_bed="subsampled.mpileup.final.txt"

scripts/remove_callable_loci.sh $callable_loci_bed $input_bed $output_bed
