#!/usr/bin/bash


python3 UMI_error_single_bam_new.py \
                --bam "mapped.grouped.sorted.bam" \
                --bed "example.bed" \
                --samplename "example_sample" \
                --min_reads_per_UMI 10 \
                --min_depth 20 \
                --max_PCR_error_rate 0.01 \
                --output "test_output.csv"
