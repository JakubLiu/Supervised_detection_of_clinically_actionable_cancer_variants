#!/usr/bin/bash

./LRB_tumor_only.sh \
                --tumor_bam "T1-DNA1-WES1.mutated.sorted.bam" \
                --negative_control_bamlist "negative_control_cohort.txt" \
                --chromosome "7" \
                --start "55259515" \
                --stop "55259515" \
                --ref_allele "T" \
                --alt_allele "G" \
                --prior "0.005" \
                --posterior_cutoff "0.5" \
                --pseudocount "0.00001" \
                --output_call_file "output_call_file.txt" \
                --min_mapQ "30" \
                --min_baseQ "30" \
                --output_read_annotation_file "output_read_annotation_file.txt" \
                --reference_genome "hs37d5.fa" \
                --padding_upstream "10" \
                --padding_downstream "10" \
                --output_genomic_context_file "output_genomic_context_file.txt"

