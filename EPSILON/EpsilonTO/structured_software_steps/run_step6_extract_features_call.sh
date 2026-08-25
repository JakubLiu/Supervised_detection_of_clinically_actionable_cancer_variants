#!/usr/bin/bash

nranks="8"

mpirun -np $nranks python3 scripts/extract_features_call.py \
    --bed calling_region.bed \
    --reference-genome reference.fa \
    --bam tumor.bam \
    --output call_data_output_dir \
    --genomic-window 100
