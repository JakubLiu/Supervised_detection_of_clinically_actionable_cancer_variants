#!/usr/bin/bash

nranks="8"


mpirun -np $nranks python3 scripts/extract_features.py \
    --bed loci.bed \
    --reference_genome reference.fa \
    --bam reads.bam \
    --output output_dir \
    --genomic_window 100
