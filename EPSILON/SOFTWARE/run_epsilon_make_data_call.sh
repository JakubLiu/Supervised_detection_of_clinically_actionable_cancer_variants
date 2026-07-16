#!/usr/bin/bash

# --------------------------------- alt specific mode -------------------------------------------
scripts/Epsilon_MakeData_call.sh \
    --bamlist tumor_bam.txt \
    --loci_list calling.bed.csv \
    --reference_genome hs37d5.fa \
    --alt_mode specific \
    --nranks 8 \
    --output_prefix tumor_data2