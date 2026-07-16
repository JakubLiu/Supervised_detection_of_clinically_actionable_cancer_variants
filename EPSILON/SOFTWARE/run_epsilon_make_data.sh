#!/usr/bin/bash

# --------------------------------- alt specific mode -------------------------------------------
scripts/Epsilon_MakeData.sh \
    --bamlist negative_control_bamlist.txt \
    --loci_list Civic_actionable_SNV_loci.csv \
    --reference_genome hs37d5.fa \
    --alt_mode specific \
    --nranks 16 \
    --output_prefix NEGATIVE_CONTROL_DATA2