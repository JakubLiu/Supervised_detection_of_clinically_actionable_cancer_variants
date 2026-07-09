#!/usr/bin/bash

# alt_mode: 'generic' or 'specific'

./Epsilon_MakeData.sh \
    --bamlist minimal_bamlist.txt \
    --loci_list loci_list_specific.minimal.txt\
    --reference_genome hs37d5.fa \
    --alt_mode specific \
    --nranks 2 \
    --output_prefix sample_specific


