#!/usr/bin/bash

# =============================== how to run ====================================
# ./LoFreqSomaric.sh $T $N $ref $bed output_ $threads
# then look for output_somatic_final.snvs.vcf.gz
# ===============================================================================


tumor="$1"
normal="$2"
ref="$3"
bed="$4"
out="$5"
threads="$6"



lofreq somatic -n $normal \
                -t $tumor \
                -f $ref \
                -l $bed \
                --threads $threads \
                -o $out