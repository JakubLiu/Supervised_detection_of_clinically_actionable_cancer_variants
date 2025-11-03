#!/usr/bin/bash

# command line arguments passed from the level of snakemake
BAM="$1"
REF="$2"
VCF="$3"
THREADS="${4:-8}"  # default to 8 threads if not provided


echo "Running Mutect2 on $BAM..."
gatk Mutect2 \
  -R "$REF" \
  -I "$BAM" \
  -O "$VCF" \
  --native-pair-hmm-threads "$THREADS"
