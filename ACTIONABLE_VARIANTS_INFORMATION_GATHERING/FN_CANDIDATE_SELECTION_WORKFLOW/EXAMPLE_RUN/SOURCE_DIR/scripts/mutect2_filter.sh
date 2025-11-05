#!/usr/bin/bash

# command line arguments passed from the level of snakemake
REF="$1"
INPUT_VCF="$2"
OUTPUT_VCF="$3"

gatk FilterMutectCalls -R $REF -V $INPUT_VCF -O $OUTPUT_VCF

