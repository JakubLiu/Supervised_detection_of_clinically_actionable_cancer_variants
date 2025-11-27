#!/usr/bin/bash

input_bam="$1"
tmp_bam="$2"
normal_sample_name="$3"  # this normal sample name will be also used for the tumor bams
threads="$4"
output_bam="$5"


gatk AddOrReplaceReadGroups \
    -I "$input_bam" \
    -O "$tmp_bam" \
    --RGID "id1" \
    --RGLB "lib1" \
    --RGPL "ILLUMINA" \
    --RGPU "unit1" \
    --RGSM "$normal_sample_name"


samtools sort -@ "$threads" -o "$output_bam" "$tmp_bam"

samtools index "$output_bam"

# usage
# ./process_bams.sh input_bam tmp_bam normal_sample_name threads output_bam