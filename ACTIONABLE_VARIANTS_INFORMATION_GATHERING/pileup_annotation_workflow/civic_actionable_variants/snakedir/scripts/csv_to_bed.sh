#!/usr/bin/bash

csv_in="$1"
bed_out="$2"

awk -F',' '
BEGIN { OFS="\t" }
NR > 1 && $1 != "" {
    chrom = $1
    start = $2 - 1
    end   = $2
    print chrom, start, end
}' "$csv_in" > "$bed_out"
