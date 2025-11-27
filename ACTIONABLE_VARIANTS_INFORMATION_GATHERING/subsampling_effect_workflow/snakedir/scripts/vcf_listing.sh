#!/usr/bib/bash

dirname="$1"
samplename="$2"
output="$3"

#ls "$dirname" | xargs realpath | grep "\.vcf\.gz$" | grep "$samplename" > "$output"
find "$dirname" -maxdepth 1 -type f -name "*${samplename}*.vcf.gz" -printf "%p\n" > "$output"