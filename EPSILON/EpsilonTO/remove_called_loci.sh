#!/usr/bin/bash

blacklist="$1"
feature_extraction_bed="$2"
output="$3"

awk '{print $2}' $blacklist > pos.txt

grep -vf pos.txt $feature_extraction_bed > $output

rm pos.txt

echo "done."
