#!/usr/bin/bash

dirname="$1"
sample_name="$2"
output="$3"


grep "callable" "$dirname/"*"$sample_name"*.stats | awk '{sum += $NF} END {printf "statistic\tvalue\ncallable\t%.1f\n", sum}' > "$output"