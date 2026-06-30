#!/usr/bin/bash


INPUT=$1
OUTPUT=$2

bcftools view -f PASS -Oz -o "$OUTPUT" "$INPUT"

bcftools index -t "$OUTPUT"