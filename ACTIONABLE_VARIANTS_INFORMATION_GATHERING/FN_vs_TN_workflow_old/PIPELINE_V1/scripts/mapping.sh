#!/usr/bin/bash

REF="$1"
R1="$2"
R2="$3"
OUT="$4"
THREADS="$5"

bwa mem -t "$THREADS" "$REF" "$R1" "$R2"  | samtools view -b -o "$OUT"
