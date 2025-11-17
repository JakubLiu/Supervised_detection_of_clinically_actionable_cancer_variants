#!/usr/bin/bash

random_seed="$1"
R1="$2"
R2="$3"
num_reads_to_keep="$4"
output_R1="$5"
output_R2="$6"

seqtk sample -s$random_seed $R1 $num_reads_to_keep | gzip -c > $output_R1
seqtk sample -s$random_seed $R2 $num_reads_to_keep | gzip -c > $output_R2