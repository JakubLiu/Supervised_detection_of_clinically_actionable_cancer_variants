#!/usr/bin/bash

# Rscript Trottellumme.R 'tumor.bam' 'bamlist.txt' '7' 100 100 'A' 'G' 0.01 0.5 0.00001 'output.txt'

tumor_bam='/data/cephfs-1/home/users/jali13_c/work/BAYESIAN_PRIOR/TUNE_PRIOR/STEP3_SIMULATE_SNV/TRAIN/MUTATED/inserted_SNV/bwa.BIH_165-T1-DNA1-WES1.mutated.sorted.bam'
normal_validation_bam='/data/cephfs-1/home/users/jali13_c/work/BAYESIAN_PRIOR/T-CELL-PROJECT-DATA/T_CELL-NEW/BIH_212-N1-DNA1-WES1/2021-02-16/out/bwa.BIH_212-N1-DNA1-WES1.bam'
negative_control_bamlist='/data/cephfs-1/home/users/jali13_c/work/BAYESIAN_PRIOR/TUNE_PRIOR/STEP3_SIMULATE_SNV/all_normals.txt'
chrom='7'
start=140453136
stop=140453136
ref='A'
alt='T'
prior=0.05
posterior_cutoff=0.5
pseudo=0.00001
out_tumor='output_tumor.txt'
out_validation='output_validation.txt'

Rscript Trottellumme.R $tumor_bam $negative_control_bamlist $chrom $start $stop $ref $alt $prior $posterior_cutoff $pseudo $out_tumor

Rscript Trottellumme.R $normal_validation_bam $negative_control_bamlist $chrom $start $stop $ref $alt $prior $posterior_cutoff $pseudo $out_validation