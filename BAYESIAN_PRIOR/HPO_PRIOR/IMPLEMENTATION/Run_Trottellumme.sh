#!/usr/bin/bash

tumor_bam='T1-DNA1-WES1.mutated.sorted.bam'
normal_validation_bam='N1-DNA1-WES1.bam'
negative_control_bamlist='/all_normals.txt'
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