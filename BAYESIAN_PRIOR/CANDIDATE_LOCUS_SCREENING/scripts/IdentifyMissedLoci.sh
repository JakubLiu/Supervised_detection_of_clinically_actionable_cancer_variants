#!/usr/bin/bash

# This script takes as input a gzvcf and a bed file containing the coordinates of a gene.
# It returns a list of loci corrdinates that are within that gene, but are not present in the VCF.

gene_bed="$1"   # cht   start   end of a gene
gzvcf="$2"      # the bgzipped vcf
output_loci_list="loci_in_${gene_bed}_not_called_in_${gzvcf}.txt"

tmp_bed_gene="${gene_bed}_per_base.bed"
tmp_bed_vcf="${gzvcf}.bed"

bedtools makewindows -b $gene_bed -w 1 > $tmp_bed_gene # list all loci in the bed

# covert the vcf into a bed
zgrep -v '^#' $gzvcf | awk '{OFS="\t"; print $1, $2-1, $2}' > $tmp_bed_vcf

# find the loci that have not been called in the vcf
bedtools intersect -v -a $tmp_bed_gene -b $tmp_bed_vcf > $output_loci_list

# clean up
rm $tmp_bed_gene
rm $tmp_bed_vcf


