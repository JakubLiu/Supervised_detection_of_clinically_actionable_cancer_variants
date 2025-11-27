#!/usr/bin/python3
import re
import os
import pysam
import argparse
import csv


def extract_basename(path):
    filename = os.path.basename(path)
    name = re.sub(r'\..*\.bam$', '', filename)
    return name

def MultiVcfVariantCompilation(vcf_list_path, output_path):
    # Read all VCF file paths
    with open(vcf_list_path) as f:
        vcf_paths = [line.strip() for line in f if line.strip()]

    n_files = len(vcf_paths)
    all_variants = {}
    
    # Process each VCF
    filenames = []
    for file_idx, path in enumerate(vcf_paths, start=1):
        filenames.append(extract_basename(path))
        vcf = pysam.VariantFile(path)
        seen = set()
        for variant in vcf.fetch():
            variant_ID = f"{variant.chrom};{variant.pos};{variant.ref};{variant.alts[0]}"
            if variant_ID not in seen:
                seen.add(variant_ID)
                if variant_ID not in all_variants:
                    all_variants[variant_ID] = ['0'] * n_files
                all_variants[variant_ID][file_idx - 1] = '1'

    # Write to CSV
    header = ['variant'] + filenames
    with open(output_path, 'w', newline='') as out_f:
        writer = csv.writer(out_f)
        writer.writerow(header)
        for variant, presence in all_variants.items():
            writer.writerow([variant] + presence)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compile variants across multiple VCF files.")
    parser.add_argument("--input", "-i", required=True, help="Path to text file containing VCF paths")
    parser.add_argument("--output", "-o", required=True, help="Path to output CSV")
    args = parser.parse_args()
    
    MultiVcfVariantCompilation(args.input, args.output)