#!/usr/bin/bash

# proportion is the number of files in which the variant was called
# so it was missed in total-proportion

import argparse

def VariantExtractor(file_in, file_out, proportion):
    
    proportion = int(proportion)

    with open(file_in, 'r') as f_in, open(file_out, 'w') as f_out:
        i = 0
        for line_in in f_in:
            line_in = line_in.strip()

            if i == 0:
                i = i + 1
                continue


            occurrences = line_in.split(',')[1:]
            occurrences = [int(x) for x in occurrences]

            if sum(occurrences) == int(proportion):
                chrom = line_in.split(',')[0].split(';')[0]
                start = str(int(line_in.split(',')[0].split(';')[1]) - 1)  # ensure 0-based positions
                end = str(int(start) + 1)
                line_out = '\t'.join([chrom, start, end])
                f_out.write(line_out + '\n')

            i = i + 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a .bed file out of the variant that have been called in m/n samples.")
    parser.add_argument("--input", "-i", required=True, help="Path to the multi vcf compilation .csv file")
    parser.add_argument("--proportion", "-p", required=True, help="In how many samples the variant should have been called")
    parser.add_argument("--output", "-o", required=True, help="Path to output .bed")
    args = parser.parse_args()

    VariantExtractor(args.input, args.output, args.proportion)
