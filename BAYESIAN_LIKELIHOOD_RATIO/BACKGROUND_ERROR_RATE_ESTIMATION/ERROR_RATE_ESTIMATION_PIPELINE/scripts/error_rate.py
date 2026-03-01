import numpy as np
import pandas as pd
import sys
from collections import Counter

processed_pileup = sys.argv[1]
output_file = sys.argv[2]

bases = {"A", "T", "G", "C", "N"}
epsilon = 1e-10  # add a pseudocunt for places with no coverage (avoid division by zero error)

with open(processed_pileup, 'r') as fin, open(output_file, 'w') as fout:

    fout.write('sampleID,chrom,start,stop,ref,coverage,n_alt_bases,read_pileup,norm_err_A,norm_err_T,norm_err_G,norm_err_C,norm_err_N' + '\n')
    next(fin)  # skip the header of the input file

    for input_line in fin:

        input_line = input_line.split(',')
        
        sampleID = str(input_line[0])
        chrom = int(input_line[1])
        start = int(input_line[2])
        stop = start
        ref = str(input_line[4])
        coverage = int(input_line[5])
        n_alt_bases = int(input_line[6])
        read_pileup = str(input_line[7])

        # remove all non-standars bases (because of this the error rate calculation will only work for SNVs (not indels))
        read_pileup = "".join(char.upper() for char in read_pileup if char.upper() in bases)
        read_pileup_freq = Counter(read_pileup)  # returns a dict

        per_base_normalized_errors = {
            "A": 0,
            "T": 0,
            "G": 0,
            "C": 0,
            "N": 0
        }

        for base in bases:
            normalized_read_count = read_pileup_freq[base]/(coverage+epsilon)  # calculate the number of reads supporting the given base, normalized by the coverage at that locus
            per_base_normalized_errors[base] = normalized_read_count

        per_base_normalized_errors[ref] = 0    # set the error rate for the reference base to 0

        output_line = [
            sampleID,chrom,start,stop,ref,coverage,n_alt_bases,read_pileup,
            per_base_normalized_errors['A'],
            per_base_normalized_errors['T'],
            per_base_normalized_errors['G'],
            per_base_normalized_errors['C'],
            per_base_normalized_errors['N']
        ]

        output_line = ','.join(str(element) for element in output_line)
        fout.write(output_line + '\n')

        

