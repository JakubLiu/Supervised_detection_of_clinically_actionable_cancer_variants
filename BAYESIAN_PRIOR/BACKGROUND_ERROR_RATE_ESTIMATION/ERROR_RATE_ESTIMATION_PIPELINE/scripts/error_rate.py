import numpy as np
import pandas as pd
import sys
from collections import Counter

processed_pileup = sys.argv[1]
output_file = sys.argv[2]

standard_bases = {"A", "T", "G", "C"}
epsilon = 1e-10  # add a pseudocunt for places with no coverage (avoid division by zero error)

with open(processed_pileup, 'r') as fin, open(output_file, 'w') as fout:

    fout.write('sampleID,chrom,start,stop,ref,coverage,n_alt_bases,read_pileup,norm_err_A,norm_err_T,norm_err_G,norm_err_C' + '\n')
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
        read_pileup = "".join(char for char in read_pileup if char in standard_bases)
        read_pileup_freq = Counter(read_pileup)  # returns a dict

        per_base_normalized_errors = {
            "A": 0,
            "T": 0,
            "G": 0,
            "C": 0
        }

        for base in standard_bases:
            normalized_read_count = read_pileup_freq[base]/(coverage+epsilon)  # calculate the number of reads supporting the given base, normalized by the coverage at that locus
            per_base_normalized_errors[base] = normalized_read_count

        per_base_normalized_errors[ref] = 0    # set the error rate for the reference base to 0

        output_line = [
            str(sampleID), str(chrom), str(start), str(stop), str(ref), str(coverage), str(n_alt_bases), str(read_pileup),
            str(per_base_normalized_errors['A']),
            str(per_base_normalized_errors['T']),
            str(per_base_normalized_errors['G']),
            str(per_base_normalized_errors['C'])
        ]

        output_line = ','.join(element for element in output_line)
        fout.write(output_line + '\n')

        

