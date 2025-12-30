#!/usr/bin/python
 
import pandas as pd
import numpy as np
import sys
 
 
# function definitions _________________________________________________________________________________________________
def split_line(line, sample_names):
    n_sample_agnostic_rows = 5
    n_entries_per_sample = 4
    alt_chrom_pos_ref_cov = line[:n_sample_agnostic_rows]
    rest = line[n_sample_agnostic_rows:]
    n_bams = len(sample_names)
    splitted_rows = []
 
    start = 0
    end = n_entries_per_sample
    for i in range(0, n_bams):
        row = rest[start:end][:2]
        row = alt_chrom_pos_ref_cov + [sample_names[i]] + row
        splitted_rows.append(row)
        start = end
        end = end + n_entries_per_sample
 
    return splitted_rows
 
 
 
def translator(seq_in, ref, civic_alt):
    seq_out = ""
    n_alt = 0   # number of non-reference bases
    n_civic_alt = 0  # number of non-referecne bases that match the Civic alternative (only count SNVs)
    i = 0
 
    while i < len(seq_in):
        symbol = seq_in[i]
 
        if symbol in ['.', ',']:
            seq_out += ref
            i += 1
 
        elif symbol == '^':  # start of read (skip symbol + next char)
            i += 2
 
        elif symbol == '$':  # end of read
            i += 1
 
        elif symbol == '+':  # insertion
            # get the size of the insertion
            j = i + 1
            num_str = ""
            while j < len(seq_in) and seq_in[j].isdigit():
                num_str += seq_in[j]
                j += 1
            insertion_size = int(num_str)
 
            insertion_seq = seq_in[j:j+insertion_size]
            seq_out += insertion_seq
 
            i = j + insertion_size  # skip over inserted bases
            n_alt += insertion_size
 
        elif symbol == '-':  # deletion
            j = i + 1
            num_str = ""
            while j < len(seq_in) and seq_in[j].isdigit():
                num_str += seq_in[j]
                j += 1
            deletion_size = int(num_str)
 
            seq_out += 'DEL' * deletion_size
            i = j + deletion_size  # skip deleted bases
            n_alt += deletion_size
 
        elif symbol == '*':  # placeholder for deletion
            seq_out += 'DEL'
            i += 1
            n_alt += 1
 
        else:
            # actual nucleotide (A, T, G, C, a, t, g, c)
            if symbol.upper() in ['A', 'T', 'G', 'C', 'N']:
                seq_out += symbol.upper()
                n_alt += 1
                if symbol.upper() == civic_alt:
                    n_civic_alt += 1
                i += 1
 
            else:
                print(f"Unrecognized symbol at position {i}: {symbol}")
                i += 1  # to avoid infinite loop
 
    return (seq_out, n_alt, n_civic_alt)
 
# ____________________________________________________________________________________________________________________________________________
 
loci = sys.argv[1]
pileup_filepath = sys.argv[2]
output_filepath = sys.argv[3]
bamlist_filepath = sys.argv[4]
multi_vcf = sys.argv[5]
  
sample_names = pd.read_csv(multi_vcf, header=0, delimiter = ',').columns.tolist()


variants = pd.DataFrame(np.loadtxt(loci, delimiter = ',', skiprows = 0, dtype = str))
variants.columns = ['chrom', 'pos', 'ref', 'alt']
 
rows = []
with open(pileup_filepath, 'r') as file_in:
    for line in file_in:
        locus = line.rstrip("\n").split("\t")     # one line in the pileup file
        rows.append(locus)
 
 
pileup = pd.DataFrame(rows)
pileup = pileup.rename(columns = {0:'chrom', 1:'pos'})
merged = pd.merge(pileup, variants, on=["chrom", "pos"], how="inner")

 
cols = [
    'chrom', 'pos', 'ref', 'alt',
    2, 3, 4, 5, 6, 7, 8, 9,
    10, 11, 12, 13, 14, 15, 16, 17, 18
]
merged = merged[cols]


with open(output_filepath, 'w') as fout:
 
    # write the header line
    fout.write('sample,chrom,pos,ref,variant_alt,coverage,n_alt_reads,n_variant_supporting_alt_reads,translated_sequence' + '\n')

    for line in merged.to_numpy().tolist():
        sample_splitted_line = split_line(line, sample_names)
        for sample in sample_splitted_line:
            chromosome = str(sample[0])
            pos = int(sample[1])
            ref = str(sample[2])
            variant_alt = str(sample[3].split(' ')[0])
            sampleID = str(sample[5])
            raw_sequence = str(sample[-1])

            translated_sequence, n_alt, n_variant_alt = translator(raw_sequence,
                                                                 ref,
                                                                 variant_alt)
            
            coverage = len(translated_sequence)
            
            output_line = [sampleID, chromosome, pos, ref, variant_alt, coverage, n_alt, n_variant_alt, translated_sequence]
            output_line = ','.join([str(element) for element in output_line])
            fout.write(output_line + '\n')
