import numpy as np
import pandas as pd
import gc
import gzip

# function definitions _________________________________________________________________________________________________
def split_line(line):
    line = line.split("\t")
    chrom_pos_ref = line[:3]
    rest = line[3:]
    n_bams = int(len(rest)/4)
    splitted_rows = []

    start = 0
    end = 4
    for i in range(0, n_bams):
        row = rest[start:end][:2]
        row = chrom_pos_ref + ['sample' + str(i+1)] + row
        splitted_rows.append(row)
        start = end
        end = end + 4
        
    return splitted_rows



def translator(seq_in, ref):
    seq_out = ""
    n_alt = 0   # number of non-reference bases
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
                i += 1
                
            else:
                print(f"Unrecognized symbol at position {i}: {symbol}")
                i += 1  # to avoid infinite loop

    return (seq_out, n_alt)



def process_pileup(raw_pileup_file_path, output_file_path):
    
    with open(raw_pileup_file_path, "r", encoding='utf-8') as infile, \
        open(output_file_path, "w", encoding="utf-8") as outfile:

        outfile.write('CHROM,POS,REF,SAMPLE,COVERAGE,READ_BASES,N_ALT_BASES' + '\n') # write the header

        # iterate over the lines of the input file, do some stuff and write that processed line to disc
        for line in infile:  # here one line has multiple samples
            
            line = line.strip()
            lines_splitted = split_line(line)

            # iterate over the samples in each line of the input file
            for sample in lines_splitted:    # here one line is one sample
                raw_sequence = sample[-1]
                translated_sequence, n_alternative_bases = translator(raw_sequence, sample[2])
                output_line = sample[:-1] + [translated_sequence, n_alternative_bases]
                output_line = ','.join([str(element) for element in output_line])  # change every element of the list to a string and then change the list to a string (sep=',')
                outfile.write(output_line + '\n') # write a single line/sample to disc

                # garbage collection just to make sure
                del raw_sequence
                del translated_sequence
                del n_alternative_bases
                del output_line
            
            # garbage collection just to make sure
            del line
            del lines_splitted
            gc.collect()



# usage ______________________________________________________________________________________________________________________________-

process_pileup(raw_pileup_file_path=snakemake.input[0],
               output_file_path=snakemake.output[0])

