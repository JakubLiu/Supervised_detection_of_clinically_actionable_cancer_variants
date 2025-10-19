import numpy as np
import pandas as pd

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




def split_multibam_pileup(pileup_f):

    with open(pileup_f, 'r') as file:
    
        output = []
        for line in file:
            
            line = line.strip()
            lines_splitted = split_line(line)
            
            for sample in lines_splitted:
                output.append(list(sample))

    output = pd.DataFrame(output)
    output.columns = ['chromosome', 'position', 'ref', 'sample', 'coverage', 'reads']
    return output





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



def process_pileup(pileup_file):
    
    splitted_pileup_file = split_multibam_pileup(pileup_file)
    reads_bases = np.zeros(splitted_pileup_file.shape[0], dtype = object)
    n_alternative_bases = np.full(splitted_pileup_file.shape[0], -9, dtype = np.int64)

    i = 0
    for sequence in splitted_pileup_file['reads']:
        reference_base = splitted_pileup_file.iloc[i,2]
        reads_bases[i], n_alternative_bases[i] = translator(sequence, reference_base)
        i += 1

    splitted_pileup_file["bases_sequence"] = reads_bases
    splitted_pileup_file["n_alt_reads"] = n_alternative_bases
    return splitted_pileup_file



# usage ______________________________________________________________________________________________________________________________-

result = process_pileup(snakemake.input[0])
result.to_csv(snakemake.output[0], index = False)