"""
Usage: python3 process_pileup.py raw_pileup output bamlist_file bamfile_extension
"""




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



def split_line2(line, sample_names):
    n_entries_per_sample = 3
    n_common_entries = 3

    chr_pos_ref = line[:n_common_entries]
    rest = line[n_common_entries:]

    n_bams = len(sample_names)
    splitted_rows = []

    start = 0
    end = n_entries_per_sample
    for i in range(0, n_bams):
        row = rest[start:end][:2]
        row = chr_pos_ref + [sample_names[i]] + row
        splitted_rows.append(row)
        start = end
        end = end + n_entries_per_sample
        
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

# ____________________________________________________________________________________________________________________________________________

pileup_filepath = sys.argv[1]
output_filepath = sys.argv[2]
bamlist_filepath = sys.argv[3]
bamfile_extension = sys.argv[4]


with open(bamlist_filepath) as f:  # extract the sample names based on the bamlist file
    sample_names = [
        line.strip().split('/')[-1].replace(bamfile_extension, '')
        for line in f
    ]


#pileup = pd.DataFrame(np.loadtxt(pileup_filepath, dtype = str))
pileup = pd.read_csv(pileup_filepath, dtype = str, sep = '\t', engine = 'python', on_bad_lines='skip')
fixed_cols = ['chrom', 'start', 'ref', 'coverage']  # sample-agnostic columns
pileup.columns = fixed_cols + list(range(len(fixed_cols),pileup.shape[1]))

with open(output_filepath, 'w') as fout:

    # write the header line
    fout.write('sample,chrom,start,stop,ref,coverage,n_alt_reads,translated_sequennce' + '\n')

    for line in pileup.to_numpy().tolist():  # iterate over the rows of the merged df
        sample_splitted_line = split_line2(line, sample_names)  # a list of lists, each inner list is the row for one sample
        
        # ['14', '105249145', 'A', '2', '.,', 'BIH_099-N1-DNA1-WES1.dedupa', 'FF', ']]']

        for sample in sample_splitted_line:  # iterate over each sample for the given actionable variant

            
            chromosome = sample[0]
            start = int(sample[1])
            stop = start
            ref = str(sample[2])
            sampleID = str(sample[3])
            coverage = int(sample[4])
            raw_sequence = str(sample[5])

            translated_sequence, n_alt = translator(raw_sequence,ref)
            
            output_line = [sampleID, chromosome, start, stop, ref, coverage, n_alt, translated_sequence]
            output_line = ','.join([str(element) for element in output_line])


            fout.write(output_line + '\n')