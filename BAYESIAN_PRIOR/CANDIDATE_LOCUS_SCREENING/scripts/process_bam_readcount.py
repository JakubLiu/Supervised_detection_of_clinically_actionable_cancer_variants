import numpy as np
import sys
 
 
input_file = sys.argv[1]
output_file = sys.argv[2]
 
with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
 
 
    fout.write('chrom,pos,ref,coverage,ref_count,total_alt_count,A_count,C_count,G_count,N_count' + '\n')
 
    for line in fin:
           
        line = line.split('\t')
 
        chrom = str(line[0])
        pos = str(line[1])
        ref = str(line[2])
        coverage = str(line[3])
        A_count = str(line[5].split(':')[1])
        C_count = str(line[6].split(':')[1])
        G_count = str(line[7].split(':')[1])
        T_count = str(line[8].split(':')[1])
        N_count = str(line[9].split(':')[1])
        read_count_dict = {'A' : A_count, 'C' : C_count, 'G' : G_count, 'T' : T_count, 'N' : N_count}
        ref_count = read_count_dict[ref]
        total_alt_count = int(coverage) - int(ref_count)
        total_alt_count = str(total_alt_count)
 
 
        output_line = [chrom, pos, ref, coverage, ref_count, total_alt_count, A_count, C_count, G_count, N_count]
        output_line = ','.join(output_line)
        fout.write(output_line + '\n')
