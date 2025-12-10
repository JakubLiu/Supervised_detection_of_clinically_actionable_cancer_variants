from numba import njit
from Bio import SeqIO
import numpy as np
import sys
from mpi4py import MPI
import os

@njit
def GC_content(sequence):  # the function accepts a sequence of bytes

    # convert the character representation of the base to its ASCII code
    A_ascii = ord('A')
    T_ascii = ord('T')
    G_ascii = ord('G')
    C_ascii = ord('C')

    A_count = T_count = G_count = C_count = 0
    
    for base in sequence:
        if base == A_ascii:
            A_count = A_count + 1
        elif base == T_ascii:
            T_count = T_count + 1
        elif base == G_ascii:
            G_count = G_count + 1
        elif base == C_ascii:
            C_count = C_count + 1

    total_count = A_count + T_count + G_count + C_count

    if total_count == 0.0:
        return 0.0
    
    GC_count = G_count + C_count
    GC_content_ = GC_count / total_count

    return GC_content_



def Calc_Shannon_entropy(ref, chrom, windowsize, output_file):

    common_world = MPI.COMM_WORLD
    rank = common_world.Get_rank()
    size = common_world.Get_size()


    # let the root rank (0) load the chromosome fasta
    if rank == 0:
        sequence = None
        
        for rec in SeqIO.parse(ref, "fasta"):

            if rec.id == str(chrom):
                sequence = str(rec.seq).upper().replace('N', '')  # remove N bases
    else:
        sequence = None

    # send the chromosome fasta sequence to all other CPU's
    sequence = common_world.bcast(sequence, root = 0)
    num_total_windows = len(sequence) - int(windowsize) + 1 # how many sliding windows we will have in total (assuming a stride of 1)
    num_windows_per_CPU = num_total_windows // size  # how many windows will each CPU need to handle
    start_window = rank * num_windows_per_CPU  # the 1st window for the given CPU

    # the last window for the given CPU (this depends on the rank of the CPU)
    if rank < (size-1):
        end_window = (rank+1) * num_windows_per_CPU
    else:
        end_window = num_total_windows

    num_windows_per_CPU = end_window - start_window
    report_interval = max(1, num_windows_per_CPU // 1000) 


    tmp = f'{output_file}.rank{rank}.tmp'
    with open(tmp, 'w') as fout_tmp:
        for idx, i in enumerate(range(start_window, end_window)):
            window = sequence[i:i+int(windowsize)] # actuall extract the window from the sequence
            window_bytes = np.frombuffer(window.encode('ascii'), dtype=np.uint8)
            GC = GC_content(window_bytes)
            fout_tmp.write(f'{chrom},{i+1},{i+int(windowsize)},{GC:.6f}\n')

            if (idx + 1) % report_interval == 0 or (idx + 1) == num_windows_per_CPU:
                percent_done = int((idx + 1) / num_windows_per_CPU * 100)
                print(f"Rank {rank}: {percent_done}% done\n", end='\r', flush = True)

    common_world.Barrier() # wait for all CPUs _______________________________________________________________________________________________________________

    # let the root rank (0) write the header and concatename all the temporary files to the final output file
    if rank == 0:
        with open(output_file, 'w') as fout:
            fout.write('CHROMOSOME,START,STOP,GC_CONTENT\n')

            for current_rank in range(0, size):

                # open the previously created tempoprary files (one file = one CPU) and write their contents to the final output file
                tmp = f'{output_file}.rank{current_rank}.tmp'
                with open(tmp, 'r') as tmp_file:
                    
                    for line in tmp_file:
                        fout.write(line)
                
                os.remove(tmp)
