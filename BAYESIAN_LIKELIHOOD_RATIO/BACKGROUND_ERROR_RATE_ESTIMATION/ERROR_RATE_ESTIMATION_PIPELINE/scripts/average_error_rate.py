import numpy as np
import sys


bamlist_file = sys.argv[1]
error_rate_file = sys.argv[2]
output_error_rate_file = sys.argv[3]


bamlist_file = np.loadtxt(bamlist_file, dtype = str)
num_samples_per_locus = bamlist_file.shape[0]  # the average error rate across these many samples will be calculated




with open(error_rate_file, 'r') as input_file, open(output_error_rate_file, 'w') as output_file:

    output_file.write('chrom,start,stop,ref,coverage,mean_errA,mean_errT,mean_errG,mean_errC,mean_errN' + '\n')
    # [chrom, start, stop, ref, mean_errA, mean_errT, mean_errG, mean_errC, mean_errN]

    sample_counter = sum_errA = sum_errT = sum_errG = sum_errC = sum_errN = 0 # set all the counters to 0
    next(input_file)

    for line in input_file:  # loop over the input file line by line

        line = line.rstrip("\n").split(",")

        # rea the necessary columns
        chrom = str(line[1])
        start = int(line[2])
        stop = int(line[3])
        ref = str(line[4])
        coverage = int(line[5])
        errA = float(line[8])
        errT = float(line[9])
        errG = float(line[10])
        errC = float(line[11])
        errN = float(line[12])

        # update the sums of the error rates
        sum_errA = sum_errA + errA
        sum_errT = sum_errT + errT
        sum_errG = sum_errG + errG
        sum_errC = sum_errC + errC
        sum_errN = sum_errN + errN

        sample_counter = sample_counter + 1

        if sample_counter == num_samples_per_locus:     # if we reach the end of a block then calcaulte the means and write a line to the output file

            mean_errA = sum_errA/num_samples_per_locus
            mean_errT = sum_errT/num_samples_per_locus
            mean_errG = sum_errG/num_samples_per_locus
            mean_errC = sum_errC/num_samples_per_locus
            mean_errN = sum_errN/num_samples_per_locus

            outline = [chrom, start, stop, ref, coverage, mean_errA, mean_errT, mean_errG, mean_errC, mean_errN]
            output_file.write(','.join(str(x) for x in outline) + '\n')

            # reset the counters
            sum_errA = sum_errT = sum_errG = sum_errC = sum_errN = 0
            sample_counter = 0
        
        

            
            

        


