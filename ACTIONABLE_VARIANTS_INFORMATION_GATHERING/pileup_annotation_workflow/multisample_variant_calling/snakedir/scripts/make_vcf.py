import sys
import numpy as np

csv_in = sys.argv[1]
vcf_out = sys.argv[2]


csv = np.loadtxt(csv_in, delimiter = ',', dtype = str, usecols = (0), skiprows = 1)

with open(vcf_out, 'w') as vcf:

    vcf.write("##fileformat=VCFv4.2\n")
    vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    for i in range(0, csv.shape[0]):

        line = str(csv[i]).split(';')
        vcf.write(f'{line[0]}\t{line[1]}\t.\t{line[2]}\t{line[3]}\t.\t.\t.\n')