import sys
import numpy as np

csv_in = sys.argv[1]
vcf_out = sys.argv[2]


csv = np.loadtxt(csv_in, delimiter = ',', dtype = str)

with open(vcf_out, 'w') as vcf:

    vcf.write("##fileformat=VCFv4.2\n")
    vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    for i in range(0, csv.shape[0]):
        vcf.write(f'{csv[i,0]}\t{csv[i,1]}\t.\t{csv[i,3]}\t{csv[i,4]}\t.\t.\t.\n')


