#!/usr/bin/python

import numpy as np
import pysam
import argparse

def MultiVcfVariantCompilation(vcf_list_path):
    with open(vcf_list_path, 'r') as f:
        count = 1
        tmp = np.zeros((1,2), dtype='U100')
        n_files = 0

        # loop over the VCF files
        for path in f:
            n_files += 1
            vcf = pysam.VariantFile(path.strip())
            variant_IDs = []
            vcf_ID = []

            # loop over the variants within a given VCF file
            for variant in vcf.fetch():
                variant_ID = f"{variant.chrom};{variant.pos};{variant.ref};{variant.alts[0]}"
                variant_IDs.append(variant_ID)
                vcf_ID.append(str(count))

            variant_IDs = list(set(variant_IDs))  # remove duplicates in VCF if any
            count += 1

            variant_IDs_updated = np.column_stack((variant_IDs, vcf_ID))
            tmp = np.concatenate((tmp, variant_IDs_updated), axis=0)

        all_files_variant_IDs = tmp[1:,:]
        unique_variants_across_files = np.unique(all_files_variant_IDs[:,0])

        # create the output
        output = np.full((len(unique_variants_across_files), n_files+1), '-9', dtype='U100')

        for k, query in enumerate(unique_variants_across_files):
            result = all_files_variant_IDs[all_files_variant_IDs[:,0] == query]
            output[k,0] = query
            output[k,1:] = '0'  # initialize presence/absence with 0
            for idx in result[:,1]:
                output[k,int(idx)] = '1'

    output = output[~np.any(output == '-9', axis=1)]
    header = ['variant'] + ['file_'+str(i) for i in range(1, n_files + 1)]
    output = np.vstack((np.array(header), output))
#    output = np.vstack((np.array(['variant', list(range(1,(n_files+1)))]), output))

    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compile variants across multiple VCF files.")
    parser.add_argument("--input", '-i',  help="Path to a text file containing VCF file paths (one per line)", required = True)
    parser.add_argument("--output", '-o', help="Name of the output file", required = True)
    args = parser.parse_args()

    output = MultiVcfVariantCompilation(args.input)
    print(output)
    print(f'saving output to: {args.output}')
    np.savetxt(args.output, output, delimiter=",", fmt="%s")
    print('done.')

