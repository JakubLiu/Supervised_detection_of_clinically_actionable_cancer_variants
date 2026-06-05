from collections import defaultdict, Counter
import pysam
import numpy as np
import pandas as pd


def UMI_error_single_bam(bam,samplename = None, min_reads_per_UMI = 10, min_coverage = 20):

    grouped_mapped_bam = pysam.AlignmentFile(bam)
    loci = []
    erorr_rates = []



    for pileupcolumn in grouped_mapped_bam.pileup():

        families = defaultdict(list)
        locus = pileupcolumn.reference_pos

        coverage = 0
        for pileupread in pileupcolumn.pileups:

            read = pileupread.alignment

            if read.has_tag("MI"):
                mi = read.get_tag("MI")

                if not pileupread.is_del:
                    base = read.query_sequence[pileupread.query_position]
                    families[mi].append(base)

            coverage += 1

        # if the min coverage is satisfied then continue with the UMI-based error rate estimation
        local_error_rates = []
        weights = []   # weights are the numer of reads per UMI tag
        if coverage >= min_coverage:
            for mi, bases in families.items():

                # if the minimum number of reads per UMI is satisfied then continue with the UMI-based error rate estimation
                if len(bases) >= min_reads_per_UMI:

                    majority_allele, majority_allele_count = Counter(bases).most_common(1)[0] # get the most common allele for the locus x UMI combination
                    non_major_alleles = len(bases) - majority_allele_count
                    per_locus_per_UMI_error_rate = non_major_alleles/len(bases) # the local error rate for the given UMI
                    local_error_rates.append(per_locus_per_UMI_error_rate)
                    weights.append(len(bases))

        if len(weights) > 0:
            erorr_rates.append(np.average(local_error_rates, weights = weights))
            loci.append(int(locus))

    if samplename == None:
        sample = [bam] * len(loci)
    else:
        sample = [samplename] * len(loci)

    results = pd.DataFrame([loci, erorr_rates, sample]).T
    results.columns = ["locus", "error_rate", "samplename"]
    
    return results

'''
==================================== usage example =========================================

x = UMI_error_single_bam(
    bam = 'minimal.grouped.sorted.bam',
    samplename='sample',
    min_reads_per_UMI=1,
    min_coverage=1
)

print(x)

  locus error_rate samplename
0    140753409        0.0     sample
1    140753410        0.0     sample
2    140753411        0.0     sample
3    140753412        0.0     sample
4    140753413        0.0     sample
..         ...        ...        ...
379  140753930        0.0     sample
380  140753931        0.0     sample
381  140753932        0.0     sample
382  140753933        0.0     sample
383  140753934        0.0     sample

'''
