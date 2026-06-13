from collections import defaultdict, Counter
import pysam
import numpy as np
import pandas as pd


def load_bed(bed_file):
    regions = []
    with open(bed_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            chrom, start, end = line.strip().split("\t")[:3]
            regions.append((chrom, int(start), int(end)))
    return regions


def UMI_error_single_bam(
    bam,
    bed,
    samplename=None,
    min_reads_per_UMI=10,
    min_coverage=20,
    max_PCR_error_rate=0.1
):

    bamfile = pysam.AlignmentFile(bam, "rb")
    regions = load_bed(bed)

    loci = []
    error_rates = []

    # iterate over the regions in the bed file
    for chrom, start, end in regions:

        for pileupcolumn in bamfile.pileup(
            chrom,
            start,
            end,
            truncate=True
        ):

            locus = pileupcolumn.reference_pos

            families = defaultdict(list)
            coverage = 0

            # for the given locus gather information about the UMI families......................................
            for pileupread in pileupcolumn.pileups:

                read = pileupread.alignment

                if read.has_tag("MI") and not pileupread.is_del:
                    mi = read.get_tag("MI")

                    base = read.query_sequence[pileupread.query_position]
                    families[mi].append(base)

                coverage += 1

            local_error_rates = []
            weights = []

            # the minimum coverage filter...................................................
            if coverage >= min_coverage:

                # iterate over reads inside the UMI family
                for mi, bases in families.items():

                    # the minimum reads per UMI filter.......................................
                    if len(bases) >= min_reads_per_UMI:

                        majority_count = Counter(bases).most_common(1)[0][1]
                        err = (len(bases) - majority_count) / len(bases)

                        # the minimum error filter .........................................
                        """
                        Basically if the UMI family has a strongly dominant major allele (i.e. err is small),
                        then the PCR error is low. Then the remaining error is saved and treated as the
                        sequencing error.
                        """
                        if err <= max_PCR_error_rate:
                            local_error_rates.append(err)
                            weights.append(len(bases))


            if len(local_error_rates) > 0:
                er = np.average(local_error_rates, weights=weights)
            else:
                er = np.nan

            loci.append(locus)
            error_rates.append(er)

    bamfile.close()

    # sample column must match loci length
    if samplename is None:
        sample = [bam] * len(loci)
    else:
        sample = [samplename] * len(loci)

    results = pd.DataFrame({
        "locus": loci,
        "error_rate": error_rates,
        "samplename": sample
    })

    return results





if __name__ == "__main__":

    x = UMI_error_single_bam(
        bam="mapped.grouped.sorted.bam",
        bed="example.bed",
        samplename="sample",
        min_reads_per_UMI=10,
        min_coverage=3,
        max_PCR_error_rate=0.00239
    )

    print(x)

    x.to_csv("umi_error_results.csv", index=False)


"""
============= structure of the bed file: ===========================
6	29910192	29911329
           ...

"""
