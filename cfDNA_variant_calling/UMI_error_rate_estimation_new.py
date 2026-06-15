from collections import defaultdict, Counter
import pysam
import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser()


parser.add_argument("--bam", required=True, type=str)
parser.add_argument("--bed", required=True, type = str)
parser.add_argument("--samplename", required=False, type = str, default = None)
parser.add_argument("--min_reads_per_UMI", required=False, type = int, default=10)
parser.add_argument("--min_depth", required=False, type = int, default=20)
parser.add_argument("--max_PCR_error_rate", required=False, type = float, default=0.01)
parser.add_argument("--output", required=False, type = str, default='output.csv')
args = parser.parse_args()


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
    min_reads_per_UMI,
    min_coverage,
    max_PCR_error_rate,
    samplename=None
):

    bamfile = pysam.AlignmentFile(bam, "rb")
    regions = load_bed(bed)

    loci = []
    error_rates = []
    PCR_errors = []

    # iterate over the regions in the bed file
    for chrom, start, end in regions:

        for pileupcolumn in bamfile.pileup(
            chrom,
            start,
            end,
            truncate=True
        ):

            locus = str(chrom) + ';' + str(pileupcolumn.reference_pos)

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
                er_PCR = np.nan

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
        "seq_error_rate": error_rates,
        "samplename": sample
    })

    return results




if __name__ == "__main__":

    x = UMI_error_single_bam(
        bam=args.bam,
        bed=args.bed,
        samplename=args.samplename,
        min_reads_per_UMI=args.min_reads_per_UMI,
        min_coverage=args.min_depth,
        max_PCR_error_rate=args.max_PCR_error_rate
    )


    x.to_csv(args.output, index=False)


"""
============= structure of the bed file: ===========================
6	29910192	29911329
           ...

"""
