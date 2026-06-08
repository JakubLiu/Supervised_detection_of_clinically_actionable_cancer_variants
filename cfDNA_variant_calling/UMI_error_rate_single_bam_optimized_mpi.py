#!/usr/bin/env python

import pysam
import pandas as pd
import numpy as np
from mpi4py import MPI
from collections import defaultdict, Counter
import argparse

CHUNK_SIZE = 10_000

# function to process each chunk of the bam file
# returns a list of loci and their corresponding error rates
def process_region(bam, chrom, start, end, min_reads_per_UMI, min_coverage):
    bamfile = pysam.AlignmentFile(bam, "rb")

    loci = []
    error_rates = []

    # here iterate over the loci
    for pileupcolumn in bamfile.pileup(
        chrom,
        start,
        end,
        truncate=True
    ):

        families = defaultdict(list)

        locus = pileupcolumn.reference_pos
        coverage = 0

        # here iterate over the the reads that align to a give locus
        for pileupread in pileupcolumn.pileups:

            read = pileupread.alignment

            if read.has_tag("MI") and not pileupread.is_del:

                mi = read.get_tag("MI")

                if pileupread.query_position is not None:
                    base = read.query_sequence[
                        pileupread.query_position
                    ]

                    families[mi].append(base)

            coverage += 1

        if coverage < min_coverage:
            continue

        local_error_rates = []
        weights = []

        # iterate over the reads that share the same UMI (and have mapped to the same locus)
        for mi, bases in families.items():

            if len(bases) < min_reads_per_UMI:
                continue

            majority_allele_count = (Counter(bases).most_common(1)[0][1])  # most common (major) allele (I dont care that is the reference allele here)
            non_major_alleles = (len(bases) - majority_allele_count)  # all other alleles (I dont care what is the reference allele here)
            umi_error_rate = (non_major_alleles / len(bases)) # just the number of non-major alleles divided by the coverage
            local_error_rates.append(umi_error_rate)
            weights.append(len(bases))   # if UMIx is represented by 30 bases and UMIy is represented by 10 bases, then the error rate estimation based on UMIx should have more weight

        if len(weights) > 0:
            weighted_error_rate = np.average(local_error_rates,weights=weights)  # this is the weighted mean
            loci.append(int(locus))
            error_rates.append(weighted_error_rate)

    bamfile.close()

    return loci, error_rates


# this function builds the regions for MPI
def build_regions(bam, chunk_size):

    bamfile = pysam.AlignmentFile(bam, "rb")

    regions = []

    for chrom, length in zip(bamfile.references, bamfile.lengths):

        for start in range(0, length, chunk_size):

            end = min(start + chunk_size, length)

            regions.append((chrom, start, end))

    bamfile.close()

    return regions


# this is a wrapper, it runs the process_regions() function for each region that is returned by the build_regions() function
def UMI_error_single_bam(bam, samplename, min_reads_per_UMI, min_coverage, chunk_size=CHUNK_SIZE):

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # build the regions on rank 0 and bcast them to all other ranks
    if rank == 0:
        regions = build_regions(bam,chunk_size)
    else:
        regions = None

    regions = comm.bcast(regions,root=0)

    # point to the region chunk that should be processed by the current rank
    regions_for_this_rank = regions[rank::size]

    local_loci = []
    local_error_rates = []

    # process the region chunk
    for chrom, start, end in regions_for_this_rank:
        loci, errs = process_region(bam=bam, chrom=chrom, start=start, end=end, min_reads_per_UMI=min_reads_per_UMI, min_coverage=min_coverage)
        local_loci.extend(loci)
        local_error_rates.extend(errs)

    # gather the results from all workers/ranks on rank 0
    gathered_loci = comm.gather(local_loci,root=0)
    gathered_errors = comm.gather(local_error_rates,root=0)

    if rank == 0:

        loci = [
            x
                for sublist in gathered_loci
                    for x in sublist
        ]

        errors = [
            x
                for sublist in gathered_errors
                    for x in sublist
        ]

        sample = (
            samplename
            if samplename is not None
            else bam
        )

        results = pd.DataFrame({
            "locus": loci,
            "error_rate": errors,
            "samplename": sample
        })

        return results

    return None


# get the coomand line arguments and run that
def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--bam",
        required=True,
        type=str
    )

    parser.add_argument(
        "--samplename",
        default=None,
        type=str
    )

    parser.add_argument(
        "--min_reads_per_UMI",
        default=5,
        type=int
    )

    parser.add_argument(
        "--min_coverage",
        default=20,
        type=int
    )

    parser.add_argument(
        "--chunk_size",
        default=10000,
        type=int
    )

    parser.add_argument(
        "--output",
        default=None,
        type=str
    )

    args = parser.parse_args()

    results = UMI_error_single_bam(
        bam=args.bam,
        samplename=args.samplename,
        min_reads_per_UMI=args.min_reads_per_UMI,
        min_coverage=args.min_coverage,
        chunk_size=args.chunk_size
    )

    rank = MPI.COMM_WORLD.Get_rank()

    if rank == 0:

        output = (
            args.output
            if args.output is not None
            else (
                f"{args.samplename}.UMI_error_rates.tsv"
                if args.samplename
                else "UMI_error_rates.tsv"
            )
        )

        results.to_csv(
            output,
            sep="\t",
            index=False
        )

        print(
            f"Wrote {len(results)} loci to {output}"
        )


if __name__ == "__main__":
    main()