#!/usr/bin/python

from mpi4py import MPI
import pysam
import numpy as np
import os
import argparse
from itertools import product

from extractor_functions import (
    read_bed,
    extract_read_tensor,
    extract_reference_context,
    calculate_error_rate,
    parse_args
)




def main():

    args = parse_args()

    comm = MPI.COMM_WORLD

    rank = comm.Get_rank()
    size = comm.Get_size()


    loci = read_bed(args.bed)

    # run the loci in parallel
    my_loci = loci[rank::size]


    fasta = pysam.FastaFile(
        args.reference_genome
    )

    bam = pysam.AlignmentFile(
        args.bam,
        "rb"
    )


    rank_dir = os.path.join(
        args.output,
        f"rank{rank:03d}"
    )

    os.makedirs(
        rank_dir,
        exist_ok=True
    )

    counter = 0
    # iterate over the loci (in the given rank)
    for chrom, start, stop, ref in my_loci:

        progress = counter/len(my_loci)*100
        print(f'rank {rank}, progress: {progress}%', flush = True)
        counter = counter + 1

        context = extract_reference_context(
            fasta,
            chrom,
            start,
            args.genomic_window
        )


        error_rate = calculate_error_rate(
            fasta,
            chrom,
            start,
            bam
        )


        reads = []

        # iterate over the reads in the given locus
        for read in bam.fetch(
            chrom,
            start,
            stop
        ):

            if read.is_unmapped:
                continue

            if read.reference_end is None:
                continue

            tensor = extract_read_tensor(
                start,
                read
            )


            reads.append(tensor)


        if len(reads) > 0:
            reads = np.array(reads, dtype=object)
        else:
            reads = np.array([], dtype=object)


        outfile = os.path.join(
            rank_dir,
            f"{chrom}_{start}.npz"
        )


        np.savez_compressed(
            outfile,
            chrom=np.array(chrom),
            pos=start,
            context=context,
            reads=reads,
            error_rate=error_rate
        )


    bam.close()
    fasta.close()


if __name__ == "__main__":
    main()