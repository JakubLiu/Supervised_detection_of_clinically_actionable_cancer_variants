import pysam
import pandas as pd
import numpy as np
from mpi4py import MPI
from collections import defaultdict
import os
import re
import sys
import argparse
WINDOW_SIZE = 100

parser = argparse.ArgumentParser()

parser.add_argument("--bamlist", type=str)
parser.add_argument("--reference_genome", type=str)
parser.add_argument('--loci_list', type = str)
parser.add_argument('--output_file_prefix', type = str, default = 'output')
parser.add_argument('--outdir', type = str, default = 'OUTPUT_DIR')

args = parser.parse_args()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


def max_homo_run(sequence, base):
    matches = re.findall(f"{base}+", sequence)
    if not matches:
        return 0
    return max(len(x) for x in matches)


def gc_fraction(seq):
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq)


def split_chunks(lst, n):
    return np.array_split(lst, n)



if rank == 0:

    print('reading loci list...')
    loci_df = pd.read_csv(args.loci_list)  # read the file with all the loci (doesnt have to be actionable)
    print('reading loci list done.')
    

    # REQUIRED COLUMNS:
    # chrom,pos,ref
    loci = loci_df[["chrom", "pos", "ref"]].values.tolist()   # turn that into a list

    with open(args.bamlist) as f:
        bamlist = [x.strip() for x in f]

    chunks = split_chunks(loci, size)  # use your function to split it into chunks, one chunk per worker

else:
    bamlist = None
    chunks = None


bamlist = comm.bcast(bamlist, root=0)  # broadcast the bamlist to all workers, because the iteration over the bams will be done sequentailly
local_loci = comm.scatter(chunks, root=0)   # the iterations over the loci will be done in parallel, one chunk per worker



# list of filehandles to open (each worker sees all), this is kinda already opening of the files
bam_handles = [
    pysam.AlignmentFile(bam, "rb")
    for bam in bamlist
]

fasta = pysam.FastaFile(args.reference_genome)


results = []


# this loop is done in parallel (loop over the positions in the current bam)
for idx, (chrom, pos, ref) in enumerate(local_loci):

    chrom = str(chrom)
    pos = int(pos)
    ref = ref.upper()

    print(f"[rank {rank}] processing {chrom}:{pos}", flush=True)

    

    # get the start, end indices for the current worker
    start = max(0, pos - WINDOW_SIZE // 2)
    end = pos + WINDOW_SIZE // 2

    context = fasta.fetch(chrom, start, end).upper()

    center = len(context) // 2

    # this is done serially (all worker sees all bams)
    for bam_path, bam in zip(bamlist, bam_handles):

        sample_id = os.path.basename(bam_path)

        # fetch reads overlapping locus
        for read in bam.fetch(chrom, pos - 1, pos):

            # filters
            if read.is_unmapped:
                continue

            if read.is_secondary or read.is_supplementary:
                continue

            if read.query_sequence is None:
                continue

            # map reference position -> read position
            ref_positions = read.get_reference_positions(full_length=True)

            try:
                read_offset = ref_positions.index(pos - 1)
            except ValueError:
                continue

            # skip deletions/skips
            if read_offset is None:
                continue

            # base quality
            base_qual = read.query_qualities[read_offset]

            # mapping quality
            mapq = read.mapping_quality

            # strand
            strand = "-" if read.is_reverse else "+"

            # distances
            dist_to_read_start = read_offset
            dist_to_read_end = read.query_length - read_offset - 1

            results.append({
                "chrom": chrom,
                "pos": pos,
                "sampleID": sample_id,
                "baseQ": base_qual,
                "mapQ": mapq,
                "strand": strand,
                "dist_to_read_start": dist_to_read_start,
                "dist_to_read_end": dist_to_read_end,
                '100bp_context' : context
            })


all_results = comm.gather(results, root=0)

if rank == 0:

    merged = []

    for x in all_results:
        merged.extend(x)

    df = pd.DataFrame(merged)

    os.makedirs(args.outdir, exist_ok=True)

    outfile = os.path.join(
        args.outdir,
        f"{args.output_file_prefix}.csv"
    )

    df.to_csv(outfile, index=False)