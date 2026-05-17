#!/usr/bin/python

"""
==================================== EXAMPLE HOW TO RUN =======================================

mpirun -n 32 python3 make_data_optimized.py \
    --bamlist "/data/cephfs-1/home/users/jali13_c/work/EPSILON/MASK/make_data/single_bamlist.txt" \
    --reference_genome "/data/cephfs-1/home/users/jali13_c/work/BAYESIAN_PRIOR/TUNE_PRIOR/STEP1_PREPROCESS_BAMS/hs37d5.fa" \
    --loci_list "loci_list_minimal.csv" \
    --outdir "testodir"


"""



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

    loci_df = pd.read_csv(args.loci_list)  # read the file with all the loci (doesnt have to be actionable)

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

# this loop is done in parallel
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

    GGC_upstream = int("GGC" in context[:center])
    CGG_downstream = int("CGG" in context[center:])
    GGT_upstream = int("GGT" in context[:center])
    TGG_downstream = int("TGG" in context[center:])

    GC_content = gc_fraction(context)

    

    # this is done serially (all worker sees all bams)
    for bam_path, bam in zip(bamlist, bam_handles):

        sample_id = os.path.basename(bam_path)

        # set the variables for the current bam
        ref_count = 0
        alt_count = 0

        mapqs = []   # this is a list of all mapqs of all reads that have mapped to the current pos in the current bam
        baseqs = []   # this is a list of all baseqs at positions of reads that have mapped to the current pos in the current bam

        mismatch_values = [] # a list for the current position for the current bam

        ref_R1 = 0
        ref_R2 = 0
        alt_R1 = 0
        alt_R2 = 0


        for read in bam.fetch(chrom, pos - 1, pos):


            if read.is_unmapped:
                continue

            if read.is_duplicate:
                continue

            if read.mapping_quality == 0:
                continue

            aligned_pairs = read.get_aligned_pairs(matches_only=False)

            query_pos = None

            for qpos, rpos in aligned_pairs:

                if rpos == pos - 1:
                    query_pos = qpos
                    break

            if query_pos is None:
                continue

            if query_pos >= len(read.query_sequence):
                continue

            base = read.query_sequence[query_pos].upper()

            if base == "N":
                continue

            baseq = read.query_qualities[query_pos]
            mapq = read.mapping_quality

            mapqs.append(mapq)
            baseqs.append(baseq)

            try:
                nm = read.get_tag("NM")
            except KeyError:
                nm = 0

            mismatch_values.append(nm)

            strand = "-" if read.is_reverse else "+"

            

            if base == ref:

                ref_count += 1

                if strand == "+":
                    ref_R1 += 1
                else:
                    ref_R2 += 1

            else:

                alt_count += 1

                if strand == "+":
                    alt_R1 += 1
                else:
                    alt_R2 += 1

      

        coverage = ref_count + alt_count

        if coverage == 0:
            continue

        # for the given position for the given bam compute the medians/means of the lists
        median_mapQ = np.median(mapqs) if mapqs else np.nan
        median_baseQ = np.median(baseqs) if baseqs else np.nan

        mismatches = np.mean(
            np.array(mismatch_values) > 0
        )

       

        row = {
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "sampleID": sample_id,

            "general_alt_counts": alt_count,
            "ref_counts": ref_count,

            "coverage": coverage,

            "mismatches": mismatches,

            "median_mapQ": median_mapQ,
            "median_baseQ": median_baseQ,

            "GGC_upstream": GGC_upstream,
            "CGG_downstream": CGG_downstream,
            "GGT_upstream": GGT_upstream,
            "TGG_downstream": TGG_downstream,

            "GC_content": GC_content,

            "ref_R1": ref_R1,
            "ref_R2": ref_R2,
            "alt_R1": alt_R1,
            "alt_R2": alt_R2
        }

        results.append(row)


os.makedirs(args.outdir, exist_ok=True)

local_df = pd.DataFrame(results)

outfile = os.path.join(
    args.outdir,
    f"{args.output_file_prefix}.rank{rank}.txt"
)

local_df.to_csv(outfile)

print(f"[rank {rank}] wrote {outfile}", flush=True)


for bam in bam_handles:
    bam.close()

fasta.close()

comm.Barrier()

if rank == 0:
    print("All ranks completed.", flush=True)
