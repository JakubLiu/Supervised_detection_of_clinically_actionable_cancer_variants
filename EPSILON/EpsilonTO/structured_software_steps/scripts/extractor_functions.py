#!/usr/bin/python
import pysam
import numpy as np
import argparse
import os
from itertools import product


def parse_args():

    parser = argparse.ArgumentParser(
        description="Extract BAM read tensors for model training"
    )

    parser.add_argument(
            "--reference_genome",
            required=True,
            type = str,
            help="Path to the reference genome"
        )

    parser.add_argument(
        "--bam",
        required=True,
        type = str,
        help="Path to the bam file"
    )

    parser.add_argument(
        "--bed",
        required=True,
        type = str,
        help="BED TSV file: chr<TAB>start<TAB>stop<TAB>ref"
    )

    parser.add_argument(
        "--output",
        required=True,
        type = str,
        help="Output directory"
    )

    parser.add_argument(
        "--genomic_window",
        type=int,
        default=1001,
        help="The size of the genomic window"
    )


    return parser.parse_args()



# reads the bam files
def read_bams(path):
    with open(path) as f:
        return [x.strip() for x in f if x.strip()]


# reads the bed file
def read_bed(path):

    regions = []

    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, start, stop, ref= line.rstrip().split("\t")[:5]

            regions.append((chrom,int(start),int(stop),ref.upper(),))

    return regions



# for the given locus and the given bam extract the read 2D tensor
'''
(locus + BAM + read) --> (d,9) tensor of that read ---> CNN
d is the variable read length (this varies, because different reads can have different lengths)
'''
def extract_read_tensor(
        start,  # genomic position from BED
        read    # current pysam read
):

    sequence = read.query_sequence
    qualities = read.query_qualities

    # handle reads without sequence
    if sequence is None or len(sequence) == 0:
        return np.zeros((0, 9), dtype=np.float32)

    distance_to_read_end = read.reference_end - start

    read_tensor = []

    # iterate over every base in the read
    for i, base in enumerate(sequence):

        feature = np.zeros(9, dtype=np.float32)

        # base one-hot encoding
        if base == "A":
            feature[0] = 1.0
        elif base == "C":
            feature[1] = 1.0
        elif base == "G":
            feature[2] = 1.0
        elif base == "T":
            feature[3] = 1.0
        # N remains [0,0,0,0]

        # read-level features
        feature[5] = read.mapping_quality          # mapping quality
        feature[6] = qualities[i]                  # base quality
        feature[7] = 1.0 if read.is_reverse else 0.0  # strand
        feature[8] = distance_to_read_end / len(sequence)  # normalized distance

        read_tensor.append(feature)

    read_tensor = np.array(read_tensor, dtype=np.float32)
    return read_tensor






    
# for the given locus extract the genomic context
'''
locus ---> context window ---> MHSA + positional encoding
'''
def extract_reference_context(fasta,
                              chrom,
                              pos,   # ---> pos is the first bed entry 100\tab101 --> pos = 100
                              window):

    center = pos
    half_window = window // 2
    region_start = center - half_window
    region_end = center + half_window

    # extract sequence
    seq = fasta.fetch(chrom, region_start, region_end).upper()

    # one-hot encode
    encoded = np.zeros((window, 4), dtype=np.float32)

    for i, base in enumerate(seq):

        if base == "A":
            encoded[i,0] = 1.0

        elif base == "C":
            encoded[i,1] = 1.0

        elif base == "G":
            encoded[i,2] = 1.0

        elif base == "T":
            encoded[i,3] = 1.0

    return encoded


def calculate_error_rate(refgen, # should aready be opened with:  refgen = pysam.FastaFile("reference.fa")
                                 chrom,
                                 pos,   # 0 based, so the first column in the bed, if the bed is 100\tab101, then pos = 100
                                 bam   # should already be opened with:  bam = pysam.AlignmentFile("sample.bam", "rb")
                                 ):


    ref_base = refgen.fetch(chrom, pos, pos + 1).upper()

    matches = 0
    mismatches = 0

    for pileupcolumn in bam.pileup(
            chrom,
            pos,
            pos + 1,
            truncate=True,
            stepper="all",
            min_base_quality=20):

        if pileupcolumn.reference_pos != pos:
            continue

        for pileupread in pileupcolumn.pileups:

            if pileupread.is_del or pileupread.is_refskip:
                continue

            read = pileupread.alignment
            qpos = pileupread.query_position
            base = read.query_sequence[qpos].upper()

            if base == ref_base:
                matches += 1
            else:
                mismatches += 1

    coverage = matches + mismatches
    error_rate = mismatches / coverage if coverage else 0
    return error_rate
