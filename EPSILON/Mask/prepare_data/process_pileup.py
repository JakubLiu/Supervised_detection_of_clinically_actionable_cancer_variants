#!/usr/bin/python3

import sys
import re

def parse_bases(bases):
    """
    Parse mpileup bases string.

    Returns:
        ref_count: number of reference observations
        alt_count: number of non-reference observations
    """

    ref_count = 0
    alt_count = 0

    i = 0
    n = len(bases)

    while i < n:
        c = bases[i]

        # start of read marker: ^ plus following mapping quality
        if c == "^":
            i += 2
            continue

        # end of read marker
        if c == "$":
            i += 1
            continue

        # reference matches
        if c in ".,":  # forward/reverse strand reference
            ref_count += 1
            i += 1
            continue

        # insertion/deletion
        if c in "+-":
            i += 1

            num = ""
            while i < n and bases[i].isdigit():
                num += bases[i]
                i += 1

            if num:
                indel_len = int(num)
                i += indel_len

            continue

        # deletion placeholder from previous line
        if c == "*":
            alt_count += 1
            i += 1
            continue

        # A/C/G/T/N (either strand)
        if c.upper() in "ACGTN":
            alt_count += 1
            i += 1
            continue

        i += 1

    return ref_count, alt_count


def main():
    print("chrom\tpos\tref\tcoverage\tnon_ref_reads\terror_rate")

    for line in sys.stdin:
        line = line.rstrip()
        if not line:
            continue

        fields = line.split("\t")

        if len(fields) < 5:
            continue

        chrom = fields[0]
        pos = fields[1]
        ref = fields[2]
        depth = int(fields[3])
        bases = fields[4]

        ref_count, alt_count = parse_bases(bases)

        coverage = ref_count + alt_count

        # use parsed coverage if it differs from reported depth
        if coverage == 0:
            error_rate = 0.0
        else:
            error_rate = alt_count / coverage

        print(
            f"{chrom}\t{pos}\t{ref}\t{coverage}\t{alt_count}\t{error_rate:.6f}"
        )


if __name__ == "__main__":
    main()
