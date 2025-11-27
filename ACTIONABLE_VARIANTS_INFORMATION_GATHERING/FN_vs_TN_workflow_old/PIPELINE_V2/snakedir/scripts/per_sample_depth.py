#!/usr/bin/python

import numpy as np

samtools_depth_output = np.loadtxt(snakemake.input[0], usecols = 2, dtype = np.uint16)
median_depth = np.median(samtools_depth_output)

with open(snakemake.output[0], "w") as out:
    out.write(f"{median_depth}\n")