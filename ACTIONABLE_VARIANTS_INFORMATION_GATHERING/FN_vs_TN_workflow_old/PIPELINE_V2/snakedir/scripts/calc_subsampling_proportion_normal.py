#!/usr/bin/python

import numpy as np

samtools_depth_output = np.loadtxt(snakemake.input[0], usecols = 2, dtype = np.uint16)
target_coverage_normal = snakemake.params.rule_target_coverage_normal
median_depth = np.median(samtools_depth_output)

if median_depth < target_coverage_normal:
    raise RuntimeError("Target coverage is higher than the original coverage. Stopping Snakemake.")


subsampling_proportion = np.round(float(target_coverage_normal)/float(median_depth),decimals = 2)


with open(snakemake.output[0], "w") as out:
    out.write(f'median_sample_depth,target_depth,subsampling_proportion\n{median_depth},{target_coverage_normal},{subsampling_proportion}')