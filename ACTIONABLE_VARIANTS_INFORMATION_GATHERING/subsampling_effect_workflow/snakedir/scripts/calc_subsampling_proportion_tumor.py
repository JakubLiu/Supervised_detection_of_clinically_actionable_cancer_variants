#!/usr/bin/python

import numpy as np

median_depths_all_tumor_samples = []
target_coverage_tumor = snakemake.params.rule_target_coverage_tumor

for tumor_sample in snakemake.input:
    median_depths_all_tumor_samples.append(np.loadtxt(tumor_sample, dtype = np.float64))

mean_depth = np.mean(median_depths_all_tumor_samples) # mean of the per sample median depths


if np.max(median_depths_all_tumor_samples) < target_coverage_tumor:
    raise RuntimeError("Target coverage is higher than the original coverage. Stopping Snakemake.")


subsampling_proportion = np.round(float(target_coverage_tumor)/float(mean_depth),decimals = 2)


with open(snakemake.output[0], "w") as out:
    out.write(f'median_sample_depth,target_depth,subsampling_proportion\n{mean_depth},{target_coverage_tumor},{subsampling_proportion}')
