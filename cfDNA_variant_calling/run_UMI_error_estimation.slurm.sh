#!/usr/bin/bash

#SBATCH --job-name=umi_error
#SBATCH --output=logs/umi_error_%j.out
#SBATCH --error=logs/umi_error_%j.err

#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=24:00:00

#SBATCH --mail-user=
#SBATCH --mail-type=ALL

mpirun -np 32 python3 UMI_error_rate_single_bam_optimized_mpi.py \
    --bam mapped.grouped.sorted.bam \
    --samplename SAMPLE1 \
    --min_reads_per_UMI 5 \
    --min_coverage 20 \
    --chunk_size 10000 \
    --output sample.error_rates.tsv