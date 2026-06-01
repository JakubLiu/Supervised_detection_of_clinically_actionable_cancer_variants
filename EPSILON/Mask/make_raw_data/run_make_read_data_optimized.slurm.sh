#!/bin/bash
#SBATCH --job-name=mask
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem=80G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=




mpirun -n 64 python3 make_read_data_mpi.py \
    --bamlist "bamlist.txt" \
    --reference_genome "hs37d5.fa" \
    --loci_list "output_chr1.halve.txt" \
    --outdir "results"
