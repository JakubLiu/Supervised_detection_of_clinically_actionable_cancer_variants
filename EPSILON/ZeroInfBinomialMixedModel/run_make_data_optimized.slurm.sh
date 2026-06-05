#!/usr/bin/bash
#SBATCH --job-name=mpi_make_data
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=100:00:00
#SBATCH --mem=80G
#SBATCH --mail-user=
#SBATCH --mail-type=ALL


mpirun -n 32 python3 make_data_optimized_mpi.py \
    --bamlist "bamlist.txt" \
    --reference_genome "hs37d5.fa" \
    --loci_list "Civic_actionable_SNV_loci.csv" \
    --outdir "data_dir"

