#!/usr/bin/bash
#SBATCH --job-name=entropy_defaultname
#SBATCH --output=entropy_defaultname_%j.out
#SBATCH --error=entropy_defualtname_%j.err
#SBATCH --ntasks=10
#SBATCH --time=300:00:00
#SBATCH --mem=80G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=******

REF="hs37d5.fa"
CHROM=$1     # variable
WINDOW="100000"
OUTPUT_FILE=$2   # variable


mpirun -np 10 python3 Calc_Shannon_entropy.mpi.numba.py "$REF" \
                                        "$CHROM" \
                                        "$WINDOW" \
                                        "$OUTPUT_FILE"


# How to run:
# sbatch --job-name="entropy_${1}" run_calc_entropy.slurm.sh 1 chr1_entropy.txt
