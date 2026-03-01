#!/usr/bin/bash
#SBATCH --job-name=parallel_HPO             # Job name
#SBATCH --output=slurm_logs/parallel_HPO_%j.out        # Standard output (%j = job ID)
#SBATCH --error=slurm_logs/parallel_HPO_%j.err         # Standard error
#SBATCH --mail-type=ALL               # Send email on BEGIN, END, FAIL, REQUEUE, etc.
#SBATCH --mail-user=trolek@gmail.com  # Your email address

#SBATCH --nodes=1                      # Number of nodes
#SBATCH --ntasks=1                     # Number of tasks (processes)
#SBATCH --cpus-per-task=4              # CPUs per task
#SBATCH --mem=70G                      # Total memory
#SBATCH --time=64:00:00                # Walltime (HH:MM:SS)

Rscript GridSearch_HPO_parallel_locus1.R
