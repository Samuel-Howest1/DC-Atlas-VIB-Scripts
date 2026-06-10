#!/usr/bin/bash

#SBATCH --job-name=Harmony_Integration_Treat_Exp_Theta_ALl_3
#SBATCH --output=/srv/data/local/samuelg/Output/Harmomy/slurm_output/%x_%j.out
#SBATCH --error=/srv/data/local/samuelg/Output/Harmomy/Errors/%x_%j.err

#SBATCH --cpus-per-task=12
#SBATCH --mem=300G
#SBATCH --time=12:00:00

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=samuel.gundes@irc.vib-ugent.be

set -euo pipefail

# Load R module
module purge


# Print some diagnostics
echo "Job started on:"
hostname

echo "Current directory:"
pwd

echo "Starting R script..."

# Run R
Rscript Integration/Integration_Harmony.R

echo "Job finished."
