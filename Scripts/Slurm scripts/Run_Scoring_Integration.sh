#!/usr/bin/bash

#SBATCH --job-name=Score_Harmony_Treat_Exp_V2
#SBATCH --output=/srv/data/local/samuelg/Output/Harmomy/slurm_output/%x_%j.out
#SBATCH --error=/srv/data/local/samuelg/Output/Harmomy/Errors/%x_%j.err

#SBATCH --cpus-per-task=6
#SBATCH --mem=100G
#SBATCH --time=8:00:00

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
Rscript Integration/Integration_Scoring.R

echo "Job finished."


