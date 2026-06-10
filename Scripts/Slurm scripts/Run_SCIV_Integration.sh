#!/usr/bin/bash

#SBATCH --job-name=SCVI_Integration_Treatment
#SBATCH --output=/srv/data/local/samuelg/Output/SCIV/slurm_output/%x_%j.out
#SBATCH --error=/srv/data/local/samuelg/Output/SCIV/Errors/%x_%j.err

#SBATCH --cpus-per-task=4
#SBATCH --mem=300G
#SBATCH --time=18:00:00

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
Rscript Integration/Integration_SCiv.R

echo "Job finished."

