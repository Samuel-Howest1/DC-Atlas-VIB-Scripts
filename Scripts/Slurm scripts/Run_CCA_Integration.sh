#!/usr/bin/bash

#SBATCH --job-name=CCA_Integration_Experiment
#SBATCH --output=/mnt/temp/srv/data/local/samuelg/Output/CCA/slurm_output/%x_%j.out
#SBATCH --error=/mnt/temp/srv/data/local/samuelg/Output/CCA/Errors/%x_%j.err

#SBATCH --cpus-per-task=12
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
Rscript Integration/Integration_CCA.R

echo "Job finished."
