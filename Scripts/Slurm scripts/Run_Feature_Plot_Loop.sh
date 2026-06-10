#!/usr/bin/bash

#SBATCH --job-name=Harmony_Feature_Loop_Theta_2
#SBATCH --output=/srv/data/local/samuelg/Output/Harmomy/slurm_output/%x_%j.out
#SBATCH --error=/srv/data/local/samuelg/Output/Harmomy/Errors/%x_%j.err

#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=2:00:00

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
Rscript Integration/Feature_plot_loop_Theta_2.R

echo "Job finished."

