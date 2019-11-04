#!/bin/bash
#SBATCH --time=30:00
#SBATCH --mem=10G
#SBATCH --output=./error_files/monocle_markers.out
#SBATCH --error=./error_files/monocle_markers.err
#SBATCH --job-name=monocle_markers

module load R

Rscript ./R_scripts/monocle_markers.R
