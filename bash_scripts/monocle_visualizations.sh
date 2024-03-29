#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --mem=40G
#SBATCH --output=./error_files/monocle_visualizations.out
#SBATCH --error=./error_files/monocle_visualizations.err
#SBATCH --job-name=monocle_visualization

module load R

Rscript ./R_scripts/monocle_visualizations.R
