#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --mem=40G
#SBATCH --output=monocle_visualization.out
#SBATCH --error=monocle_visualization.err
#SBATCH --job-name=monocle_visualization

module load R

Rscript monocle_visualizations.R
