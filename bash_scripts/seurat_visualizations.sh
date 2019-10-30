#!/bin/bash
#SBATCH --time=20:00
#SBATCH --output=./error_files/seurat_visualizations.out
#SBATCH --error=./error_files/seurat_visualizations.err
#SBATCH --job-name=seurat_visualizations

module load R

Rscript ./R_scripts/seurat_visualizations.R
