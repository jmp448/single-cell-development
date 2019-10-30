#!/bin/bash
#SBATCH --time=06:00:00
#SBATCH --mem=40G
#SBATCH --output=./error_files/seurat_dimreduc.out
#SBATCH --error=./error_files/seurat_dimreduc.err
#SBATCH --job-name=seurat_dimreduc

module load R

Rscript ./R_scripts/seurat_dimreduc.R
