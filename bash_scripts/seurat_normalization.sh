#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --mem=40G
#SBATCH --output=./error_files/seurat_normalization.out
#SBATCH --error=./error_files/seurat_normalization.err
#SBATCH --job-name=seurat_normalization

source_files="${1}"

module load R python/3.5.2

Rscript ./R_scripts/seurat_normalization.R $source_files
