#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --mem=40G
#SBATCH --output=./error_files/create_seurat_object.out
#SBATCH --error=./error_files/create_seurat_object.err
#SBATCH --job-name=create_seurat_object

module load R

Rscript ./R_scripts/create_seurat_object.R
