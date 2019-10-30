#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --mem=40G
#SBATCH --output=./error_files/seurat2monocle.out
#SBATCH --error=./error_files/seurat2monocle.err
#SBATCH --job-name=seurat2monocle

Rscript ./R_scripts/seurat2monocle.R
