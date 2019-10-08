#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --partition=bigmem2
#SBATCH --mem=100G
#SBATCH --error=sctransform_error.err
#SBATCH --output=sctransform_output.out
#SBATCH --job-name=sctransform_proper


module load rstudio
module load python/3.5.2
pip install umap-learn

Rscript sctransform_pca_umap.R
