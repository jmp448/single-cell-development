#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --mem=30G
#SBATCH --output=./error_files/doublet_data.out
#SBATCH --error=./error_files/doublet_data.err
#SBATCH --job-name=doublet_data

module load R

Rscript ./R_scripts/doublet_dataframes.R
