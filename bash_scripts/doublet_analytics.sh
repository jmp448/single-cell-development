#!/bin/bash
#SBATCH --time=20:00
#SBATCH --mem=10G
#SBATCH --output=./error_files/doublet_analytics.out
#SBATCH --error=./error_files/doublet_analytics.err
#SBATCH --job-name=doublet_analytics

module load R

Rscript doublet_analytics.R
