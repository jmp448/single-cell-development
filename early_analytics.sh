#!/bin/bash
#SBATCH --time=20:00
#SBATCH --mem=5G
#SBATCH --output=early_analytics.out
#SBATCH --error=early_analytics.err
#SBATCH --job-name=early_analytics

module load R

Rscript early_analytics.R
