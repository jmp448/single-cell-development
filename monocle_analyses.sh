#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --mem=40G
#SBATCH --output=monocle_analyses.out
#SBATCH --error=monocle_analyses.err
#SBATCH --job-name=monocle_analyses

module load R

Rscript monocle_analyses.R
