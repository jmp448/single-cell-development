#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --mem=40G
#SBATCH --output=./error_files/monocle_analyses.out
#SBATCH --error=./error_files/monocle_analyses.err
#SBATCH --job-name=monocle_analyses

module load R

Rscript monocle_analyses.R
