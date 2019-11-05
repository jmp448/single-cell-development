#!/bin/bash
#SBATCH --time=06:00:00
#SBATCH --mem=40G
#SBATCH --output=./error_files/monocle_dimreduc.out
#SBATCH --error=./error_files/monocle_dimreduc.err
#SBATCH --job-name=monocle_dimreduc

module load R

source_files="$1"

Rscript ./R_scripts/monocle_dimreduc.R source_files
