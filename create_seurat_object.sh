#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --mem=40G
#SBATCH --output=create_seurat_object.out
#SBATCH --error=create_seurat_object.err
#SBATCH --jobname=create_seurat_object

cd1col1="$1"
cd1col2="$2"
cd1col3="$3"
cd1col4="$4"
cd1col5="$5"
cd1col6="$6"
cd2col1="$7"
cd2col2="$8"
cd2col3="$9"
cd2col4="$10"
cd2col5="$11"
cd2col6="$12"
cd3col1="$13"
cd3col2="$14"
cd3col3="$15"
cd3col4="$16"
cd3col5="$17"
cd3col6="$18"

Rscript create_seurat_object.R $cd1col1 $cd1col2 $cd1col3 $cd1col4 $cd1col5 $cd1col6
                        $cd2col1 $cd2col2 $cd2col3 $cd2col4 $cd2col5 $cd2col6
                        $cd3col1 $cd3col2 $cd3col3 $cd3col4 $cd3col5 $cd3col6
