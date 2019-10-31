#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --mem=40G
#SBATCH --output=./error_files/create_seurat_object.out
#SBATCH --error=./error_files/create_seurat_object.err
#SBATCH --job-name=create_seurat_object

module load R

cd1col1="$0"
cd1col2="$1"
cd1col3="$2"
cd1col4="$3"
cd1col5="$4"
cd1col6="$5"
cd2col1="$6"
cd2col2="$7"
cd2col3="$8"
cd2col4="$9"
cd2col5="${10}"
cd2col6="${11}"
cd3col1="${12}"
cd3col2="${13}"
cd3col3="${14}"
cd3col4="${15}"
cd3col5="${16}"
cd3col6="${17}"

min_cells_per_gene=${18}  # minimum num cells in which a gene must appear
min_genes_per_cell=${19}  # minimum num genes for a cell to be included
mito_threshold=${20}  # what is the threshold for mito cutoff?
source_files="${21}"

Rscript ./R_scripts/create_seurat_object.R $cd1col1 $cd1col2 $cd1col3 $cd1col4 $cd1col5 $cd1col6 \
                        $cd2col1 $cd2col2 $cd2col3 $cd2col4 $cd2col5 $cd2col6 \
                        $cd3col1 $cd3col2 $cd3col3 $cd3col4 $cd3col5 $cd3col6 \
                        $min_cells_per_gene $min_genes_per_cell $mito_threshold \
                        $source_files
