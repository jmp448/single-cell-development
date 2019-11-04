#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=40G
#SBATCH --output=./error_files/all_in_one.out
#SBATCH --error=./error_files/all_in_one.err
#SBATCH --job-name=sc_pipeline

cd1col1="$1"
cd1col2="$2"
cd1col3="$3"
cd1col4="$4"
cd1col5="$5"
cd1col6="$6"
cd2col1="$7"
cd2col2="$8"
cd2col3="$9"
cd2col4="${10}"
cd2col5="${11}"
cd2col6="${12}"
cd3col1="${13}"
cd3col2="${14}"
cd3col3="${15}"
cd3col4="${16}"
cd3col5="${17}"
cd3col6="${18}"

min_cells_per_gene=${19}  # minimum num cells in which a gene must appear
min_genes_per_cell=${20}  # minimum num genes for a cell to be included
mito_threshold=${21}  # what is the threshold for mito cutoff?
source_files="${22}"

module load R

# srun Rscript ./R_scripts/create_seurat_object.R $cd1col1 $cd1col2 $cd1col3 $cd1col4 $cd1col5 $cd1col6 \
#                         $cd2col1 $cd2col2 $cd2col3 $cd2col4 $cd2col5 $cd2col6 \
#                         $cd3col1 $cd3col2 $cd3col3 $cd3col4 $cd3col5 $cd3col6 \
#                         $min_cells_per_gene $min_genes_per_cell $mito_threshold \
#                         $source_files 0
# wait
srun Rscript ./R_scripts/seurat_normalization.R $source_files 1
wait
srun Rscript ./R_scripts/seurat2monocle.R 2
wait
srun Rscript ./R_scripts/seurat_dimreduc.R  3 &
srun Rscript ./R_scripts/monocle_dimreduc.R 4
wait
# srun Rscript ./R_scripts/seurat_visualizations.R 5 &
# srun Rscript ./R_scripts/monocle_visualizations.R 6
# wait
