##### INPUT FILES #####

source_files="lowpass"

# LOW PASS DATA
# If you want to use the low pass data, set this to true, otherwise set to false
if [ "$source_files" == "lowpass" ]; then
   cd1col1="/project2/gilad/reem/singlecellCM/round1/lowpass/CD1/CD1col1/output/dge_data/YG-RE1-Drop-CD1col1_S1_gene_counts.tsv.gz"
   cd1col2="/project2/gilad/reem/singlecellCM/round1/lowpass/CD1/CD1col2/output/dge_data/YG-RE1-Drop-CD1col2_S2_gene_counts.tsv.gz"
   cd1col3="/project2/gilad/reem/singlecellCM/round1/lowpass/CD1/CD1col3/output/dge_data/YG-RE1-Drop-CD1col3_S3_gene_counts.tsv.gz"
   cd1col4="/project2/gilad/reem/singlecellCM/round1/lowpass/CD1/CD1col4/output/dge_data/YG-RE2-Drop-CD1col4_S1_gene_counts.tsv.gz"
   cd1col5="/project2/gilad/reem/singlecellCM/round1/lowpass/CD1/CD1col5/output/dge_data/YG-RE2-Drop-CD1col5_S2_gene_counts.tsv.gz"
   cd1col6="/project2/gilad/reem/singlecellCM/round1/lowpass/CD1/CD1col6/output/dge_data/YG-RE2-Drop-CD1col6_S3_gene_counts.tsv.gz"
   cd2col1="/project2/gilad/reem/singlecellCM/round1/lowpass/CD2/CD2col1/output/dge_data/YG-RE1-Drop-CD2col1_S4_gene_counts.tsv.gz"
   cd2col2="/project2/gilad/reem/singlecellCM/round1/lowpass/CD2/CD2col2/output/dge_data/YG-RE2-Drop-CD2col2_S4_gene_counts.tsv.gz"
   cd2col3="/project2/gilad/reem/singlecellCM/round1/lowpass/CD2/CD2col3/output/dge_data/YG-RE2-Drop-CD2col3_S5_gene_counts.tsv.gz"
   cd2col4="/project2/gilad/reem/singlecellCM/round1/lowpass/CD2/CD2col4/output/dge_data/YG-RE1-Drop-CD2col4_S5_gene_counts.tsv.gz"
   cd2col5="/project2/gilad/reem/singlecellCM/round1/lowpass/CD2/CD2col5/output/dge_data/YG-RE2-Drop-CD2col5_S6_gene_counts.tsv.gz"
   cd2col6="/project2/gilad/reem/singlecellCM/round1/lowpass/CD2/CD2col6/output/dge_data/YG-RE1-Drop-CD2col6_S6_gene_counts.tsv.gz"
   cd3col1="/project2/gilad/reem/singlecellCM/round1/lowpass/CD3/CD3col1/output/dge_data/YG-RE2-Drop-CD3col1_S7_gene_counts.tsv.gz"
   cd3col2="/project2/gilad/reem/singlecellCM/round1/lowpass/CD3/CD3col2/output/dge_data/YG-RE2-Drop-CD3col2_S8_gene_counts.tsv.gz"
   cd3col3="/project2/gilad/reem/singlecellCM/round1/lowpass/CD3/CD3col3/output/dge_data/YG-RE1-Drop-CD3col3_S7_gene_counts.tsv.gz"
   cd3col4="/project2/gilad/reem/singlecellCM/round1/lowpass/CD3/CD3col4/output/dge_data/YG-RE1-Drop-CD3col4_S8_gene_counts.tsv.gz"
   cd3col5="/project2/gilad/reem/singlecellCM/round1/lowpass/CD3/CD3col5/output/dge_data/YG-RE1-Drop-CD3col5_S9_gene_counts.tsv.gz"
   cd3col6="/project2/gilad/reem/singlecellCM/round1/lowpass/CD3/CD3col6/output/dge_data/YG-RE2-Drop-CD3col6_S9_gene_counts.tsv.gz"

# FULL DATA
# If you want to use the full data, set this to true, otherwise set to false
elif [ "$source_files" == "fulldata" ]; then
   cd1col1='/project2/gilad/reem/singlecellCM/round1/fulldata/CD1/CD1col1/output/dge_data/YG-RE-RE1-hpCD1col1_S1_gene_counts.tsv.gz'
   cd1col2='/project2/gilad/reem/singlecellCM/round1/fulldata/CD1/CD1col2/output/dge_data/YG-RE-RE2-hpCD1col2_S1_gene_counts.tsv.gz'
   cd1col3='/project2/gilad/reem/singlecellCM/round1/fulldata/CD1/CD1col3/output/dge_data/YG-RE-RE3-hpCD1col3_S1_gene_counts.tsv.gz'
   cd1col4='/project2/gilad/reem/singlecellCM/round1/fulldata/CD1/CD1col4/output/dge_data/YG-RE-RE4-hpCD1col4_S1_gene_counts.tsv.gz'
   cd1col5='/project2/gilad/reem/singlecellCM/round1/fulldata/CD1/CD1col5/output/dge_data/YG-RE-RE5-hpCD1col5_S1_gene_counts.tsv.gz'
   cd1col6='/project2/gilad/reem/singlecellCM/round1/fulldata/CD1/CD1col6/output/dge_data/YG-RE-RE6-hpCD1col6_S1_gene_counts.tsv.gz'
   cd2col1='/project2/gilad/reem/singlecellCM/round1/fulldata/CD2/CD2col1/output/dge_data/YG-RE-RE3-hpCD2col1_S2_gene_counts.tsv.gz'
   cd2col2='/project2/gilad/reem/singlecellCM/round1/fulldata/CD2/CD2col2/output/dge_data/YG-RE-RE4-hpCD2col2_S2_gene_counts.tsv.gz'
   cd2col3='/project2/gilad/reem/singlecellCM/round1/fulldata/CD2/CD2col3/output/dge_data/YG-RE-RE5-hpCD2col3_S2_gene_counts.tsv.gz'
   cd2col4='/project2/gilad/reem/singlecellCM/round1/fulldata/CD2/CD2col4/output/dge_data/YG-RE-RE6-hpCD2col4_S2_gene_counts.tsv.gz'
   cd2col5='/project2/gilad/reem/singlecellCM/round1/fulldata/CD2/CD2col5/output/dge_data/YG-RE-RE2-hpCD2col5_S2_gene_counts.tsv.gz'
   cd2col6='/project2/gilad/reem/singlecellCM/round1/fulldata/CD2/CD2col6/output/dge_data/YG-RE-RE1-hpCD2col6_S2_gene_counts.tsv.gz'
   cd3col1='/project2/gilad/reem/singlecellCM/round1/fulldata/CD3/CD3col1/output/dge_data/YG-RE-RE4-hpCD3col1_S3_gene_counts.tsv.gz'
   cd3col2='/project2/gilad/reem/singlecellCM/round1/fulldata/CD3/CD3col2/output/dge_data/YG-RE-RE3-hpCD3col2_S3_gene_counts.tsv.gz'
   cd3col3='/project2/gilad/reem/singlecellCM/round1/fulldata/CD3/CD3col3/output/dge_data/YG-RE-RE2-CD3col3_Unk1_S4_gene_counts.tsv.gz'
   cd3col4='/project2/gilad/reem/singlecellCM/round1/fulldata/CD3/CD3col4/output/dge_data/YG-RE-RE1-hpCD3col4_S3_gene_counts.tsv.gz'
   cd3col5='/project2/gilad/reem/singlecellCM/round1/fulldata/CD3/CD3col5/output/dge_data/YG-RE-RE6-hpCD3col5_S3_gene_counts.tsv.gz'
   cd3col6='/project2/gilad/reem/singlecellCM/round1/fulldata/CD3/CD3col6/output/dge_data/YG-RE-RE5-hpCD3col6_S3_gene_counts.tsv.gz'
fi

# PART 1: Create Seurat Object from raw data files
# Takes around an hour and a half to run on full data
min_cells_per_gene=3  # minimum num cells in which a gene must appear
min_genes_per_cell=200  # minimum num genes for a cell to be included
#TODO incorporate mito threshold boolean into pipeline
mito_threshold=30  # what is the percentage threshold cutoff for mito counts?

if false; then
   rm ./error_files/create_seurat_object.out
   rm ./error_files/create_seurat_object.err
   sbatch ./bash_scripts/create_seurat_object.sh $cd1col1 $cd1col2 $cd1col3 $cd1col4 $cd1col5 $cd1col6 \
                                    $cd2col1 $cd2col2 $cd2col3 $cd2col4 $cd2col5 $cd2col6 \
                                    $cd3col1 $cd3col2 $cd3col3 $cd3col4 $cd3col5 $cd3col6 \
                                    $min_cells_per_gene $min_genes_per_cell $mito_threshold \
                                    $source_files
fi

# PART 2: Normalize Seurat data
if false; then
   rm ./error_files/seurat_normalization.err
   rm ./error_files/seurat_normalization.out
   sbatch ./bash_scripts/seurat_normalization.sh $source_files
fi

# PART 3: Prepare monocle objects from seurat objects
if true; then
   rm ./error_files/seurat2monocle.err
   rm ./error_files/seurat2monocle.out
   sbatch ./bash_scripts/seurat2monocle.sh
fi

if false; then
   rm ./error_files/monocle_dimreduc.err
   rm ./error_files/monocle_dimreduc.out
   sbatch ./bash_scripts/monocle_dimreduc.sh
fi

if false; then
   rm ./error_files/seurat_dimreduc.err
   rm ./error_files/seurat_dimreduc.out
   sbatch ./bash_scripts/seurat_dimreduc.sh
fi

if false; then
   rm ./error_files/seurat_visualizations.err
   rm ./error_files/seurat_visualizations.out
   sbatch ./bash_scripts/seurat_visualizations.sh
fi

if false; then
   rm ./error_files/monocle_visualizations.err
   rm ./error_files/monocle_visualizations.out
   sbatch ./bash_scripts/monocle_visualizations.sh
fi

# Monocle Analysis
if false; then
   # bash install_monocle.sh
   sbatch ./bash_scripts/create_monocle_object.sh $cd1col1 $cd1col2 $cd1col3 $cd1col4 $cd1col5 $cd1col6 \
                                    $cd2col1 $cd2col2 $cd2col3 $cd2col4 $cd2col5 $cd2col6 \
                                    $cd3col1 $cd3col2 $cd3col3 $cd3col4 $cd3col5 $cd3col6 \
                                    $min_cells_per_gene $min_genes_per_cell $mito_threshold \
                                    $source_files
fi
