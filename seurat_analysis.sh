##### INPUT FILES #####

# LOW PASS DATA


# FULL DATA
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

# PART 1: Create Seurat Object from raw data files
if true; then
   sbatch create_seurat_object.sh $cd1col1 $cd1col2 $cd1col3 $cd1col4 $cd1col5 $cd1col6
                                    $cd2col1 $cd2col2 $cd2col3 $cd2col4 $cd2col5 $cd2col6
                                    $cd3col1 $cd3col2 $cd3col3 $cd3col4 $cd3col5 $cd3col6
fi
