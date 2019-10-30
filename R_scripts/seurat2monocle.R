library(Seurat)
library(monocle3)

# Create monocle object from just the straight seurat data (will need to be preprocessed in monocle)
seurat_file <- "./rds_objects/seurat_obj_lowpass.RDS"
sc_seurat <- readRDS(seurat_file)
expression_matrix <- sc_seurat@assays[["RNA"]]@counts
cell_metadata <- sc_seurat[[]]
gene_metadata <- as.data.frame(rownames(expression_matrix))
rownames(gene_metadata) <- rownames(expression_matrix)
colnames(gene_metadata) <- "gene_short_name"
monocle <- new_cell_data_set(expression_data = expression_matrix,
                             gene_metadata = gene_metadata,
                             cell_metadata = cell_metadata)
monocle <- monocle[,Matrix::colSums(exprs(monocle)) != 0]
monocle <- estimate_size_factors(monocle)
saveRDS(monocle, "./rds_objects/monocle_from_seurat_lowpass_unprocessed.RDS")
rm(list=ls())

# Create monocle object from the log normalized seurat data
seurat_file <- "./rds_objects/seurat_obj_lowpass_basic_norm.RDS"
sc_seurat <- readRDS(seurat_file)
expression_matrix <- sc_seurat@assays[["RNA"]]@scale.data
cell_metadata <- sc_seurat[[]]
gene_metadata <- as.data.frame(rownames(expression_matrix))
rownames(gene_metadata) <- rownames(expression_matrix)
colnames(gene_metadata) <- "gene_short_name"
monocle <- new_cell_data_set(expression_data = expression_matrix,
                             gene_metadata = gene_metadata,
                             cell_metadata = cell_metadata)
monocle <- monocle[,Matrix::colSums(exprs(monocle)) != 0]
saveRDS(monocle, "./rds_objects/monocle_from_seurat_lowpass_basic_norm.RDS")
rm(list=ls())

# Create monocle object from the log normalized seurat data but with no scale factor
seurat_file <- "./rds_objects/seurat_obj_lowpass_noscale_norm.RDS"
sc_seurat <- readRDS(seurat_file)
expression_matrix <- sc_seurat@assays[["RNA"]]@scale.data
cell_metadata <- sc_seurat[[]]
gene_metadata <- as.data.frame(rownames(expression_matrix))
rownames(gene_metadata) <- rownames(expression_matrix)
colnames(gene_metadata) <- "gene_short_name"
monocle <- new_cell_data_set(expression_data = expression_matrix,
                             gene_metadata = gene_metadata,
                             cell_metadata = cell_metadata)
monocle <- monocle[,Matrix::colSums(exprs(monocle)) != 0]
saveRDS(monocle, "./rds_objects/monocle_from_seurat_lowpass_noscale_norm.RDS")
rm(list=ls())

# Create monocle object from SC transformed seurat object
seurat_file <- "./rds_objects/seurat_obj_lowpass_sctransform.RDS"
sc_seurat <- readRDS(seurat_file)
expression_matrix <- sc_seurat@assays[["SCT"]]@scale.data
cell_metadata <- sc_seurat[[]]
gene_metadata <- as.data.frame(rownames(expression_matrix))
rownames(gene_metadata) <- rownames(expression_matrix)
colnames(gene_metadata) <- "gene_short_name"
monocle <- new_cell_data_set(expression_data = expression_matrix,
                             gene_metadata = gene_metadata,
                             cell_metadata = cell_metadata)
monocle <- monocle[,Matrix::colSums(exprs(monocle)) != 0]
saveRDS(monocle, "./rds_objects/monocle_from_seurat_lowpass_sctransform.RDS")
rm(list=ls())

# Create monocle object from SC transformed seurat object with mito regressed out
seurat_file <- "./rds_objects/seurat_obj_lowpass_sctransform_mitoregress.RDS"
sc_seurat <- readRDS(seurat_file)
expression_matrix <- sc_seurat@assays[["RNA"]]@scale.data
cell_metadata <- sc_seurat[[]]
gene_metadata <- as.data.frame(rownames(expression_matrix))
rownames(gene_metadata) <- rownames(expression_matrix)
colnames(gene_metadata) <- "gene_short_name"
monocle <- new_cell_data_set(expression_data = expression_matrix,
                             gene_metadata = gene_metadata,
                             cell_metadata = cell_metadata)
monocle <- monocle[,Matrix::colSums(exprs(monocle)) != 0]
saveRDS(monocle, "./rds_objects/monocle_from_seurat_lowpass_sctransform_mitoregress.RDS")
rm(list=ls())

# Create monocle object from SC transformed seurat object with mito and colday regressed out 
seurat_file <- "./rds_objects/seurat_obj_lowpass_sctransform_mitocoldayregress.RDS"
sc_seurat <- readRDS(seurat_file)
expression_matrix <- sc_seurat@assays[["RNA"]]@scale.data
cell_metadata <- sc_seurat[[]]
gene_metadata <- as.data.frame(rownames(expression_matrix))
rownames(gene_metadata) <- rownames(expression_matrix)
colnames(gene_metadata) <- "gene_short_name"
monocle <- new_cell_data_set(expression_data = expression_matrix,
                             gene_metadata = gene_metadata,
                             cell_metadata = cell_metadata)
monocle <- monocle[,Matrix::colSums(exprs(monocle)) != 0]
saveRDS(monocle, "./rds_objects/monocle_from_seurat_sctransform_mitocoldayregress.RDS")
rm(list=ls())
