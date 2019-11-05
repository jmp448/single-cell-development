library(monocle3)

args = commandArgs(trailingOnly = TRUE)
source_files = args[1]

# m <- readRDS("./rds_objects/monocle_from_seurat_lowpass_unprocessed.RDS")
# m <- preprocess_cds(m, method="PCA", norm_method="log", scaling=T)
# m <- reduce_dimension(m, reduction_method="UMAP", preprocess_method="PCA", umap.min_dist=0.1, umap.n_neighbors=30L, umap.fast_sgd=F)
# saveRDS(m, "./rds_objects/monocle_from_seurat_lowpass_processed.RDS")
# rm(list=ls())
#
# m <- readRDS("./rds_objects/monocle_from_seurat_lowpass_basic_norm.RDS")
# m <- preprocess_cds(m, method="PCA", norm_method="none", scaling=F)
# m <- reduce_dimension(m, reduction_method="UMAP", preprocess_method="PCA", umap.min_dist=0.1, umap.n_neighbors=30L, umap.fast_sgd=F)
# saveRDS(m, "./rds_objects/monocle_from_seurat_lowpass_basic_norm_processed.RDS")
# rm(list=ls())
#
# m <- readRDS("./rds_objects/monocle_from_seurat_lowpass_noscale_norm.RDS")
# m <- preprocess_cds(m, method="PCA", norm_method="none", scaling=F)
# m <- reduce_dimension(m, reduction_method="UMAP", preprocess_method="PCA", umap.min_dist=0.1, umap.n_neighbors=30L, umap.fast_sgd=F)
# saveRDS(m, "./rds_objects/monocle_from_seurat_lowpass_noscale_norm_processed.RDS")
# rm(list=ls())
#
# m <- readRDS("./rds_objects/monocle_from_seurat_lowpass_sctransform.RDS")
# m <- preprocess_cds(m, method="PCA", norm_method="none", scaling=F)
# m <- reduce_dimension(m, reduction_method="UMAP", preprocess_method="PCA", umap.min_dist=0.1, umap.n_neighbors=30L, umap.fast_sgd=F)
# saveRDS(m, "./rds_objects/monocle_from_seurat_lowpass_sctransform_processed.RDS")
# rm(list=ls())

m <- readRDS(paste0("./rds_objects/monocle_from_seurat_", source_files, "_sctransform_mitoregress.RDS"))
m <- preprocess_cds(m, method="PCA", norm_method="none", scaling=F)
m <- reduce_dimension(m, reduction_method="UMAP", preprocess_method="PCA", umap.min_dist=0.3, umap.n_neighbors=30L, umap.fast_sgd=F)
saveRDS(m, paste0("./rds_objects/monocle_from_seurat_", source_files, "_sctransform_mitoregress_processed.RDS"))
rm(list=ls())
#
# m <- readRDS("./rds_objects/monocle_from_seurat_lowpass_sctransform_mitocoldayregress.RDS")
# m <- preprocess_cds(m, method="PCA", norm_method="none", scaling=F)
# m <- reduce_dimension(m, reduction_method="UMAP", preprocess_method="PCA", umap.min_dist=0.1, umap.n_neighbors=30L, umap.fast_sgd=F)
# saveRDS(m, "./rds_objects/monocle_from_seurat_lowpass_sctransform_mitocoldayregress_processed.RDS")
