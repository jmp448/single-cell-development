library(Seurat)

args = commandArgs(trailingOnly = TRUE)
source_files = args[1]

# sc <- readRDS("./rds_objects/seurat_obj_lowpass_basic_norm.RDS")
# sc <- FindVariableFeatures(sc)
# sc <- RunPCA(sc)
# sc <- RunUMAP(sc, dims=1:30)
# saveRDS(sc, "./rds_objects/seurat_obj_lowpass_basic_norm_processed.RDS")
# rm(list=ls())
#
# sc <- readRDS("./rds_objects/seurat_obj_lowpass_noscale_norm.RDS")
# sc <- FindVariableFeatures(sc)
# sc <- RunPCA(sc)
# sc <- RunUMAP(sc, dims=1:30)
# saveRDS(sc, "./rds_objects/seurat_obj_lowpass_noscale_norm_processed.RDS")
# rm(list=ls())
#
# sc <- readRDS("./rds_objects/seurat_obj_lowpass_sctransform.RDS")
# sc <- RunPCA(sc)
# sc <- RunUMAP(sc, dims=1:30)
# saveRDS(sc, "./rds_objects/seurat_obj_lowpass_sctransform_processed.RDS")
# rm(list=ls())

sc <- readRDS(paste0("./rds_objects/seurat_obj_", source_files, "_sctransform_mitoregress.RDS"))
sc <- RunPCA(sc)
sc <- RunUMAP(sc, dims=1:30)
saveRDS(sc, paste0("./rds_objects/seurat_obj_", source_files, "_sctransform_mitoregress_processed.RDS"))
rm(list=ls())

# sc <- readRDS("./rds_objects/seurat_obj_lowpass_sctransform_mitocoldayregress.RDS")
# sc <- RunPCA(sc)
# sc <- RunUMAP(sc, dims=1:30)
# saveRDS(sc, "./rds_objects/seurat_obj_lowpass_sctransform_mitocoldayregress_processed.RDS")
