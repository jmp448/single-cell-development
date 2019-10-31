library(Seurat)

args = commandArgs(trailingOnly = TRUE)
source_files = args[1]

sc <- readRDS("./rds_objects/old_seurat_obj_lowpass.rds")
sc <- SCTransform(sc, vars.to.regress=c("percent.mito"))
saveRDS(sc, "./rds_objects/old_seurat_obj_lowpass_processed.rds")
# sc <- readRDS(paste0("./rds_objects/seurat_obj_", source_files, ".RDS"))
#
# # Seurat basic normalization
# sc1 <- NormalizeData(sc)
# sc1 <- ScaleData(sc1)
# saveRDS(sc1, paste0("./rds_objects/seurat_obj_", source_files, "_basic_norm.RDS"))
# rm(sc1)
#
# # Proposed monocle normalization?
# sc1 <- NormalizeData(sc, scale.factor=1)
# sc1 <- ScaleData(sc1)
# saveRDS(sc1, paste0("./rds_objects/seurat_obj_", source_files, "_noscale_norm.RDS"))
# rm(sc1)
#
# # SCTransform
# sc1 <- SCTransform(sc)
# saveRDS(sc1, paste0("./rds_objects/seurat_obj_", source_files, "_sctransform.RDS"))
# rm(sc1)
#
# # SCTransform with mito regress
# sc1 <- SCTransform(sc, vars.to.regress = c("percent.mito"))
# saveRDS(sc1, paste0("./rds_objects/seurat_obj_", source_files, "_sctransform_mitoregress.RDS"))
# rm(sc1)
#
# # SCTransform with mito regress
# sc1 <- SCTransform(sc, vars.to.regress = c("percent.mito", "colday"))
# saveRDS(sc1, paste0("./rds_objects/seurat_obj_", source_files, "_sctransform_mitocoldayregress.RDS"))
