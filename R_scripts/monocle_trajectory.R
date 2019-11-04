library(monocle3)

m <- readRDS("./rds_objects/monocle_from_seurat_lowpass_sctransform_mitoregress_processed.RDS")
saveRDS(m, "./rds_objects/monocle_from_seurat_lowpass_sctransform_mitoregress_graphed.RDS")
rm(list=ls())

m <- readRDS("./rds_objects/monocle_from_seurat_lowpass_sctransform_mitocoldayregress_processed.RDS")
saveRDS(m, "./rds_objects/monocle_from_seurat_lowpass_sctransform_mitocoldayregress_graphed.RDS")
