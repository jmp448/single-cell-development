library(Seurat)

# plots_loc <- "./seurat_plots/"
#
# ### Basic norm
# sc <- readRDS("./rds_objects/seurat_obj_lowpass_basic_norm_processed.RDS")
#
# png(paste0(plots_loc, "lowpass_basic_norm_PCA_diffday.png"), width=800, height=600)
# PCAPlot(sc, group.by="diffday")
# dev.off()
# png(paste0(plots_loc, "lowpass_basic_norm_PCA_indiv.png"), width=800, height=600)
# PCAPlot(sc, group.by="individual")
# dev.off()
# png(paste0(plots_loc, "lowpass_basic_norm_PCA_col.png"), width=800, height=600)
# PCAPlot(sc, group.by="orig.ident")
# dev.off()
# png(paste0(plots_loc, "lowpass_basic_norm_PCA_colday.png"), width=800, height=600)
# PCAPlot(sc, group.by="colday")
# dev.off()
#
# png(paste0(plots_loc, "lowpass_basic_norm_UMAP_diffday.png"), width=800, height=600)
# DimPlot(sc, group.by="diffday")
# dev.off()
# png(paste0(plots_loc, "lowpass_basic_norm_UMAP_indiv.png"), width=800, height=600)
# DimPlot(sc, group.by="individual")
# dev.off()
# png(paste0(plots_loc, "lowpass_basic_norm_UMAP_col.png"), width=800, height=600)
# DimPlot(sc, group.by="orig.ident")
# dev.off()
# png(paste0(plots_loc, "lowpass_basic_norm_UMAP_colday.png"), width=800, height=600)
# DimPlot(sc, group.by="colday")
# dev.off()
# rm(list=ls())
#
# plots_loc <- "./seurat_plots/"
# ### Noscale norm
# sc <- readRDS("./rds_objects/seurat_obj_lowpass_noscale_norm_processed.RDS")
#
# # PCA diffday
# png(paste0(plots_loc, "lowpass_noscale_norm_PCA_diffday.png"), width=800, height=600)
# PCAPlot(sc, group.by="diffday")
# dev.off()
#
# # PCA individual
# png(paste0(plots_loc, "lowpass_noscale_norm_PCA_indiv.png"), width=800, height=600)
# PCAPlot(sc, group.by="individual")
# dev.off()
#
# # PCA collection
# png(paste0(plots_loc, "lowpass_noscale_norm_PCA_col.png"), width=800, height=600)
# PCAPlot(sc, group.by="orig.ident")
# dev.off()
#
# # PCA colday
# png(paste0(plots_loc, "lowpass_noscale_norm_PCA_colday.png"), width=800, height=600)
# PCAPlot(sc, group.by="colday")
# dev.off()
#
# # UMAP diffday
# png(paste0(plots_loc, "lowpass_noscale_norm_UMAP_diffday.png"), width=800, height=600)
# DimPlot(sc, group.by="diffday")
# dev.off()
#
# # UMAP individual
# png(paste0(plots_loc, "lowpass_noscale_norm_UMAP_indiv.png"), width=800, height=600)
# DimPlot(sc, group.by="individual")
# dev.off()
#
# # UMAP collection
# png(paste0(plots_loc, "lowpass_noscale_norm_UMAP_col.png"), width=800, height=600)
# DimPlot(sc, group.by="orig.ident")
# dev.off()
#
# # UMAP colday
# png(paste0(plots_loc, "lowpass_noscale_norm_UMAP_colday.png"), width=800, height=600)
# DimPlot(sc, group.by="colday")
# dev.off()
#
# rm(list=ls())
#
# plots_loc <- "./seurat_plots/"
#
#
# ### SCTransform
# sc <- readRDS("./rds_objects/seurat_obj_lowpass_sctransform_processed.RDS")
#
# # PCA diffday
# png(paste0(plots_loc, "lowpass_sctransform_PCA_diffday.png"), width=800, height=600)
# PCAPlot(sc, group.by="diffday")
# dev.off()
#
# # PCA individual
# png(paste0(plots_loc, "lowpass_sctransform_PCA_indiv.png"), width=800, height=600)
# PCAPlot(sc, group.by="individual")
# dev.off()
#
# # PCA collection
# png(paste0(plots_loc, "lowpass_sctransform_PCA_col.png"), width=800, height=600)
# PCAPlot(sc, group.by="orig.ident")
# dev.off()
#
# # PCA colday
# png(paste0(plots_loc, "lowpass_sctransform_PCA_colday.png"), width=800, height=600)
# PCAPlot(sc, group.by="colday")
# dev.off()
#
# # UMAP diffday
# png(paste0(plots_loc, "lowpass_sctransform_UMAP_diffday.png"), width=800, height=600)
# DimPlot(sc, group.by="diffday")
# dev.off()
#
# # UMAP individual
# png(paste0(plots_loc, "lowpass_sctransform_UMAP_indiv.png"), width=800, height=600)
# DimPlot(sc, group.by="individual")
# dev.off()
#
# # UMAP collection
# png(paste0(plots_loc, "lowpass_sctransform_UMAP_col.png"), width=800, height=600)
# DimPlot(sc, group.by="orig.ident")
# dev.off()
#
# # UMAP colday
# png(paste0(plots_loc, "lowpass_sctransform_UMAP_colday.png"), width=800, height=600)
# DimPlot(sc, group.by="colday")
# dev.off()
#
# rm(list=ls())
#
plots_loc <- "./seurat_plots/"

### SCTransform Mitoregress
sc <- readRDS("./rds_objects/seurat_obj_fulldata_sctransform_mitoregress_processed.RDS")

# PCA diffday
png(paste0(plots_loc,"fulldata_sctransform_mitoregress_PCA_diffday.png"), width=800, height=600)
PCAPlot(sc, group.by="diffday")
dev.off()

# PCA individual
png(paste0(plots_loc, "fulldata_sctransform_mitoregress_PCA_indiv.png"), width=800, height=600)
PCAPlot(sc, group.by="individual")
dev.off()

# PCA collection
png(paste0(plots_loc, "fulldata_sctransform_mitoregress_PCA_col.png"), width=800, height=600)
PCAPlot(sc, group.by="orig.ident")
dev.off()

# PCA colday
png(paste0(plots_loc, "fulldata_sctransform_mitoregress_PCA_colday.png"), width=800, height=600)
PCAPlot(sc, group.by="colday")
dev.off()

# UMAP diffday
png(paste0(plots_loc, "fulldata_sctransform_mitoregress_UMAP_diffday.png"), width=800, height=600)
DimPlot(sc, group.by="diffday")
dev.off()

# UMAP individual
png(paste0(plots_loc, "fulldata_sctransform_mitoregress_UMAP_indiv.png"), width=800, height=600)
DimPlot(sc, group.by="individual")
dev.off()

# UMAP collection
png(paste0(plots_loc, "fulldata_sctransform_mitoregress_UMAP_col.png"), width=800, height=600)
DimPlot(sc, group.by="orig.ident")
dev.off()

# UMAP colday
png(paste0(plots_loc, "fulldata_sctransform_mitoregress_UMAP_colday.png"), width=800, height=600)
DimPlot(sc, group.by="colday")
dev.off()

rm(list=ls())

# plots_loc <- "./seurat_plots/"
#
#
# ### SCTransform Mito Colday Regress
# sc <- readRDS("./rds_objects/seurat_obj_lowpass_sctransform_mitocoldayregress_processed.RDS")
#
# # PCA diffday
# png(paste0(plots_loc, "lowpass_sctransform_mitocoldayregress_PCA_diffday.png"), width=800, height=600)
# PCAPlot(sc, group.by="diffday")
# dev.off()
#
# # PCA individual
# png(paste0(plots_loc, "lowpass_sctransform_mitocoldayregress_PCA_indiv.png"), width=800, height=600)
# PCAPlot(sc, group.by="individual")
# dev.off()
#
# # PCA collection
# png(paste0(plots_loc, "lowpass_sctransform_mitocoldayregress_PCA_col.png"), width=800, height=600)
# PCAPlot(sc, group.by="orig.ident")
# dev.off()
#
# # PCA colday
# png(paste0(plots_loc, "lowpass_sctransform_mitocoldayregress_PCA_colday.png"), width=800, height=600)
# PCAPlot(sc, group.by="colday")
# dev.off()
#
# # UMAP diffday
# png(paste0(plots_loc, "lowpass_sctransform_mitocoldayregress_UMAP_diffday.png"), width=800, height=600)
# DimPlot(sc, group.by="diffday")
# dev.off()
#
# # UMAP individual
# png(paste0(plots_loc, "lowpass_sctransform_mitocoldayregress_UMAP_indiv.png"), width=800, height=600)
# DimPlot(sc, group.by="individual")
# dev.off()
#
# # UMAP collection
# png(paste0(plots_loc, "lowpass_sctransform_mitocoldayregress_UMAP_col.png"), width=800, height=600)
# DimPlot(sc, group.by="orig.ident")
# dev.off()
#
# # UMAP colday
# png(paste0(plots_loc, "lowpass_sctransform_mitocoldayregress_UMAP_colday.png"), width=800, height=600)
# DimPlot(sc, group.by="colday")
# dev.off()
#
# rm(list=ls())
