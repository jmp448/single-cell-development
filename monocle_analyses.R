library(monocle3)

# Create method for producing and saving plots
plots_loc <- "./monocle_plots/"
create_plot_png <- function(plt, filename, w=800, h=600, plt_title=NA) {
  if (!is.na(plt_title)) {
    plt <- plt + title(main=plt_title)
  }
  print(plt)
  dev.copy(png, paste0(plots_loc, filename), width = w, height = h)
  dev.off()
}

sc <- readRDS("./rds_objects/monocle_cds_fulldata.RDS")
# colday_levels <- c("1", "2", "3")
# sc@colData$colday <- factor(x=sc@colData$colday, levels=colday_levels, ordered=T)
# diffday_levels <- c("Day 0", "Day 1", "Day 3", "Day 5", "Day 7", "Day 11", "Day 15")
# sc@colData$diffday <- factor(x=sc@colData$diffday, levels=diffday_levels, ordered=T)
# indiv_levels <- c("NA18508", "NA18511", "NA18520", "NA18858", "NA18912", "NA19093")
# sc@colData$individual <- factor(x=sc@colData$individual, levels=indiv_levels, ordered=T)
# col_levels <- c("CD1col1", "CD1col2", "CD1col3", "CD1col4", "CD1col5", "CD1col6", "CD2col1", "CD2col2", "CD2col3", "CD2col4", "CD2col5", "CD2col6", "CD3col1", "CD3col2", "CD3col3", "CD3col4", "CD3col5", "CD3col6")
# sc@colData$collection <- factor(x=sc@colData$collection, levels=col_levels, ordered=T)
# saveRDS(sc, "./rds_objects/monocle_cds_fulldata.RDS")

# Run UMAP with same params as seurat
sc2 <- preprocess_cds(sc, method="PCA", norm_method="log")
sc2 <- reduce_dimension(sc2, preprocess_method="PCA", umap.min_dist=0.3, umap.n_neighbors=30L)
plt <- plot_cells(sc2, label_cell_groups=F, color_cells_by="collection")
create_plot_png(plt, filename="umap_plot_sample.png")
plt <- plot_cells(sc2, label_cell_groups=F, color_cells_by="diffday")
create_plot_png(plt, filename="umap_plot_diffday.png")
plt <- plot_cells(sc2, label_cell_groups=F, color_cells_by="individual")
create_plot_png(plt, filename="umap_plot_indiv.png")
saveRDS(sc2, "./rds_objects/monocle_cds_fulldata_processed_seuratlike.RDS")
rm(sc2)

# Try to regress out sample batch effects
# Attempt 1: regress out collection
sc2 <- preprocess_cds(sc, method="PCA", norm_method="log")
sc2 <- align_cds(sc2, preprocess_method="PCA", alignment_group="collection")
sc2 <- reduce_dimension(sc2, preprocess_method="PCA", umap.min_dist=0.3, umap.n_neighbors=30L)
plt <- plot_cells(sc2, label_cell_groups=F, color_cells_by="collection")
create_plot_png(plt, filename="umap_reg_sample_plot_sample.png")
plt <- plot_cells(sc2, label_cell_groups=F, color_cells_by="diffday")
create_plot_png(plt, filename="umap_reg_sample_plot_diffday.png")
plt <- plot_cells(sc2, label_cell_groups=F, color_cells_by="individual")
create_plot_png(plt, filename="umap_reg_sample_plot_indiv.png")
saveRDS(sc2, "./rds_objects/monocle_cds_fulldata_processed_collection.RDS")
rm(sc2)

# Attempt 2: regress out collection day
sc2 <- preprocess_cds(sc, method="PCA", norm_method="log")
sc2 <- align_cds(sc2, preprocess_method="PCA", alignment_group="colday")
sc2 <- reduce_dimension(sc2, preprocess_method="PCA", umap.min_dist=0.3, umap.n_neighbors=30L)
plt <- plot_cells(sc2, label_cell_groups=F, color_cells_by="collection")
create_plot_png(plt, filename="umap_reg_colday_plot_sample.png")
plt <- plot_cells(sc2, label_cell_groups=F, color_cells_by="diffday")
create_plot_png(plt, filename="umap_reg_colday_plot_diffday.png")
plt <- plot_cells(sc2, label_cell_groups=F, color_cells_by="individual")
create_plot_png(plt, filename="umap_reg_colday_plot_indiv.png")
saveRDS(sc2, "./rds_objects/monocle_cds_fulldata_processed_colday.RDS")
rm(sc2)
