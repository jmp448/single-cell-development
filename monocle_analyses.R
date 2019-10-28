library(monocle3)

# Create method for producing and saving plots
plots_loc <- "./monocle_plots/"
create_plot_png <- function(plt, filename, w = 800, h = 600, plt_title=NA) {
  if (!is.na(plt_title)) {
    plt <- plt + title(main=plt_title)
  }
  print(plt)
  dev.copy(png, paste0(plots_loc, filename), width = w, height = h)
  dev.off()
}

sc <- readRDS("./rds_objects/monocle_cds_fulldata.RDS")

# Try to regress out sample batch effects
# Attempt 1: regress out collection
sc2 <- preprocess_cds(sc, method="PCA", norm_method="log", alignment_group="collection")
sc2 <- reduce_dimension(sc, preprocess_method="PCA")
plt <- plot_cells(cds, label_cell_groups=F, color_cells_by="collection")
create_plot_png(plt, filename="umap_reg_sample_plot_sample.png")
plt <- plot_cells(cds, label_cell_groups=F, color_cells_by="diffday")
create_plot_png(plt, filename="umap_reg_sample_plot_diffday.png")
plt <- plot_cells(cds, label_cell_groups=F, color_cells_by="individual")
create_plot_png(plt, filename="umap_reg_sample_plot_indiv.png")
saveRDS(sc2, "monocle_cds_fulldata_processed_collection.RDS")
rm(sc2)

# Attempt 2: regress out collection day
sc2 <- preprocess_cds(sc, method="PCA", norm_method="log", alignment_group="colday")
sc2 <- reduce_dimension(sc, preprocess_method="PCA")
plt <- plot_cells(cds, label_cell_groups=F, color_cells_by="collection")
create_plot_png(plt, filename="umap_reg_colday_plot_sample.png")
plt <- plot_cells(cds, label_cell_groups=F, color_cells_by="diffday")
create_plot_png(plt, filename="umap_reg_colday_plot_diffday.png")
plt <- plot_cells(cds, label_cell_groups=F, color_cells_by="individual")
create_plot_png(plt, filename="umap_reg_colday_plot_indiv.png")
saveRDS(sc2, "monocle_cds_fulldata_processed_colday.RDS")
rm(sc2)
