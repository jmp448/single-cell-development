library(monocle3)

plots_loc <- "./monocle_plots/"
create_plot_png <- function(plt, filename, w = 800, h = 600) {
  print(plt)
  dev.copy(png, paste0(plots_loc, filename), width = w, height = h)
  dev.off()
}

sc <- readRDS("./rds_objects/sc_fulldata_monocle.RDS")

sc <- preprocess_cds(sc, method="PCA", norm_method="log", scaling=T, verbose=T)
sc <- reduce_dimension(sc, verbose=T, preprocess_method="PCA")

plt <- plot_cells(sc, label_cell_groups=F, reduction_method="PCA", color_cells_by="diffday")
create_plot_png(plt, "pca_diffday.png")

plt <- plot_cells(sc, label_cell_groups=F, reduction_method="PCA", color_cells_by="individual")
create_plot_png(plt, "pca_individual.png")

plt <- plot_cells(sc, label_cell_groups=F, reduction_method="UMAP", color_cells_by="diffday")
create_plot_png(plt, "umap_diffday_2.png")

plt <- plot_cells(sc, label_cell_groups=F, reduction_method="UMAP", color_cells_by="individual")
create_plot_png(plt, "umap_individual_2.png")

scp@colData["sample"] <- NA
for (i in length(rownames(scp@colData))) {
  scp@colData["sample"][[1]][[i]] <- substr(rownames(scp@colData)[[i]], 14,20)
}
