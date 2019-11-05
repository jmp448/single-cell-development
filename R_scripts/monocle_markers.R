# MARKER GENE ANALYSIS

library(monocle3)

cds <- readRDS("./rds_objects/monocle_from_seurat_fulldata_sctransform_mitoregress_processed.RDS")

# Check out expression for specific genes
# Stem cell markers
stem_cell_markers <- c("EPCAM")
png("./monocle_plots/fulldata_markers_stem_cell.png", width=800, height=600)
plot_cells(cds, genes=stem_cell_markers,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
dev.off()

# Mesoderm
mesoderm_markers <- c("MIXL1")
png("./monocle_plots/fulldata_markers_mesoderm.png", width=800, height=600)
plot_cells(cds, genes=mesoderm_markers,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
dev.off()

# Cardiomyocyte
cardiomyocyte_markers <- c("TNNT2", "TNNI1", "TNNC1")
png("./monocle_plots/fulldata_markers_cardiomyocyte.png", width=1800, height=600)
plot_cells(cds, genes=cardiomyocyte_markers,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
dev.off()

# Epicardium
epicardium_markers = "IGF2"
png("./monocle_plots/fulldata_markers_epicardium.png", width=800, height=600)
plot_cells(cds, genes=epicardium_markers,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
dev.off()

# Fibroblast
fibroblast_markers <- c("VIM", "S100A4", "COL1A2", "COL3A1", "TSLP")
png("./monocle_plots/fulldata_markers_fibroblast.png", width=1800, height=600)
plot_cells(cds, genes=fibroblast_markers,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
dev.off()
