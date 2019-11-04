# MARKER GENE ANALYSIS

cds <- readRDS("./rds_objects/monocle_from_seurat_lowpass_sctransform_mitoregress_processed.RDS")

# Check out expression for specific genes
# Stem cell markers
stem_cell_markers <- c("EPCAM", "LIN28A", "ZFP42")
png("./monocle_plots/stem_cell_markers", width=800, height=600)
plot_cells(cds, genes=stem_cell_markers,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
dev.off()

# Mesoderm
mesoderm_markers <- c("CXCR4", "KDR", "EOMES", "FOXF1", "MIXL1", "TBX6")
png("./monocle_plots/mesoderm_markers", width=800, height=600)
plot_cells(cds, genes=mesoderm_markers,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
dev.off()

# Cardiac progenitor to cardiomyocyte transition
progenitor_markers <- c("KDR", "ATF2", "GATA2", "MEF2C", "TBX5", "WT1")
png("./monocle_plots/progenitor_markers", width=800, height=600)
plot_cells(cds, genes=progenitor_markers,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
dev.off()

# Cardiomyocyte
cardiomyocyte_markers <- c("TNNT2", "TNNI1", "GJA1", "VCAM1", "TNNC1", "MYH7")
png("./monocle_plots/cardiomyocyte_markers", width=800, height=600)
plot_cells(cds, genes=cardiomyocyte_markers,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
dev.off()

# Myofibroblast
myofibroblast_markers <- "ACTA2"
png("./monocle_plots/myofibroblast_markers", width=800, height=600)
plot_cells(cds, genes=myofibroblast_markers,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
dev.off()

# Epicardium
epicardium_markers = "IGF2"
png("./monocle_plots/epicardium_markers", width=800, height=600)
plot_cells(cds, genes=epicardium_markers,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
dev.off()

# Fibroblast
fibroblast_markers <- c("VIM", "S100A4", "COL1A2", "COL3A1", "TSLP")
png("./monocle_plots/fibroblast_markers", width=800, height=600)
plot_cells(cds, genes=fibroblast_markers,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
dev.off()
