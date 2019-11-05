library(monocle3)

cds <- readRDS("./rds_objects/monocle_from_seurat_fulldata_sctransform_mitoregress_processed.RDS")

cds <- cluster_cells(cds, resolution=1e-7)
plt <- plot_cells(cds)
print(plt)
plt <- plot_cells(cds, color_cells_by="partition", group_cells_by="partition")
print(plt)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plt <- plot_cells(cds,
           color_cells_by = "diffday",
           label_cell_groups=F,
           label_leaves=F,
           label_branch_points=F)
print(plt)
plt <- plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=F,
           label_leaves=F,
           label_branch_points=F)
print(plt)
