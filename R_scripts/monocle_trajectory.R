library(monocle3)

cds <- readRDS("./rds_objects/monocle_from_seurat_lowpass_sctransform_mitoregress_processed.RDS")

library(monocle3)

cds=cluster_cells(cds, resolution=1e-7)
plt <- plot_cells(cds)
print(plt)
plt <- plot_cells(cds, color_cells_by="partition", group_cells_by="partition")
print(plt)
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "diffday",
           label_groups_by_cluster=F,
           label_leaves=F,
           label_branch_points=F)
cds <- order_cells(cds)
