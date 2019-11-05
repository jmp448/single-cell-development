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
plt <- plot_genes_in_pseudotime(cds["EPCAM",], min_expr = NULL, cell_size = 0.75,nrow = NULL, ncol = 1, panel_order = NULL,color_cells_by = "pseudotime",trend_formula = "~ splines::ns(pseudotime, df=3)",label_by_short_name = TRUE, vertical_jitter = NULL,horizontal_jitter = NULL)
print(plt)
plt <- plot_genes_in_pseudotime(cds["TNNT2",], min_expr = NULL, cell_size = 0.75,nrow = NULL, ncol = 1, panel_order = NULL,color_cells_by = "pseudotime",trend_formula = "~ splines::ns(pseudotime, df=3)",label_by_short_name = TRUE, vertical_jitter = NULL,horizontal_jitter = NULL)
print(plt)
plt <- plot_genes_in_pseudotime(cds["IGF2",], min_expr = NULL, cell_size = 0.75,nrow = NULL, ncol = 1, panel_order = NULL,color_cells_by = "pseudotime",trend_formula = "~ splines::ns(pseudotime, df=3)",label_by_short_name = TRUE, vertical_jitter = NULL,horizontal_jitter = NULL)
print(plt)