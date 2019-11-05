library(ggplot2)

cds <- readRDS("./rds_objects/monocle_from_seurat_fulldata_sctransform_mitoregress_graphed.RDS")

plt <- ggplot(as.data.frame(cds@colData), aes(x=pseudotime(cds), y=cds$diffday)) + 
  geom_point(aes(color=cds@assays[["counts"]]["TNNT2",]))

print(plt)