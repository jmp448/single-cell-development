library(ggplot2)
library(monocle3)

cds <- readRDS("./rds_objects/monocle_from_seurat_fulldata_sctransform_mitoregress_graphed.RDS")

plt <- ggplot(as.data.frame(cds@colData), aes(x=pseudotime(cds), y=cds$diffday, color=cds@assays[["counts"]]["TNNT2",], size=cds@assays[["counts"]]["TNNT2",])) + 
  geom_point() +
  scale_color_gradient(low="white", high="red")

print(plt)