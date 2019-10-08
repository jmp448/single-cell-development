# Upload seurat file
sc <- readRDS("sc_fulldata.RDS")

# See (in bulk) gene counts, mito counts
Idents(sc) <- "diffday"
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size=0)
Idents(sc) <- "individual"
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size=0)

# See (in bulk) cell counts per diff day, individual, and col day
plt <- barplot(table(sc$"diffday"), main=line_name)
print(plt)
dev.copy(png, paste0('diffday_barplot_bulk.png'))
dev.off()
plt <- barplot(table(sc$"individual"), main=line_name)
print(plt)
dev.copy(png, paste0('indiv_barplot_bulk.png'))
dev.off()
plt <- barplot(table(sc$"colday"), main=line_name)
print(plt)
dev.copy(png, paste0('colday_barplot_bulk.png'))
dev.off()

# Create a new seurat object for each cell line
Idents(sc) <- "individual"
na19093 <- subset(sc, idents="NA19093")
na18858 <- subset(sc, idents="NA18858")
na18520 <- subset(sc, idents="NA18520")
na18912 <- subset(sc, idents="NA18912")
na18508 <- subset(sc, idents="NA18508")
na18511 <- subset(sc, idents="NA18511")

# Create list of seurat objects, one for each cell line
lines <- c(na19093, na18858, na18912, na18520, na18508, na18511)

# Make cell count plots over time 
for (line in lines) {
  # Identify line name even when subset is not working
  for (r in rownames(table(line$"individual"))) {
    line_count <- table(line$"individual")[r]
    if (line_count > 0) {
      line_name = r
    } 
  }
  
  # Create plots for individual cell lines
  plt <- barplot(table(line$"diffday"), main=line_name)
  print(plt)
  dev.copy(png, paste0('diffday_barplot_', line_name, '.png'))
  dev.off()
  plt <- barplot(table(line$"colday"), main=line_name)
  print(plt)
  dev.copy(png, paste0('colday_barplot_', line_name, '.png'))
  dev.off()
  plt <- VlnPlot(line, features=c("nFeature_RNA","nCount_RNA","percent.mito"),ncol=3, group.by = "diffday", pt.size=0)
  print(plt)
  dev.copy(png, paste0('counts_vlnplot_', line_name, '.png'))
  dev.off()
}

# Set reticulate library
library(reticulate)
reticulate::use_python('/software/python-3.5.2-el7-x86_64/bin/python', required=TRUE)

# Run PCA and UMAP
sc <- SCTransform(sc, vars.to.regress=c("percent.mito", "colday"))
sc <- saveRDS('sc_fulldata.RDS')
sc <- RunPCA(sc)
sc <- RunUMAP(sc, dims=1:30)

# Play it safe, run normalization the long way
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
sc <- FindVariableFeatures(sc)
sc <- ScaleData(sc, vars.to.regress = c("percent.mito", "colday"))
sc <- RunPCA(sc, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
sc <- RunUMAP(sc, dims = 1:75, min.dist = 0.75)
saveRDS(sc, file="sc_crappy_transform.RDS")

# Define some new colors for easier visualization
reds <- hsv(h=1,s=seq(0.1, 1, length.out=7),v=1)
heat_cols <- heat.colors(7)
# Visualize UMAP
Idents(sc_full) <- "diffday"
DimPlot(sc_full, reduction="pca", cols=reds, group.by="diffday")
DimPlot(sc)
