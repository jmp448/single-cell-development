library(Seurat)
# Set reticulate library
library(reticulate)
reticulate::use_python("/software/python-3.5.2-el7-x86_64/bin/python", required = TRUE)

sc <- readRDS("sc_fulldata.RDS")

# Run PCA and UMAP
sc <- SCTransform(sc, vars.to.regress = c("percent.mito", "colday"))
sc <- RunPCA(sc)
sc <- RunUMAP(sc, dims = 1:30)

sc <- saveRDS("sc_fulldata_sctransformed.RDS")
