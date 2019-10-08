# Playing with normalization

FeatureScatter(sc, feature1="nCount_SCT", feature2="nFeature_SCT")
FeatureScatter(sc, feature1="nCount_RNA", feature2="nFeature_RNA")

# Set reticulate library
library(reticulate)
reticulate::use_python('/software/python-3.5.2-el7-x86_64/bin/python', required=TRUE)

sc <- readRDS("scRNA_lowmito.rds")
sc$colday <- "colday"
sc$colday <- substr(sc$orig.ident, 3, 3)

sc <- SCTransform(sc, vars.to.regress=c("colday", "percent.mito"))
sc <- RunPCA(sc)
sc <- RunUMAP(sc, dims=1:30)

#comparison of normalization methods
sc_mito <- readRDS("scRNA_lowmito_processed.rds")
sc_colday <- readRDS("scRNA_lowmito_processed_colday.rds")

