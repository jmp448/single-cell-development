# MARKER GENE ANALYSIS

# Get Marker Genes
# After clustering, Idents(sc) is set to the cluster number
cluster4.markers <- FindMarkers(sc, ident.1=4, min.pct=0.25)
head(cluster4.markers, n=5)
cluster0.markers <- FindMarkers(sc, ident.1=0, min.pct=0.25)
head(cluster0.markers, n=5)
cluster1.markers <- FindMarkers(sc, ident.1=1, min.pct=0.25)
head(cluster1.markers, n=5)
# Markers for diffday
Idents(sc) <- "diffday"
day15.markers <- FindMarkers(sc, ident.1="Day 15", min.pct=0.25)
day0.markers <- FindMarkers(sc, ident.1="Day 0", min.pct=0.25)
day11.markers <- FindMarkers(sc, ident.1="Day 11", min.pct=0.25)
day1.markers <- FindMarkers(sc, ident.1="Day 1", min.pct=0.25)
day0_vs_day15.markers <- FindMarkers(sc, ident.1="Day 0", ident.2 = "Day 15", min.pct = 0.25)
day0_vs_day11.markers <- FindMarkers(sc, ident.1="Day 0", ident.2 = "Day 11", min.pct = 0.25)
day1_vs_day15.markers <- FindMarkers(sc, ident.1="Day 1", ident.2 = "Day 15", min.pct = 0.25)  # 
day1_vs_day11.markers <- FindMarkers(sc, ident.1="Day 1", ident.2 = "Day 11", min.pct = 0.25)

# Try to create subsets for later analysis
gene_list <- rownames(day11.markers)
day11marker.subset <- subset(sc, features = gene_list)
day11marker.subset[["RNA"]]@data

# Check out expression for specific genes
# Stem cell markers
stem_cell_markers <- c("EPCAM", "LIN28A", "ZFP42")
VlnPlot(sc, features=stem_cell_markers, group.by="diffday", split.by="individual", pt.size=0)
VlnPlot(sc, features="EPCAM", group.by="diffday", split.by="individual", pt.size=0) # aka TROP1
VlnPlot(sc, features="LIN28A", group.by="diffday", split.by="individual", pt.size=0)
VlnPlot(sc, features="ZFP42", group.by="diffday", split.by="individual", pt.size=0)
VlnPlot(sc, features="FOXD3", group.by="diffday", split.by="individual")
VlnPlot(sc, features="POU5F1", group.by="diffday", split.by="individual", pt.size=0)  #OCT3 OCT4
VlnPlot(sc, features="NANOG", group.by="diffday", split.by="individual", pt.size=0)
VlnPlot(sc, features="PROM1", group.by="diffday", split.by="individual", pt.size=0) # aka CD133
VlnPlot(sc, features="KLF5", group.by="diffday", split.by="individual", pt.size=0) 
VlnPlot(sc, features="KLF2", group.by="diffday", split.by="individual", pt.size=0)
VlnPlot(sc, features="STAT3", group.by="diffday", split.by="individual", pt.size=0)
VlnPlot(sc, features="UTF1", group.by="diffday", split.by="individual", pt.size=0)
# Not present: ALPP
# Couldn't find genes: TRA-1-60

# Mesoderm 
mesoderm_markers <- c("CXCR4", "KDR", "EOMES", "FOXF1", "MIXL1", "TBX6")
VlnPlot(sc, features=mesoderm_markers, group.by="diffday", split.by="individual")
VlnPlot(sc, features="MIXL1", group.by="diffday", split.by="individual", pt.size=0)

# Cardiac progenitor to cardiomyocyte transition
progenitor_markers <- c("KDR", "ATF2", "GATA2", "MEF2C", "TBX5", "WT1")
VlnPlot(sc, features=progenitor_markers, group.by="diffday", split.by="individual", pt.size=0)

# Cardiomyocyte 
cardiomyocyte_markers <- c("TNNT2", "TNNI1", "GJA1", "VCAM1", "TNNC1", "MYH7")
VlnPlot(sc, features=cardiomyocyte_markers, group.by="diffday", split.by="individual", pt.size=0)
VlnPlot(sc, features="TNNT2", group.by="diffday", split.by="individual", pt.size=0)

# Cardiac Fibroblast
VlnPlot(sc, features="DDR2", group.by="diffday", split.by="individual")

# Myofibroblast
VlnPlot(sc, features="ACTA2", group.by="diffday", split.by="individual")

# Epicardium
VlnPlot(sc, features="IGF2", group.by="diffday", split.by="individual", pt.size=0)

# Fibroblast
fibroblast_markers <- c("VIM", "S100A4", "COL1A2", "COL3A1", "TSLP")
VlnPlot(sc, features=fibroblast_markers, group.by="diffday", split.by="individual", pt.size=0)
VlnPlot(sc, features="VIM", group.by="diffday", split.by="individual", pt.size=0)

FeaturePlot(sc, features=c("EPCAM", "MIXL1", "TNNT2", "TNNI1"))