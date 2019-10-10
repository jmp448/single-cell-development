# Create a series of barplots and violin plots to visualize the data before
# processing
library(Seurat)
library(stringr)

# Specify folder to put all plots
plots_loc <- "./plots/"

create_plot_png <- function(plt, filename, w = 480, h = 480) {
  print(plt)
  dev.copy(png, paste0(plots_loc, filename), width = w, height = h)
  dev.off()
}

# Upload seurat file
sc <- readRDS("./rds_objects/sc_fulldata.RDS")

# Create a new seurat object for each cell line
Idents(sc) <- "individual"
na19093 <- subset(sc, idents = "NA19093")
na18858 <- subset(sc, idents = "NA18858")
na18520 <- subset(sc, idents = "NA18520")
na18912 <- subset(sc, idents = "NA18912")
na18508 <- subset(sc, idents = "NA18508")
na18511 <- subset(sc, idents = "NA18511")

# Create list of seurat objects, one for each cell line
cell_lines <- c(na19093, na18858, na18912, na18520, na18508, na18511)

# See (in bulk) gene counts, mito counts
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, 
  pt.size = 0)

# See (in bulk) cell counts per diff day, individual, and col day
plt <- barplot(table(sc$diffday), main = line_name)
create_plot_png(plt, "diffday_barplot_bulk.png")

plt <- barplot(table(sc$individual), main = line_name)
create_plot_png(plt, "indiv_barplot_bulk.png")

plt <- barplot(table(sc$colday), main = line_name)
create_plot_png(plt, "colday_barplot_bulk.png")

# # Make cell count plots over time for (line in cell_lines) { # Identify line
# name even when subset is not working for (r in
# rownames(table(line$individual))) { line_count <- table(line$individual)[r] if
# (line_count > 0) { line_name = r } } # Create plots for individual cell lines
# plt <- barplot(table(line$diffday), main = line_name) create_plot_png(plt,
# paste0('diffday_barplot_', line_name, '.png'), w = 1200, h = 500) plt <-
# barplot(table(line$colday), main = line_name) create_plot_png(plt,
# paste0('colday_barplot_', line_name, '.png'), w = 1200, h = 500) Idents(line)
# <- 'diffday' plt <- VlnPlot(line, features = 'percent.mito', pt.size = 0)
# create_plot_png(plt, paste0('percentmito_vlnplot_', line_name, '.png'), w =
# 1200, h = 500) }
