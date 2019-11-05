library(biomaRt)
library(stringr)
library(Seurat)

### The following is, until the end, copied directly from Reem Elorbany's work
### found at
### https://relorbany.github.io/Exp1_prelim_analysis/scCM_round1_lowpass_combined.html

## My goal is to combine all 18 collections into a single Seurat object.  First, I
## need to read in the raw data for all collections.
load_from_bash <- T
if (load_from_bash) {
  args = commandArgs(trailingOnly = TRUE)
  rawdata <- c()
  for (i in 1:18) {
    rawdata <- c(rawdata, args[i])
  }
  min_cells_per_gene <- args[19]  # minimum num cells in which a gene must appear
  min_genes_per_cell <- args[20]  # minimum num genes for a cell to be included
  mito_threshold <- args[21]  # what is the threshold for mito cutoff?
  source_files <- args[22]
} else {
  print("Going with R file access")
  rawdata <- c("/project2/gilad/reem/singlecellCM/round1/lowpass/CD1/CD1col1/output/dge_data/YG-RE1-Drop-CD1col1_S1_gene_counts.tsv.gz",
  "/project2/gilad/reem/singlecellCM/round1/lowpass/CD1/CD1col2/output/dge_data/YG-RE1-Drop-CD1col2_S2_gene_counts.tsv.gz",
  "/project2/gilad/reem/singlecellCM/round1/lowpass/CD1/CD1col3/output/dge_data/YG-RE1-Drop-CD1col3_S3_gene_counts.tsv.gz",
  "/project2/gilad/reem/singlecellCM/round1/lowpass/CD1/CD1col4/output/dge_data/YG-RE2-Drop-CD1col4_S1_gene_counts.tsv.gz",
  "/project2/gilad/reem/singlecellCM/round1/lowpass/CD1/CD1col5/output/dge_data/YG-RE2-Drop-CD1col5_S2_gene_counts.tsv.gz",
  "/project2/gilad/reem/singlecellCM/round1/lowpass/CD1/CD1col6/output/dge_data/YG-RE2-Drop-CD1col6_S3_gene_counts.tsv.gz",
  "/project2/gilad/reem/singlecellCM/round1/lowpass/CD2/CD2col1/output/dge_data/YG-RE1-Drop-CD2col1_S4_gene_counts.tsv.gz",
  "/project2/gilad/reem/singlecellCM/round1/lowpass/CD2/CD2col2/output/dge_data/YG-RE2-Drop-CD2col2_S4_gene_counts.tsv.gz",
  "/project2/gilad/reem/singlecellCM/round1/lowpass/CD2/CD2col3/output/dge_data/YG-RE2-Drop-CD2col3_S5_gene_counts.tsv.gz",
  "/project2/gilad/reem/singlecellCM/round1/lowpass/CD2/CD2col4/output/dge_data/YG-RE1-Drop-CD2col4_S5_gene_counts.tsv.gz",
  "/project2/gilad/reem/singlecellCM/round1/lowpass/CD2/CD2col5/output/dge_data/YG-RE2-Drop-CD2col5_S6_gene_counts.tsv.gz",
  "/project2/gilad/reem/singlecellCM/round1/lowpass/CD2/CD2col6/output/dge_data/YG-RE1-Drop-CD2col6_S6_gene_counts.tsv.gz",
  "/project2/gilad/reem/singlecellCM/round1/lowpass/CD3/CD3col1/output/dge_data/YG-RE2-Drop-CD3col1_S7_gene_counts.tsv.gz",
  "/project2/gilad/reem/singlecellCM/round1/lowpass/CD3/CD3col2/output/dge_data/YG-RE2-Drop-CD3col2_S8_gene_counts.tsv.gz",
  "/project2/gilad/reem/singlecellCM/round1/lowpass/CD3/CD3col3/output/dge_data/YG-RE1-Drop-CD3col3_S7_gene_counts.tsv.gz",
  "/project2/gilad/reem/singlecellCM/round1/lowpass/CD3/CD3col4/output/dge_data/YG-RE1-Drop-CD3col4_S8_gene_counts.tsv.gz",
  "/project2/gilad/reem/singlecellCM/round1/lowpass/CD3/CD3col5/output/dge_data/YG-RE1-Drop-CD3col5_S9_gene_counts.tsv.gz",
  "/project2/gilad/reem/singlecellCM/round1/lowpass/CD3/CD3col6/output/dge_data/YG-RE2-Drop-CD3col6_S9_gene_counts.tsv.gz")
  min_cells_per_gene <- 3
  min_genes_per_cell <- 200
  mito_threshold <- 30
}

## Next, I convert Ensembl gene numbers to gene names/symbols to make it more
## readable.  To do that, first I need to get the matched Ensembl IDs and gene
## symbols through biomart.  get data from biomart
geneinfo <- readRDS("./rds_objects/geneinfo.rds")

for (i in 1:3) {
  for (j in 1:6) {
    # read in raw data matrix
    expression_matrix <- read.table(rawdata[6*(i-1)+j], header = T, stringsAsFactors = F, row.names = 1)

    # remove version numbers from gene IDs
    genes_observed <- str_replace(row.names(expression_matrix), pattern = ".[0-9]+$",
        replacement = "")
    # get rid of duplicate ensembl IDs in expression matrix
    expression_matrix <- expression_matrix[!duplicated(genes_observed),]
    row.names(expression_matrix) <- genes_observed
    rm(genes_observed)

    # find genes in geneinfo that are also in cells, and sort by ensemblID
    genes_present <- geneinfo[geneinfo$ensembl_gene_id %in% row.names(expression_matrix), ]
    genes_present <- genes_present[!duplicated(genes_present$ensembl_gene_id),]
    genes_present <- genes_present[order(genes_present$ensembl_gene_id, decreasing = F), ]

    # subset expression matrix to only the genes with info on biomart
    expression_matrix <- expression_matrix[row.names(expression_matrix) %in% genes_present$ensembl_gene_id,]

    # deal with duplicate gene name/symbols (mostly ''s) by identifying the
    # duplicates and then creating a new gene name for them that is
    # genesymbol.ensemblID
    duplicates <- unique(genes_present$hgnc_symbol[duplicated(genes_present$hgnc_symbol)])
    for (k in duplicates) {
      dupl_genes <- which(genes_present$hgnc_symbol == k)
      for (m in 1:length(dupl_genes)) {
        genes_present$hgnc_symbol[dupl_genes[m]] <- paste(genes_present$hgnc_symbol[dupl_genes[m]],
          genes_present$ensembl_gene_id[dupl_genes[m]], sep = ".")
      }
    }
    row.names(expression_matrix) <- genes_present$hgnc_symbol
    rm(genes_present, duplicates)
    sc <- CreateSeuratObject(expression_matrix, min.cells = min_cells_per_gene, min.features = min_genes_per_cell)
    rm(expression_matrix)

    if (source_files == "lowpass") {
      demux <- read.table(paste0("/project2/gilad/reem/singlecellCM/round1/lowpass/CD",
      i, "/CD", i, "col", j, "/demux/CD", i, "col", j, "_demux.best"), header = T,
      stringsAsFactors = F)
    } else {
      demux <- read.table(paste0("/project2/gilad/reem/singlecellCM/round1/fulldata/CD",
      i, "/CD", i, "col", j, "/demux/hpCD", i, "col", j, "_demux.best"), header = T,
      stringsAsFactors = F)
    }

    m <- match(row.names(sc@meta.data), demux$BARCODE)
    if (any(is.na(m))) {
      cat(paste0("Not all barcodes are in demuxlet data. Something is wrong in CD",
        i, "col", j, "!\n"))
    }

    demux <- demux[m, ]

    demux$individual <- NA
    demux$individual[which(demux$BEST == "SNG-NA19093")] <- "19093"
    demux$individual[which(demux$BEST == "SNG-NA18912")] <- "18912"
    demux$individual[which(demux$BEST == "SNG-NA18858")] <- "18858"
    demux$individual[which(demux$BEST == "SNG-NA18520")] <- "18520"
    demux$individual[which(demux$BEST == "SNG-NA18511")] <- "18511"
    demux$individual[which(demux$BEST == "SNG-NA18508")] <- "18508"

    ind <- demux$individual
    names(ind) <- demux$BARCODE

    demux$diffday <- NA
    # col 1 for all CDs
    if (j == 1) {
      if (i == 1) {
        demux$diffday[which(demux$individual == "19093")] <- "Day 7"
        demux$diffday[which(demux$individual == "18912")] <- "Day 3"
        demux$diffday[which(demux$individual == "18520")] <- "Day 1"
      } else if (i == 2) {
        demux$diffday[which(demux$individual == "19093")] <- "Day 11"
        demux$diffday[which(demux$individual == "18912")] <- "Day 7"
        demux$diffday[which(demux$individual == "18520")] <- "Day 5"
      } else if (i == 3) {
        demux$diffday[which(demux$individual == "19093")] <- "Day 15"
        demux$diffday[which(demux$individual == "18912")] <- "Day 11"
        # add day 0 col1 here
        demux$diffday[which(demux$individual == "18858")] <- "Day 0"
      }
    }

    # col 2 for all CDs
    if (j == 2) {
      if (i == 1) {
        demux$diffday[which(demux$individual == "18858")] <- "Day 7"
        demux$diffday[which(demux$individual == "18520")] <- "Day 3"
        demux$diffday[which(demux$individual == "18508")] <- "Day 1"
      } else if (i == 2) {
        demux$diffday[which(demux$individual == "18858")] <- "Day 11"
        demux$diffday[which(demux$individual == "18520")] <- "Day 7"
        demux$diffday[which(demux$individual == "18508")] <- "Day 5"
      } else if (i == 3) {
        demux$diffday[which(demux$individual == "18858")] <- "Day 15"
        demux$diffday[which(demux$individual == "18520")] <- "Day 11"
        ## add day 0 col 2
        demux$diffday[which(demux$individual == "18912")] <- "Day 0"
      }
    }

    # col 3 for all CDs
    if (j == 3) {
      if (i == 1) {
        demux$diffday[which(demux$individual == "18912")] <- "Day 7"
        demux$diffday[which(demux$individual == "18508")] <- "Day 3"
        demux$diffday[which(demux$individual == "18511")] <- "Day 1"
      } else if (i == 2) {
        demux$diffday[which(demux$individual == "18912")] <- "Day 11"
        demux$diffday[which(demux$individual == "18508")] <- "Day 7"
        demux$diffday[which(demux$individual == "18511")] <- "Day 5"
      } else if (i == 3) {
        demux$diffday[which(demux$individual == "18912")] <- "Day 15"
        demux$diffday[which(demux$individual == "18508")] <- "Day 11"
        # add day 0 col 3
        demux$diffday[which(demux$individual == "18520")] <- "Day 0"
      }
    }

    # col 4 for all CDs
    if (j == 4) {
      if (i == 1) {
        demux$diffday[which(demux$individual == "18520")] <- "Day 7"
        demux$diffday[which(demux$individual == "18511")] <- "Day 3"
        demux$diffday[which(demux$individual == "19093")] <- "Day 1"
      } else if (i == 2) {
        demux$diffday[which(demux$individual == "18520")] <- "Day 11"
        demux$diffday[which(demux$individual == "18511")] <- "Day 7"
        demux$diffday[which(demux$individual == "19093")] <- "Day 5"
      } else if (i == 3) {
        demux$diffday[which(demux$individual == "18520")] <- "Day 15"
        demux$diffday[which(demux$individual == "18511")] <- "Day 11"
        # add day 0 col 4
        demux$diffday[which(demux$individual == "18508")] <- "Day 0"
      }
    }

    # col 5 for all CDs
    if (j == 5) {
      if (i == 1) {
        demux$diffday[which(demux$individual == "18508")] <- "Day 7"
        demux$diffday[which(demux$individual == "19093")] <- "Day 3"
        demux$diffday[which(demux$individual == "18858")] <- "Day 1"
      } else if (i == 2) {
        demux$diffday[which(demux$individual == "18508")] <- "Day 11"
        demux$diffday[which(demux$individual == "19093")] <- "Day 7"
        demux$diffday[which(demux$individual == "18858")] <- "Day 5"
      } else if (i == 3) {
        demux$diffday[which(demux$individual == "18508")] <- "Day 15"
        demux$diffday[which(demux$individual == "19093")] <- "Day 11"
        # add day 0 col 5
        demux$diffday[which(demux$individual == "18511")] <- "Day 0"
      }
    }

    # col 6 for all CDs
    if (j == 6) {
      if (i == 1) {
        demux$diffday[which(demux$individual == "18511")] <- "Day 7"
        demux$diffday[which(demux$individual == "18858")] <- "Day 3"
        demux$diffday[which(demux$individual == "18912")] <- "Day 1"
      } else if (i == 2) {
        demux$diffday[which(demux$individual == "18511")] <- "Day 11"
        demux$diffday[which(demux$individual == "18858")] <- "Day 7"
        demux$diffday[which(demux$individual == "18912")] <- "Day 5"
      } else if (i == 3) {
        demux$diffday[which(demux$individual == "18511")] <- "Day 15"
        demux$diffday[which(demux$individual == "18858")] <- "Day 11"
        # add day 0 col 6
        demux$diffday[which(demux$individual == "19093")] <- "Day 0"
      }
    }

    dday <- demux$diffday
    names(dday) <- demux$BARCODE

    sc <- AddMetaData(sc, ind, col.name = "individual")
    sc <- AddMetaData(sc, dday, col.name = "diffday")

    assign(paste0("CD", i, "col", j, "_lbld"), sc)
    rm(sc)
  }
}

rm(i, j, m, dday, ind, demux, geneinfo)

# merge AFTER giving them individual and diffday labels

sc <- merge(CD1col1_lbld, y = c(CD1col2_lbld, CD1col3_lbld, CD1col4_lbld,
  CD1col5_lbld, CD1col6_lbld, CD2col1_lbld, CD2col2_lbld, CD2col3_lbld, CD2col4_lbld,
  CD2col5_lbld, CD2col6_lbld, CD3col1_lbld, CD3col2_lbld, CD3col3_lbld, CD3col4_lbld,
  CD3col5_lbld, CD3col6_lbld), add.cell.ids = c("CD1col1", "CD1col2", "CD1col3",
  "CD1col4", "CD1col5", "CD1col6", "CD2col1", "CD2col2", "CD2col3", "CD2col4",
  "CD2col5", "CD2col6", "CD3col1", "CD3col2", "CD3col3", "CD3col4", "CD3col5",
  "CD3col6"))

# Assign colday to the cells
sc$colday <- "colday"
sc$colday <- substr(sc$orig.ident, 3, 3)

# get rid of high mito and subset to the cells demuxlet approved
sc$percent.mito <- PercentageFeatureSet(sc, pattern = "^MT-")
sc <- subset(sc, subset = percent.mito < mito_threshold)

diffdays <- c("Day 0", "Day 1", "Day 3", "Day 5", "Day 7", "Day 11", "Day 15")
sc <- subset(sc, subset=diffday %in% diffdays)
individuals <- c("19093", "18912", "18858", "18520", "18511", "18508")
sc <- subset(sc, subset=individual %in% individuals)
coldays <- c("CD1", "CD2", "CD3")
collections <- c("CD1col1", "CD1col2", "CD1col3", "CD1col4", "CD1col5", "CD1col6",
                       "CD2col1", "CD2col2", "CD2col3", "CD2col4", "CD2col5", "CD2col6",
                       "CD3col1", "CD3col2", "CD3col3", "CD3col4", "CD3col5", "CD3col6")


# Make sure that all metadata are in factor form for later analysis
sc$colday <- factor(x=sc$colday, levels=coldays, ordered=T)
sc$diffday <- factor(x=sc$diffday, levels=diffdays, ordered=T)
sc$individual <- factor(x=sc$individual, levels=individuals, ordered=T)
sc$orig.ident <- factor(x=sc$orig.ident, levels=collections, ordered=T)


# (Josh) Save to rds file
saveRDS(sc, file = "./rds_objects/old_seurat_obj_lowpass.rds")
