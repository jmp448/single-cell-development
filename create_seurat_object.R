library(stringr)
library(Seurat)
library(gtools)

### The following is largely copied directly from Reem Elorbany's work found at
### https://relorbany.github.io/Exp1_prelim_analysis/scCM_round1_lowpass_combined.html

load_from_bash <- TRUE

if (load_from_bash) {
  # If using the bash files
  args = commandArgs(trailingOnly = TRUE)
  rawdata <- c()
  for (i in 1:18) {
    rawdata <- c(rawdata, args[i])
  }

  min_cells_per_gene = args[19]  # minimum num cells in which a gene must appear
  min_genes_per_cell = args[20]  # minimum num genes for a cell to be included
  cutoff_mito = args[21]  # wanna cut off cells w a certain percent mito? (bool)
  mito_threshold = args[22]  # what is the threshold for mito cutoff?
} else {
  # If running this file directly in R, ie not using the bash files
  rawdata <- c('/project2/gilad/reem/singlecellCM/round1/fulldata/CD1/CD1col1/output/dge_data/YG-RE-RE1-hpCD1col1_S1_gene_counts.tsv.gz',
                '/project2/gilad/reem/singlecellCM/round1/fulldata/CD1/CD1col2/output/dge_data/YG-RE-RE2-hpCD1col2_S1_gene_counts.tsv.gz',
                '/project2/gilad/reem/singlecellCM/round1/fulldata/CD1/CD1col3/output/dge_data/YG-RE-RE3-hpCD1col3_S1_gene_counts.tsv.gz',
                '/project2/gilad/reem/singlecellCM/round1/fulldata/CD1/CD1col4/output/dge_data/YG-RE-RE4-hpCD1col4_S1_gene_counts.tsv.gz',
                '/project2/gilad/reem/singlecellCM/round1/fulldata/CD1/CD1col5/output/dge_data/YG-RE-RE5-hpCD1col5_S1_gene_counts.tsv.gz',
                '/project2/gilad/reem/singlecellCM/round1/fulldata/CD1/CD1col6/output/dge_data/YG-RE-RE6-hpCD1col6_S1_gene_counts.tsv.gz',
                '/project2/gilad/reem/singlecellCM/round1/fulldata/CD2/CD2col1/output/dge_data/YG-RE-RE3-hpCD2col1_S2_gene_counts.tsv.gz',
                '/project2/gilad/reem/singlecellCM/round1/fulldata/CD2/CD2col2/output/dge_data/YG-RE-RE4-hpCD2col2_S2_gene_counts.tsv.gz',
                '/project2/gilad/reem/singlecellCM/round1/fulldata/CD2/CD2col3/output/dge_data/YG-RE-RE5-hpCD2col3_S2_gene_counts.tsv.gz',
                '/project2/gilad/reem/singlecellCM/round1/fulldata/CD2/CD2col4/output/dge_data/YG-RE-RE6-hpCD2col4_S2_gene_counts.tsv.gz',
                '/project2/gilad/reem/singlecellCM/round1/fulldata/CD2/CD2col5/output/dge_data/YG-RE-RE2-hpCD2col5_S2_gene_counts.tsv.gz',
                '/project2/gilad/reem/singlecellCM/round1/fulldata/CD2/CD2col6/output/dge_data/YG-RE-RE1-hpCD2col6_S2_gene_counts.tsv.gz',
                '/project2/gilad/reem/singlecellCM/round1/fulldata/CD3/CD3col1/output/dge_data/YG-RE-RE4-hpCD3col1_S3_gene_counts.tsv.gz',
                '/project2/gilad/reem/singlecellCM/round1/fulldata/CD3/CD3col2/output/dge_data/YG-RE-RE3-hpCD3col2_S3_gene_counts.tsv.gz',
                '/project2/gilad/reem/singlecellCM/round1/fulldata/CD3/CD3col3/output/dge_data/YG-RE-RE2-CD3col3_Unk1_S4_gene_counts.tsv.gz',
                '/project2/gilad/reem/singlecellCM/round1/fulldata/CD3/CD3col4/output/dge_data/YG-RE-RE1-hpCD3col4_S3_gene_counts.tsv.gz',
                '/project2/gilad/reem/singlecellCM/round1/fulldata/CD3/CD3col5/output/dge_data/YG-RE-RE6-hpCD3col5_S3_gene_counts.tsv.gz',
                '/project2/gilad/reem/singlecellCM/round1/fulldata/CD3/CD3col6/output/dge_data/YG-RE-RE5-hpCD3col6_S3_gene_counts.tsv.gz')
  min_cells_per_gene = 3  # minimum num cells in which a gene must appear
  min_genes_per_cell = 200  # minimum num genes for a cell to be included
  cutoff_mito = FALSE  # wanna cut off cells w a certain percent mito? (bool)
  mito_threshold = 30  # what is the threshold for mito cutoff?
}


## Combine all 18 collections into a single Seurat object
## First, read in the raw data for all collections

## Get gene ID info from biomart
geneinfo <- readRDS("./rds_objects/geneinfo.rds")

for (i in 1:3) {
  for (j in 1:6) {
    # read in raw data matrix
    expression_matrix <- read.table(rawdata[6 * (i - 1) + j], header = T, stringsAsFactors = F,
      row.names = 1)

    # remove version numbers from those gene IDs
    rownames(expression_matrix) <- str_replace(rownames(expression_matrix), pattern = ".[0-9]+$",
      replacement = "")
    # subset geneinfo to only those genes from the full list of genes
    genes_present <- geneinfo[geneinfo$ensembl_gene_id %in% rownames(expression_matrix),]
    # subset expression matrix to only those genes with info on biomart
    expression_matrix <- expression_matrix[rownames(expression_matrix) %in% genes_present$ensembl_gene_id,]

    # sort genes by ensemblID (as they already are in genes_present)
    expression_matrix <- expression_matrix[order(mygeneinfo_temp$ensembl_gene_id, decreasing = F),]
    # get rid of duplicate ensembl IDs
    expression_matrix <- expression_matrix[!duplicated(rownames(expression_matrix)),]

    if !identical(rownames(expression_matrix), genes_present$ensemblID) {
      print("Look at line 75: expression matrix rownames and genes_present don't match")
      print("Here's the beginning of expression matrix:")
      print(head(expression_matrix))
      print("Here's the beginning of the genes_present ensemblIDs:")
      print(head(genes_present$ensemblID))
    }

    # deal with duplicate gene name/symbols (mostly ''s) by identifying the
    # duplicates and then creating a new gene name for them that is
    # genesymbol.ensemblID
    duplicates <- unique(genes_present$hgnc_symbol[duplicated(genes_present$hgnc_symbol)])
    for (symb in duplicates) {
      dupl_genes <- which(genes_present$hgnc_symbol == symb)
      for (g in 1:length(dupl_genes)) {
        genes_present$hgnc_symbol[dupl_genes[g]] <- paste0(genes_present$hgnc_symbol[dupl_genes[g]],
          ".", genes_present$ensembl_gene_id[dupl_genes[g]])
      }
    }
    # add those new gene names/symbols as rownames in the raw data matrix
    rownames(expression_matrix) <- genes_present$hgnc_symbol
    rm(genes_present)

    seurat_cd <- CreateSeuratObject(expression_matrix, min.cells = min_cells_per_gene,
      min.features = min_genes_per_cell, project = paste0("CD", i, "col", j))

    assign(paste0("Seurat_CD", i, "col", j), seurat_cd)

    rm(expression_matrix, seurat_cd)
  }
}

# Get rid of very large geneinfo from workspace
rm(geneinfo)

for (i in 1:3) {
  for (j in 1:6) {
    seurat_cd <- eval(as.name(paste0("Seurat_CD", i, "col", j)))

    # Get demuxlet data for this collection day
    demux <- read.table(paste0("/project2/gilad/reem/singlecellCM/round1/fulldata/CD",
      i, "/CD", i, "col", j, "/demux/hpCD", i, "col", j, "_demux.best"), header = T,
      stringsAsFactors = F)

    # Make sure all cells are represented in the demuxlet data
    m <- match(rownames(seurat_cd@meta.data), demux$BARCODE)
    if (any(is.na(m))) {
      cat(paste0("Not all barcodes are in demuxlet data. Something is wrong in CD",
        i, "col", j, "!\n"))
    }

    # Order the demuxlet data according to the matches
    demux <- demux[m, ]

    # Assign individual based on demuxlet assignments
    demux$individual <- NA
    demux$individual[which(demux$BEST == "SNG-NA19093")] <- "19093"
    demux$individual[which(demux$BEST == "SNG-NA18912")] <- "18912"
    demux$individual[which(demux$BEST == "SNG-NA18858")] <- "18858"
    demux$individual[which(demux$BEST == "SNG-NA18520")] <- "18520"
    demux$individual[which(demux$BEST == "SNG-NA18511")] <- "18511"
    demux$individual[which(demux$BEST == "SNG-NA18508")] <- "18508"

    individual_metadata <- demux$individual
    names(individual_metadata) <- demux$BARCODE

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

    diffday_metadata <- demux$diffday
    names(diffday_metadata) <- demux$BARCODE

    seurat_cd <- AddMetaData(seurat_cd, individual_metadata, col.name = "individual")
    seurat_cd <- AddMetaData(seurat_cd, diffday_metadata, col.name = "diffday")

    assign(paste0("Seurat_CD", i, "col", j), seurat_cd)
    rm(seurat_cd)
  }
}

# merge AFTER giving them individual and diffday labels
seurat_comb <- merge(Seurat_CD1col1, y = c(Seurat_CD1col2, Seurat_CD1col3,
  Seurat_CD1col4, Seurat_CD1col5, Seurat_CD1col6, Seurat_CD2col1, Seurat_CD2col2,
  Seurat_CD2col3, Seurat_CD2col4, Seurat_CD2col5, Seurat_CD2col6, Seurat_CD3col1,
  Seurat_CD3col2, Seurat_CD3col3, Seurat_CD3col4, Seurat_CD3col5, Seurat_CD3col6),
  add.cell.ids = c("CD1col1", "CD1col2", "CD1col3",
  "CD1col4", "CD1col5", "CD1col6", "CD2col1", "CD2col2", "CD2col3", "CD2col4",
  "CD2col5", "CD2col6", "CD3col1", "CD3col2", "CD3col3", "CD3col4", "CD3col5",
  "CD3col6"), project = "seurat_obj_fulldata")

# remove the individual Seurat objects for each col
rm(list = ls(pattern = "CD"))

# subset so only the cells with an individual and a diffday are used
# make diffday the relevant identity marker
Idents(seurat_comb) <- "diffday"
seurat_comb <- subset(seurat_comb, idents = c("Day 1", "Day 3", "Day 5", "Day 7",
  "Day 11", "Day 15", "Day 0"))

# make a 'sample' label in the metadata tag, which has the combined
# individual and diffday information in the form: NAInd.Day# make a sample label
seurat_comb$sample <- "collection"

seurat_comb$sample <- paste0(seurat_comb$individual, ".", seurat_comb$diffday)

# i will also add metadata info for percent mitochondrial genes
seurat_comb[["percent.mito"]] <- PercentageFeatureSet(seurat_comb, pattern = "^MT-")

# first, let's order the levels in a more intuitive way for diffday
ordered_diffdays <- c("Day 0", "Day 1", "Day 3", "Day 5", "Day 7", "Day 11", "Day 15")
seurat_comb@meta.data$diffday <- factor(x = seurat_comb@meta.data$diffday,
  levels = ordered_diffdays)

# for individual
ordered_individuals <- c("19093", "18858", "18912", "18520", "18508", "18511")
seurat_comb@meta.data$individual <- factor(x = seurat_comb@meta.data$individual,
  levels = ordered_individuals)

# for sample
ordered_samples <- unique(as.character(seurat_comb$sample))
ordered_samples <- mixedsort(ordered_samples, decreasing = F)
seurat_comb@meta.data$sample <- factor(seurat_comb@meta.data$sample, levels = ordered_samples)

# Assign colday to the cells
seurat_comb$colday <- "colday"
seurat_comb$colday <- substr(all_cols_noNA_S$orig.ident, 3, 3)

if (cutoff_mito) {
  seurat_comb <- subset(seurat_comb, subset = percent.mito < mito_threshold)
}

saveRDS(seurat_comb, "./rds_objects/seurat_obj_fulldata_myversion.RDS")
