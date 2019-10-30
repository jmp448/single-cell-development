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
  mito_threshold = args[21]  # what is the threshold for mito cutoff?
  source_files = args[22]
} else {
  # If running this file directly in R, ie not using the bash files, here's the full data
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
  # cutoff_mito = FALSE  # wanna cut off cells w a certain percent mito? (bool)
  mito_threshold = 30  # what is the threshold for mito cutoff?
}


## Combine all 18 collections into a single Seurat object
## First, read in the raw data for all collections

## Get gene ID info from biomart
geneinfo <- readRDS("./rds_objects/geneinfo.rds")

for (i in 1:3) {
  for (j in 1:6) {
    # read in every raw data matrix
    raw_dat_temp <- read.table(rawdata[6 * (i - 1) + j], header = T, stringsAsFactors = F,
      row.names = 1)

    # remove the version numbers from those gene IDs
    gene_id_temp <- str_replace(rownames(raw_dat_temp), pattern = ".[0-9]+$",
      replacement = "")
    # and subset to only those genes from the full list of genes
    mygeneinfo_temp <- geneinfo[geneinfo$ensembl_gene_id %in% gene_id_temp, ]

    # make a temporary raw data matrix with the non-version gene IDs
    bm_rawdat_temp <- raw_dat_temp
    rownames(bm_rawdat_temp) <- gene_id_temp
    rm(gene_id_temp, raw_dat_temp)
    # subset temp raw data matrix to only the ones with info on biomart
    bm_rawdat_temp <- bm_rawdat_temp[rownames(bm_rawdat_temp) %in% mygeneinfo_temp$ensembl_gene_id,
      ]

    # get all the index numbers for a gene list sorted by number
    mygeneinfo_idsort_temp <- order(mygeneinfo_temp$ensembl_gene_id, decreasing = F)
    # and sort the genes with those index numbers
    mygeneinfo_sort_temp <- mygeneinfo_temp[mygeneinfo_idsort_temp, ]
    rm(mygeneinfo_temp, mygeneinfo_idsort_temp)
    # get rid of duplicate ensembl IDs
    mygeneinfo_s_uni_temp <- mygeneinfo_sort_temp[order(mygeneinfo_sort_temp$ensembl_gene_id)[!duplicated(mygeneinfo_sort_temp$ensembl_gene_id)],
      ]
    rm(mygeneinfo_sort_temp)
    # deal with duplicate gene name/symbols (mostly ''s) by identifying the
    # duplicates and then creating a new gene name for them that is
    # genesymbol.ensemblID
    dupl_temp <- unique(mygeneinfo_s_uni_temp$hgnc_symbol[duplicated(mygeneinfo_s_uni_temp$hgnc_symbol)])
    for (k in dupl_temp) {
      rows_temp <- which(mygeneinfo_s_uni_temp$hgnc_symbol == k)
      for (m in 1:length(rows_temp)) {
        mygeneinfo_s_uni_temp$hgnc_symbol[rows_temp[m]] <- paste(mygeneinfo_s_uni_temp$hgnc_symbol[rows_temp[m]],
          ".", mygeneinfo_s_uni_temp$ensembl_gene_id[rows_temp[m]], sep = "")
      }
    }
    # add those new gene names/symbols as rownames in the raw data matrix
    rownames(bm_rawdat_temp) <- mygeneinfo_s_uni_temp$hgnc_symbol
    rm(mygeneinfo_s_uni_temp)

    SObject_temp <- CreateSeuratObject(bm_rawdat_temp, min.cells = min_cells_per_gene,
      min.features = min_genes_per_cell, project = paste0("CD", i, "col", j))

    assign(paste0("CD", i, "col", j, "SObj"), SObject_temp)

    rm(bm_rawdat_temp, SObject_temp)
  }
}

# Get rid of very large geneinfo from workspace
rm(geneinfo)

for (i in 1:3) {
  for (j in 1:6) {
    SObject <- eval(as.name(paste0("CD", i, "col", j, "SObj")))

    if (source_files == "fulldata") {
      demux_temp <- read.table(paste0("/project2/gilad/reem/singlecellCM/round1/fulldata/CD",
        i, "/CD", i, "col", j, "/demux/hpCD", i, "col", j, "_demux.best"), header = T,
        stringsAsFactors = F)
    } else if (source_files == "lowpass") {
      demux_temp <- read.table(paste0("/project2/gilad/reem/singlecellCM/round1/lowpass/CD",
        i, "/CD", i, "col", j, "/demux/CD", i, "col", j, "_demux.best"), header = T,
        stringsAsFactors = F)
    }

    m <- match(rownames(SObject@meta.data), demux_temp$BARCODE)
    if (any(is.na(m))) {
      cat(paste0("Not all barcodes are in demuxlet data. Something is wrong in CD",
        i, "col", j, "!\n"))
    }

    demux_temp <- demux_temp[m, ]

    demux_temp$individual <- "doublet"
    demux_temp$individual[which(demux_temp$BEST == "SNG-NA19093")] <- "NA19093"
    demux_temp$individual[which(demux_temp$BEST == "SNG-NA18912")] <- "NA18912"
    demux_temp$individual[which(demux_temp$BEST == "SNG-NA18858")] <- "NA18858"
    demux_temp$individual[which(demux_temp$BEST == "SNG-NA18520")] <- "NA18520"
    demux_temp$individual[which(demux_temp$BEST == "SNG-NA18511")] <- "NA18511"
    demux_temp$individual[which(demux_temp$BEST == "SNG-NA18508")] <- "NA18508"

    tmp_ind <- demux_temp$individual
    names(tmp_ind) <- demux_temp$BARCODE

    demux_temp$diffday <- "NA"

    # col 1 for all CDs
    if (j == 1) {
      if (i == 1) {
        demux_temp$diffday[which(demux_temp$individual == "NA19093")] <- "Day 7"
        demux_temp$diffday[which(demux_temp$individual == "NA18912")] <- "Day 3"
        demux_temp$diffday[which(demux_temp$individual == "NA18520")] <- "Day 1"
      } else if (i == 2) {
        demux_temp$diffday[which(demux_temp$individual == "NA19093")] <- "Day 11"
        demux_temp$diffday[which(demux_temp$individual == "NA18912")] <- "Day 7"
        demux_temp$diffday[which(demux_temp$individual == "NA18520")] <- "Day 5"
      } else if (i == 3) {
        demux_temp$diffday[which(demux_temp$individual == "NA19093")] <- "Day 15"
        demux_temp$diffday[which(demux_temp$individual == "NA18912")] <- "Day 11"
        # add day 0 col1 here
        demux_temp$diffday[which(demux_temp$individual == "NA18858")] <- "Day 0"
      }
    }

    # col 2 for all CDs
    if (j == 2) {
      if (i == 1) {
        demux_temp$diffday[which(demux_temp$individual == "NA18858")] <- "Day 7"
        demux_temp$diffday[which(demux_temp$individual == "NA18520")] <- "Day 3"
        demux_temp$diffday[which(demux_temp$individual == "NA18508")] <- "Day 1"
      } else if (i == 2) {
        demux_temp$diffday[which(demux_temp$individual == "NA18858")] <- "Day 11"
        demux_temp$diffday[which(demux_temp$individual == "NA18520")] <- "Day 7"
        demux_temp$diffday[which(demux_temp$individual == "NA18508")] <- "Day 5"
      } else if (i == 3) {
        demux_temp$diffday[which(demux_temp$individual == "NA18858")] <- "Day 15"
        demux_temp$diffday[which(demux_temp$individual == "NA18520")] <- "Day 11"
        ## add day 0 col 2
        demux_temp$diffday[which(demux_temp$individual == "NA18912")] <- "Day 0"
      }
    }

    # col 3 for all CDs
    if (j == 3) {
      if (i == 1) {
        demux_temp$diffday[which(demux_temp$individual == "NA18912")] <- "Day 7"
        demux_temp$diffday[which(demux_temp$individual == "NA18508")] <- "Day 3"
        demux_temp$diffday[which(demux_temp$individual == "NA18511")] <- "Day 1"
      } else if (i == 2) {
        demux_temp$diffday[which(demux_temp$individual == "NA18912")] <- "Day 11"
        demux_temp$diffday[which(demux_temp$individual == "NA18508")] <- "Day 7"
        demux_temp$diffday[which(demux_temp$individual == "NA18511")] <- "Day 5"
      } else if (i == 3) {
        demux_temp$diffday[which(demux_temp$individual == "NA18912")] <- "Day 15"
        demux_temp$diffday[which(demux_temp$individual == "NA18508")] <- "Day 11"
        # add day 0 col 3
        demux_temp$diffday[which(demux_temp$individual == "NA18520")] <- "Day 0"
      }
    }

    # col 4 for all CDs
    if (j == 4) {
      if (i == 1) {
        demux_temp$diffday[which(demux_temp$individual == "NA18520")] <- "Day 7"
        demux_temp$diffday[which(demux_temp$individual == "NA18511")] <- "Day 3"
        demux_temp$diffday[which(demux_temp$individual == "NA19093")] <- "Day 1"
      } else if (i == 2) {
        demux_temp$diffday[which(demux_temp$individual == "NA18520")] <- "Day 11"
        demux_temp$diffday[which(demux_temp$individual == "NA18511")] <- "Day 7"
        demux_temp$diffday[which(demux_temp$individual == "NA19093")] <- "Day 5"
      } else if (i == 3) {
        demux_temp$diffday[which(demux_temp$individual == "NA18520")] <- "Day 15"
        demux_temp$diffday[which(demux_temp$individual == "NA18511")] <- "Day 11"
        # add day 0 col 4
        demux_temp$diffday[which(demux_temp$individual == "NA18508")] <- "Day 0"
      }
    }

    # col 5 for all CDs
    if (j == 5) {
      if (i == 1) {
        demux_temp$diffday[which(demux_temp$individual == "NA18508")] <- "Day 7"
        demux_temp$diffday[which(demux_temp$individual == "NA19093")] <- "Day 3"
        demux_temp$diffday[which(demux_temp$individual == "NA18858")] <- "Day 1"
      } else if (i == 2) {
        demux_temp$diffday[which(demux_temp$individual == "NA18508")] <- "Day 11"
        demux_temp$diffday[which(demux_temp$individual == "NA19093")] <- "Day 7"
        demux_temp$diffday[which(demux_temp$individual == "NA18858")] <- "Day 5"
      } else if (i == 3) {
        demux_temp$diffday[which(demux_temp$individual == "NA18508")] <- "Day 15"
        demux_temp$diffday[which(demux_temp$individual == "NA19093")] <- "Day 11"
        # add day 0 col 5
        demux_temp$diffday[which(demux_temp$individual == "NA18511")] <- "Day 0"
      }
    }

    # col 6 for all CDs
    if (j == 6) {
      if (i == 1) {
        demux_temp$diffday[which(demux_temp$individual == "NA18511")] <- "Day 7"
        demux_temp$diffday[which(demux_temp$individual == "NA18858")] <- "Day 3"
        demux_temp$diffday[which(demux_temp$individual == "NA18912")] <- "Day 1"
      } else if (i == 2) {
        demux_temp$diffday[which(demux_temp$individual == "NA18511")] <- "Day 11"
        demux_temp$diffday[which(demux_temp$individual == "NA18858")] <- "Day 7"
        demux_temp$diffday[which(demux_temp$individual == "NA18912")] <- "Day 5"
      } else if (i == 3) {
        demux_temp$diffday[which(demux_temp$individual == "NA18511")] <- "Day 15"
        demux_temp$diffday[which(demux_temp$individual == "NA18858")] <- "Day 11"
        # add day 0 col 6
        demux_temp$diffday[which(demux_temp$individual == "NA19093")] <- "Day 0"
      }
    }

    tmp_dday <- demux_temp$diffday
    names(tmp_dday) <- demux_temp$BARCODE

    SObject <- AddMetaData(SObject, tmp_ind, col.name = "individual")
    SObject <- AddMetaData(SObject, tmp_dday, col.name = "diffday")

    assign(paste0("CD", i, "col", j, "_lbld"), SObject)
  }
}

# merge AFTER giving them individual and diffday labels
all_cols_S <- merge(CD1col1_lbld, y = c(CD1col2_lbld, CD1col3_lbld, CD1col4_lbld,
  CD1col5_lbld, CD1col6_lbld, CD2col1_lbld, CD2col2_lbld, CD2col3_lbld, CD2col4_lbld,
  CD2col5_lbld, CD2col6_lbld, CD3col1_lbld, CD3col2_lbld, CD3col3_lbld, CD3col4_lbld,
  CD3col5_lbld, CD3col6_lbld), add.cell.ids = c("CD1col1", "CD1col2", "CD1col3",
  "CD1col4", "CD1col5", "CD1col6", "CD2col1", "CD2col2", "CD2col3", "CD2col4",
  "CD2col5", "CD2col6", "CD3col1", "CD3col2", "CD3col3", "CD3col4", "CD3col5",
  "CD3col6"), project = "scCM_round1_fulldata")

# finally, i'll remove the individual Seurat objects for each col
rm(list = ls(pattern = "lbld"))

# now i want to subset so only the cells with an individual and a diffday are
# being used (i.e. no mislabeled individuals on each collection and no doublets)
Idents(all_cols_S) <- "diffday"
# make diffday the relevant identity marker

all_cols_noNA_S = subset(all_cols_S, idents = c("Day 1", "Day 3", "Day 5", "Day 7",
  "Day 11", "Day 15", "Day 0"))

# I will make a 'sample' label in the metadata tag, which has the combined
# individual and diffday information in the form: NAInd.Day# make a sample label

# let's make a column that's a combo of individual and diffday
all_cols_noNA_S$sample <- "IND.DAY"

all_cols_noNA_S$sample <- paste0(all_cols_noNA_S$individual, ".", all_cols_noNA_S$diffday)

# i will also add metadata info for percent mitochondrial genes
all_cols_noNA_S[["percent.mito"]] <- PercentageFeatureSet(all_cols_noNA_S, pattern = "^MT-")

# first, let's order the levels in a more intuitive way for diffday
mydaylevels <- c("Day 0", "Day 1", "Day 3", "Day 5", "Day 7", "Day 11", "Day 15")
all_cols_noNA_S@meta.data$diffday <- factor(x = all_cols_noNA_S@meta.data$diffday,
  levels = mydaylevels)
# for individual
myindlevels <- c("NA19093", "NA18858", "NA18912", "NA18520", "NA18508", "NA18511")
all_cols_noNA_S@meta.data$individual <- factor(x = all_cols_noNA_S@meta.data$individual,
  levels = myindlevels)
# for sample
mysamplelevels <- as.character(all_cols_noNA_S$sample)
mysamplelevels <- mixedsort(mysamplelevels, decreasing = F)
mysamplelevels <- unique(mysamplelevels)
all_cols_noNA_S@meta.data$sample <- factor(all_cols_noNA_S@meta.data$sample, levels = mysamplelevels)

# Assign colday to the cells
all_cols_noNA_S$colday <- "colday"
all_cols_noNA_S$colday <- substr(all_cols_noNA_S$orig.ident, 3, 3)

all_cols_noNA_S <- subset(all_cols_noNA_S, subset = percent.mito < mito_threshold)

saveRDS(all_cols_noNA_S, paste0("./rds_objects/seurat_obj_", source_files, ".RDS"))