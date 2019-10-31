library(biomaRt)
library(stringr)
library(Seurat)

### The following is, until the end, copied directly from Reem Elorbany's work
### found at
### https://relorbany.github.io/Exp1_prelim_analysis/scCM_round1_lowpass_combined.html

## My goal is to combine all 18 collections into a single Seurat object.  First, I
## need to read in the raw data for all collections.
load_from_bash <- F
if (load_from_bash) {
  args = commandArgs(trailingOnly = TRUE)
  rawdata <- c()
  for (i in 1:18) {
    rawdata <- c(rawdata, args[i])
  }
  min_cells_per_gene <- args[19]  # minimum num cells in which a gene must appear
  min_genes_per_cell <- args[20]  # minimum num genes for a cell to be included
  mito_threshold <- args[21]  # what is the threshold for mito cutoff?
} else {
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

for (i in 1:1) {
  for (j in 1:1) {
    # read in every raw data matrix
    raw_dat_temp <- read.table(rawdata[6*(i-1)+j], header = T, stringsAsFactors = F, row.names = 1)

    # remove the version numbers from those gene IDs
    gene_id_temp <- str_replace(rownames(raw_dat_temp), pattern = ".[0-9]+$",
      replacement = "")
    # and subset to only those genes from the full list of genes
    mygeneinfo_temp <- geneinfo[geneinfo$ensembl_gene_id %in% gene_id_temp, ]

    # make a temporary raw data matrix with the non-version gene IDs
    bm_rawdat_temp <- raw_dat_temp
    rownames(bm_rawdat_temp) <- gene_id_temp
    # subset temp raw data matrix to only the ones with info on biomart
    bm_rawdat_temp <- bm_rawdat_temp[rownames(bm_rawdat_temp) %in% mygeneinfo_temp$ensembl_gene_id,
      ]

    # get all the index numbers for a gene list sorted by number
    mygeneinfo_idsort_temp <- order(mygeneinfo_temp$ensembl_gene_id, decreasing = F)
    # and sort the genes with those index numbers
    mygeneinfo_sort_temp <- mygeneinfo_temp[mygeneinfo_idsort_temp, ]

    # get rid of duplicate ensembl IDs
    mygeneinfo_s_uni_temp <- mygeneinfo_sort_temp[order(mygeneinfo_sort_temp$ensembl_gene_id)[!duplicated(mygeneinfo_sort_temp$ensembl_gene_id)],
      ]

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

    # make a gene info matrix with the now unique gene symbols
    mygeneinfo_gsym_temp <- mygeneinfo_s_uni_temp[order(mygeneinfo_s_uni_temp$hgnc_symbol)[!duplicated(mygeneinfo_s_uni_temp$hgnc_symbol)],
      ]

    # add those new gene names/symbols as rownames in the raw data matrix
    rownames(bm_rawdat_temp) <- mygeneinfo_s_uni_temp$hgnc_symbol

    # output the raw data matrix with the proper collection name
    assign(paste0("bm_rawdat_C", i, "c", j), bm_rawdat_temp)

  }
}

# remove all those temp variables so it doesn't clutter my env
rm(i, j, k, m, rows_temp, mygeneinfo_gsym_temp, mygeneinfo_idsort_temp, mygeneinfo_sort_temp,
  mygeneinfo_s_uni_temp, mygeneinfo_temp, gene_id_temp, raw_dat_temp, dupl_temp,
  bm_rawdat_temp)

## Now let's create the individual Seurat objects for each collection.  I do this
## first so I can assign individual and diffday labels more easily, before merging
## and subsetting the object.

for (i in 1:3) {
  for (j in 1:6) {

    bmrawdat_temp <- eval(as.name(paste0("bm_rawdat_C", i, "c", j)))

    SObject_temp <- CreateSeuratObject(bmrawdat_temp, min.cells = 3, min.features = 200,
      project = paste0("CD", i, "col", j))

    assign(paste0("CD", i, "col", j, "SObj"), SObject_temp)

  }
}

rm(list = ls(pattern = "bm_rawdat"))
rm(bmrawdat_temp, SObject_temp, i, j)

# these are Seurat objects for each collection with gene symbols!  filtered by
# only genes in at least 3 cells and only cells with at least 200 genes

## The next step is assigning labels for 'individual' and 'diffday' to each cell
## barcode.  I do this using the results from demuxlet, which was run previously
## on the data from each collection.

for (i in 1:3) {
  for (j in 1:6) {

    SObject <- eval(as.name(paste0("CD", i, "col", j, "SObj")))

    demux_temp <- read.table(paste0("/project2/gilad/reem/singlecellCM/round1/lowpass/CD",
      i, "/CD", i, "col", j, "/demux/CD", i, "col", j, "_demux.best"), header = T,
      stringsAsFactors = F)

    # demux_temp <- demux_temp[-1,]

    m <- match(rownames(SObject@meta.data), demux_temp$BARCODE)
    if (any(is.na(m)))
      cat(paste0("Not all barcodes are in demuxlet data. Something is wrong in CD",
        i, "col", j, "!\n"))

    # demux_temp <- demux_temp[m, ]
    #
    # demux_temp$individual <- "doublet"
    # demux_temp$individual[which(demux_temp$BEST == "SNG-NA19093")] <- "NA19093"
    # demux_temp$individual[which(demux_temp$BEST == "SNG-NA18912")] <- "NA18912"
    # demux_temp$individual[which(demux_temp$BEST == "SNG-NA18858")] <- "NA18858"
    # demux_temp$individual[which(demux_temp$BEST == "SNG-NA18520")] <- "NA18520"
    # demux_temp$individual[which(demux_temp$BEST == "SNG-NA18511")] <- "NA18511"
    # demux_temp$individual[which(demux_temp$BEST == "SNG-NA18508")] <- "NA18508"
    #
    # tmp_ind <- demux_temp$individual
    # names(tmp_ind) <- demux_temp$BARCODE
    #
    # demux_temp$diffday <- "NA"
    #
    # # col 1 for all CDs
    # if (j == 1) {
    #   if (i == 1) {
    #     demux_temp$diffday[which(demux_temp$individual == "NA19093")] <- "Day 7"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18912")] <- "Day 3"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18520")] <- "Day 1"
    #   } else if (i == 2) {
    #     demux_temp$diffday[which(demux_temp$individual == "NA19093")] <- "Day 11"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18912")] <- "Day 7"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18520")] <- "Day 5"
    #   } else if (i == 3) {
    #     demux_temp$diffday[which(demux_temp$individual == "NA19093")] <- "Day 15"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18912")] <- "Day 11"
    #     # add day 0 col1 here
    #     demux_temp$diffday[which(demux_temp$individual == "NA18858")] <- "Day 0"
    #   }
    # }
    #
    # # col 2 for all CDs
    # if (j == 2) {
    #   if (i == 1) {
    #     demux_temp$diffday[which(demux_temp$individual == "NA18858")] <- "Day 7"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18520")] <- "Day 3"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18508")] <- "Day 1"
    #   } else if (i == 2) {
    #     demux_temp$diffday[which(demux_temp$individual == "NA18858")] <- "Day 11"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18520")] <- "Day 7"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18508")] <- "Day 5"
    #   } else if (i == 3) {
    #     demux_temp$diffday[which(demux_temp$individual == "NA18858")] <- "Day 15"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18520")] <- "Day 11"
    #     ## add day 0 col 2
    #     demux_temp$diffday[which(demux_temp$individual == "NA18912")] <- "Day 0"
    #   }
    # }
    #
    # # col 3 for all CDs
    # if (j == 3) {
    #   if (i == 1) {
    #     demux_temp$diffday[which(demux_temp$individual == "NA18912")] <- "Day 7"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18508")] <- "Day 3"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18511")] <- "Day 1"
    #   } else if (i == 2) {
    #     demux_temp$diffday[which(demux_temp$individual == "NA18912")] <- "Day 11"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18508")] <- "Day 7"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18511")] <- "Day 5"
    #   } else if (i == 3) {
    #     demux_temp$diffday[which(demux_temp$individual == "NA18912")] <- "Day 15"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18508")] <- "Day 11"
    #     # add day 0 col 3
    #     demux_temp$diffday[which(demux_temp$individual == "NA18520")] <- "Day 0"
    #   }
    # }
    #
    # # col 4 for all CDs
    # if (j == 4) {
    #   if (i == 1) {
    #     demux_temp$diffday[which(demux_temp$individual == "NA18520")] <- "Day 7"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18511")] <- "Day 3"
    #     demux_temp$diffday[which(demux_temp$individual == "NA19093")] <- "Day 1"
    #   } else if (i == 2) {
    #     demux_temp$diffday[which(demux_temp$individual == "NA18520")] <- "Day 11"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18511")] <- "Day 7"
    #     demux_temp$diffday[which(demux_temp$individual == "NA19093")] <- "Day 5"
    #   } else if (i == 3) {
    #     demux_temp$diffday[which(demux_temp$individual == "NA18520")] <- "Day 15"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18511")] <- "Day 11"
    #     # add day 0 col 4
    #     demux_temp$diffday[which(demux_temp$individual == "NA18508")] <- "Day 0"
    #   }
    # }
    #
    # # col 5 for all CDs
    # if (j == 5) {
    #   if (i == 1) {
    #     demux_temp$diffday[which(demux_temp$individual == "NA18508")] <- "Day 7"
    #     demux_temp$diffday[which(demux_temp$individual == "NA19093")] <- "Day 3"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18858")] <- "Day 1"
    #   } else if (i == 2) {
    #     demux_temp$diffday[which(demux_temp$individual == "NA18508")] <- "Day 11"
    #     demux_temp$diffday[which(demux_temp$individual == "NA19093")] <- "Day 7"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18858")] <- "Day 5"
    #   } else if (i == 3) {
    #     demux_temp$diffday[which(demux_temp$individual == "NA18508")] <- "Day 15"
    #     demux_temp$diffday[which(demux_temp$individual == "NA19093")] <- "Day 11"
    #     # add day 0 col 5
    #     demux_temp$diffday[which(demux_temp$individual == "NA18511")] <- "Day 0"
    #   }
    # }
    #
    # # col 6 for all CDs
    # if (j == 6) {
    #   if (i == 1) {
    #     demux_temp$diffday[which(demux_temp$individual == "NA18511")] <- "Day 7"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18858")] <- "Day 3"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18912")] <- "Day 1"
    #   } else if (i == 2) {
    #     demux_temp$diffday[which(demux_temp$individual == "NA18511")] <- "Day 11"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18858")] <- "Day 7"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18912")] <- "Day 5"
    #   } else if (i == 3) {
    #     demux_temp$diffday[which(demux_temp$individual == "NA18511")] <- "Day 15"
    #     demux_temp$diffday[which(demux_temp$individual == "NA18858")] <- "Day 11"
    #     # add day 0 col 6
    #     demux_temp$diffday[which(demux_temp$individual == "NA19093")] <- "Day 0"
    #   }
    # }
    #
    # tmp_dday <- demux_temp$diffday
    # names(tmp_dday) <- demux_temp$BARCODE
    #
    # SObject <- AddMetaData(SObject, tmp_ind, col.name = "individual")
    # SObject <- AddMetaData(SObject, tmp_dday, col.name = "diffday")
    #
    # assign(paste0("CD", i, "col", j, "_lbld"), SObject)
  }
}

# rm(i, j, m, tmp_dday, tmp_ind, SObject, demux_temp)
# rm(list = ls(pattern = "SObj"))
#
# # this is a mess lol but i'm leaving it until I need to find a better way
#
# # if we want to explore the data, can do: table(CD1col1_lbld$individual)
# # table(CD1col1_lbld$diffday)
#
# # merge AFTER giving them individual and diffday labels
#
# all_cols_S <- merge(CD1col1_lbld, y = c(CD1col2_lbld, CD1col3_lbld, CD1col4_lbld,
#   CD1col5_lbld, CD1col6_lbld, CD2col1_lbld, CD2col2_lbld, CD2col3_lbld, CD2col4_lbld,
#   CD2col5_lbld, CD2col6_lbld, CD3col1_lbld, CD3col2_lbld, CD3col3_lbld, CD3col4_lbld,
#   CD3col5_lbld, CD3col6_lbld), add.cell.ids = c("CD1col1", "CD1col2", "CD1col3",
#   "CD1col4", "CD1col5", "CD1col6", "CD2col1", "CD2col2", "CD2col3", "CD2col4",
#   "CD2col5", "CD2col6", "CD3col1", "CD3col2", "CD3col3", "CD3col4", "CD3col5",
#   "CD3col6"), project = "CMdiff_round1_lowpass")
#
# # Assign colday to the cells
# all_cols_S$colday <- "colday"
# all_cols_S$colday <- substr(all_cols_S$orig.ident, 3, 3)
# table(all_cols_S$colday)
#
# # get rid of high mito
# all_cols_S[["percent.mito"]] <- PercentageFeatureSet(all_cols_S, pattern = "^MT-")
# all_cols_S <- subset(all_cols_S, subset = percent.mito < 30)
#
# # (Josh) Save to rds file
# saveRDS(all_cols_S, file = "./rds_objects/old_seurat_obj_lowpass.rds")
