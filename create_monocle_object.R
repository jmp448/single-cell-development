library(stringr)
library(Matrix)
library(monocle3)
library(tidyverse)

load_from_bash<-T

if (load_from_bash) {
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

update_masters <- function(combined_expression, combined_genes, combined_cells, expression_matrix, genes_matrix, cell_metadata) {
  # create new genes matrix
  new_genes <- subset(genes_matrix, !(rownames(genes_matrix) %in% rownames(combined_genes)))
  overlap_genes <- subset(rownames(genes_matrix), rownames(genes_matrix) %in% rownames(combined_genes))
  old_genes <- subset(rownames(combined_genes), !(rownames(combined_genes) %in% rownames(genes_matrix)))
  num_new_genes <- length(rownames(new_genes))
  num_overlap_genes <- length(overlap_genes)
  num_old_genes <- length(old_genes)
  combined_genes <- rbind(combined_genes, new_genes)

  # sort appropriately for the new genes
  combined_genes %>% arrange(rownames(combined_genes))

  # open up some space
  rm(genes_matrix)
  new_genes <- rownames(new_genes)

  # create new cells matrix
  num_new_cells <- length(rownames(cell_metadata))
  num_old_cells <- length(rownames(combined_cells))
  combined_cells <- rbind(combined_cells, cell_metadata)
  rm(cell_metadata)

  # add new genes into combined expression matrix
  combined_expression[(num_old_genes+num_overlap_genes+1):(num_old_genes+num_overlap_genes+num_new_genes),] <- 0
  rownames(combined_expression)[(num_old_genes+num_overlap_genes+1):(num_old_genes+num_overlap_genes+num_new_genes)] <- new_genes

  # add old genes into new expression matrix
  expression_matrix[(num_new_genes+num_overlap_genes+1):(num_new_genes+num_overlap_genes+num_old_genes),] <- 0
  rownames(expression_matrix)[(num_new_genes+num_overlap_genes+1):(num_new_genes+num_overlap_genes+num_old_genes)] <- old_genes

  # sort rows on both new and old before combination
  expression_matrix %>% arrange(rownames(expression_matrix))
  combined_expression %>% arrange(rownames(combined_expression))

  # combine expression matrices
  combined_expression <- cbind(combined_expression, expression_matrix)
  masterlist <- list(combined_expression, combined_genes, combined_cells)
  return(masterlist)
}

## First, read in the raw data for all collections

## Get gene ID info from biomart
geneinfo <- readRDS("./rds_objects/geneinfo.rds")

for (i in 1:3) {
  for (j in 1:6) {
    # Read in raw data matrix
    expression_matrix <- read.table(rawdata[6 * (i - 1) + j], header = T, stringsAsFactors = F, row.names = 1)

    # remove the version numbers from those gene IDs
    non_version_genes <- str_replace(rownames(expression_matrix), pattern = ".[0-9]+$", replacement = "")
    rownames(expression_matrix) <- non_version_genes

    # remove any duplicates (later versions, etc)
    expression_matrix <- expression_matrix[!duplicated(rownames(expression_matrix)),]

    # subset the genes from Biomart to those genes in the sample
    genes_present <- geneinfo[geneinfo$ensembl_gene_id %in% non_version_genes, ]
    genes_present <- genes_present[!duplicated(genes_present$ensembl_gene_id), ]

    # subset the genes in the sample to the genes in Biomart
    expression_matrix <- expression_matrix[rownames(expression_matrix) %in% genes_present$ensembl_gene_id, ]

    rm(non_version_genes)

    ### Sort genes by ensemblID
    id_sort <- order(genes_present$ensembl_gene_id, decreasing=F)
    genes_present <- genes_present[id_sort, ]
    rm(id_sort)

    ### Handle duplicates in gene names (due to unnamed transcripts)
    duplicates <- unique(genes_present$hgnc_symbol[duplicated(genes_present$hgnc_symbol) | duplicated(genes_present$hgnc_symbol, fromLast = T)])
    for (g in duplicates) {
      copies_of_g <- which(genes_present$hgnc_symbol == g)
      if (g == "") {
        for (m in 1:length(copies_of_g)) {
          genes_present$hgnc_symbol[copies_of_g[m]] <- paste0("unnamed.", genes_present$ensembl_gene_id[copies_of_g[m]])
        }
      } else {
        for (m in 1:length(copies_of_g)) {
            genes_present$hgnc_symbol[copies_of_g[m]] <- paste0(g, ".", genes_present$ensembl_gene_id[copies_of_g[m]])
        }
      }
    }
    rm(duplicates, copies_of_g, g, m)

    ### Filter for quality
    cell_counts <- rowSums(expression_matrix)
    expression_matrix <- expression_matrix[cell_counts >= min_cells_per_gene,]
    genes_present <- genes_present[cell_counts>=min_cells_per_gene,]

    gene_counts <- colSums(expression_matrix)
    expression_matrix <- expression_matrix[,gene_counts>=min_genes_per_cell]
    rm(cell_counts, gene_counts)

    ### Prep genes for CDS construction
    genes_matrix <- data.frame(genes_present$hgnc_symbol)
    rownames(genes_matrix) <- genes_present$ensembl_gene_id
    colnames(genes_matrix) <- "gene_short_name"
    rm(genes_present)

    ### Put together cell metadata
    demux_file <- paste0('/project2/gilad/reem/singlecellCM/round1/fulldata/CD', i, '/CD', i, 'col', j, '/demux/hpCD', i, 'col', j, '_demux.best')
    demux <- read.table(demux_file, header = T, stringsAsFactors = F, row.names = 1)

    m <- match(colnames(expression_matrix), rownames(demux))
    demux <- demux[m,]

    # Check that all cells were run through demuxlet
    if (any(is.na(m))) {
      cat(paste0("Not all barcodes are in demuxlet data. ", sum(is.na(m)), " barcodes are missing for CD", i, "col", j, "!\n"))
      expression_matrix <- expression_matrix[,!is.na(m)]
      demux <- demux[!is.na(m),]
    }
    rm(m)

    # Assign individuals to the metadata
    demux$individual <- "doublet"
    demux$individual[which(demux$BEST == "SNG-NA19093")] <- "NA19093"
    demux$individual[which(demux$BEST == "SNG-NA18912")] <- "NA18912"
    demux$individual[which(demux$BEST == "SNG-NA18858")] <- "NA18858"
    demux$individual[which(demux$BEST == "SNG-NA18520")] <- "NA18520"
    demux$individual[which(demux$BEST == "SNG-NA18511")] <- "NA18511"
    demux$individual[which(demux$BEST == "SNG-NA18508")] <- "NA18508"

    # Create cell metadata object
    cell_metadata <- data.frame(demux$individual, row.names=rownames(demux))
    colnames(cell_metadata) <- "individual"
    rm(demux)

    # Assign differentiation days
    cell_metadata$diffday <- "NA"
    # col 1 for all CDs
    if (j == 1) {
      if (i == 1) {
        cell_metadata$diffday[which(cell_metadata$individual == "NA19093")] <- "Day 7"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18912")] <- "Day 3"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18520")] <- "Day 1"
      } else if (i == 2) {
        cell_metadata$diffday[which(cell_metadata$individual == "NA19093")] <- "Day 11"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18912")] <- "Day 7"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18520")] <- "Day 5"
      } else if (i == 3) {
        cell_metadata$diffday[which(cell_metadata$individual == "NA19093")] <- "Day 15"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18912")] <- "Day 11"
        # add day 0 col1 here
        cell_metadata$diffday[which(cell_metadata$individual == "NA18858")] <- "Day 0"
      }
    }

    # col 2 for all CDs
    if (j == 2) {
      if (i == 1) {
        cell_metadata$diffday[which(cell_metadata$individual == "NA18858")] <- "Day 7"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18520")] <- "Day 3"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18508")] <- "Day 1"
      } else if (i == 2) {
        cell_metadata$diffday[which(cell_metadata$individual == "NA18858")] <- "Day 11"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18520")] <- "Day 7"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18508")] <- "Day 5"
      } else if (i == 3) {
        cell_metadata$diffday[which(cell_metadata$individual == "NA18858")] <- "Day 15"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18520")] <- "Day 11"
        ## add day 0 col 2
        cell_metadata$diffday[which(cell_metadata$individual == "NA18912")] <- "Day 0"
      }
    }

    # col 3 for all CDs
    if (j == 3) {
      if (i == 1) {
        cell_metadata$diffday[which(cell_metadata$individual == "NA18912")] <- "Day 7"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18508")] <- "Day 3"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18511")] <- "Day 1"
      } else if (i == 2) {
        cell_metadata$diffday[which(cell_metadata$individual == "NA18912")] <- "Day 11"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18508")] <- "Day 7"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18511")] <- "Day 5"
      } else if (i == 3) {
        cell_metadata$diffday[which(cell_metadata$individual == "NA18912")] <- "Day 15"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18508")] <- "Day 11"
        # add day 0 col 3
        cell_metadata$diffday[which(cell_metadata$individual == "NA18520")] <- "Day 0"
      }
    }

    # col 4 for all CDs
    if (j == 4) {
      if (i == 1) {
        cell_metadata$diffday[which(cell_metadata$individual == "NA18520")] <- "Day 7"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18511")] <- "Day 3"
        cell_metadata$diffday[which(cell_metadata$individual == "NA19093")] <- "Day 1"
      } else if (i == 2) {
        cell_metadata$diffday[which(cell_metadata$individual == "NA18520")] <- "Day 11"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18511")] <- "Day 7"
        cell_metadata$diffday[which(cell_metadata$individual == "NA19093")] <- "Day 5"
      } else if (i == 3) {
        cell_metadata$diffday[which(cell_metadata$individual == "NA18520")] <- "Day 15"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18511")] <- "Day 11"
        # add day 0 col 4
        cell_metadata$diffday[which(cell_metadata$individual == "NA18508")] <- "Day 0"
      }
    }

    # col 5 for all CDs
    if (j == 5) {
      if (i == 1) {
        cell_metadata$diffday[which(cell_metadata$individual == "NA18508")] <- "Day 7"
        cell_metadata$diffday[which(cell_metadata$individual == "NA19093")] <- "Day 3"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18858")] <- "Day 1"
      } else if (i == 2) {
        cell_metadata$diffday[which(cell_metadata$individual == "NA18508")] <- "Day 11"
        cell_metadata$diffday[which(cell_metadata$individual == "NA19093")] <- "Day 7"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18858")] <- "Day 5"
      } else if (i == 3) {
        cell_metadata$diffday[which(cell_metadata$individual == "NA18508")] <- "Day 15"
        cell_metadata$diffday[which(cell_metadata$individual == "NA19093")] <- "Day 11"
        # add day 0 col 5
        cell_metadata$diffday[which(cell_metadata$individual == "NA18511")] <- "Day 0"
      }
    }

    # col 6 for all CDs
    if (j == 6) {
      if (i == 1) {
        cell_metadata$diffday[which(cell_metadata$individual == "NA18511")] <- "Day 7"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18858")] <- "Day 3"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18912")] <- "Day 1"
      } else if (i == 2) {
        cell_metadata$diffday[which(cell_metadata$individual == "NA18511")] <- "Day 11"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18858")] <- "Day 7"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18912")] <- "Day 5"
      } else if (i == 3) {
        cell_metadata$diffday[which(cell_metadata$individual == "NA18511")] <- "Day 15"
        cell_metadata$diffday[which(cell_metadata$individual == "NA18858")] <- "Day 11"
        # add day 0 col 6
        cell_metadata$diffday[which(cell_metadata$individual == "NA19093")] <- "Day 0"
      }
    }

    # Add collection day
    cell_metadata$colday <- str(i)

    # Filter out doublets and nonsense demuxlet assignments
    non_doublets <- cell_metadata$individual != "doublet"
    properly_labeled <- cell_metadata$diffday != "NA"
    cell_metadata <- cell_metadata[non_doublets & properly_labeled,]
    expression_matrix <- expression_matrix[,non_doublets & properly_labeled]
    rm(non_doublets, properly_labeled)

    # Add collection identifier to each cell
    for (k in 1:length(rownames(cell_metadata))) {
      rownames(cell_metadata)[k] <- paste0(rownames(cell_metadata)[k], "_CD", i, "col", j)
      colnames(expression_matrix)[k] <- paste0(colnames(expression_matrix)[k], "_CD", i, "col", j)
    }

    if (i == 1 & j == 1) {
      expression_master <- expression_matrix
      genes_master <- genes_matrix
      cells_master <- cell_metadata
    } else {
      save.image("./workspaces/monocle_construction_beforemerge.RData")
      master_list <- update_masters(expression_master,
        genes_master, cells_master, expression_matrix, genes_matrix, cell_metadata)
      expression_master <- master_list[[1]]
      genes_master <- master_list[[2]]
      cells_master <- master_list[[3]]
    }
    rm(expression_matrix, cell_metadata, genes_matrix)
  }
}
rm(i,j)

expression_sparse <- Matrix(as.matrix(expression_master), sparse=TRUE)

save.image("./workspaces/monocle_construction_cds_all.RData")

# create cell_data_set object
cds <- new_cell_data_set(expression_sparse,
  cell_metadata = cells_master,
  gene_metadata = genes_master)

saveRDS(cds, "./rds_objects/sc_fulldata_monocle.RDS")
