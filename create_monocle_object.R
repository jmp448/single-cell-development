args = commandArgs(trailingOnly = TRUE)
library(stringr)
library(Matrix)
library(monocle3)

rawdata <- c()
for (i in 1:18) {
  rawdata <- c(rawdata, args[i])
}

min_cells_per_gene = args[19]  # minimum num cells in which a gene must appear
min_genes_per_cell = args[20]  # minimum num genes for a cell to be included
cutoff_mito = args[21]  # wanna cut off cells w a certain percent mito? (bool)
mito_threshold = args[22]  # what is the threshold for mito cutoff?

# Write a function for updating the matrix
update_masters <- function (expression_master,
  genes_master, cells_master, expression_matrix, genes_matrix, cell_metadata) {

    # create new genes matrix
    unseen_genes <- subset(genes_matrix, !(rownames(genes_matrix) %in% rownames(genes_master)))
    overlap_genes <- subset(genes_matrix, rownames(genes_matrix) %in% rownames(genes_master))
    new_genes_master <- rbind(genes_master, unseen_genes)

    # create new cells matrix
    new_cells_master <- rbind(cells_master, cell_metadata)

    # create new expression matrix
    new_expression_master <- expression_master

    # add new genes into expression matrix
    num_new_genes <- length(rownames(unseen_genes))
    num_old_genes <- length(rownames(genes_master))
    new_expression_master[(num_old_genes+1):(num_old_genes+num_new_genes),] <- 0
    rownames(new_expression_master)[(num_old_genes+1):(num_old_genes+num_new_genes)]
      <- rownames(unseen_genes)

    # add new cells into expression matrix
    num_new_cells <- length(rownames(cell_metadata))
    num_old_cells <- length(rownames(cells_master))
    new_expression_master[,(num_old_cells+1):(num_old_cells+num_new_cells)] <- 0
    for (g in rownames(overlap_genes)) {
      new_expression_master[g, (num_old_cells+1):(num_old_cells+num_new_cells)]
        <- expression_matrix[g,]
    }
    for (g in rownames(new_genes)) {
      new_expression_master[g, (num_old_cells+1):(num_old_cells+num_new_cells)]
        <- expression_matrix[g,]
    }

    return(new_expression_master, new_genes_master, new_cells_master)
}

## First, read in the raw data for all collections

## Get gene ID info from biomart
geneinfo <- readRDS("./rds_objects/geneinfo.rds")

for (i in 1:1) {
  for (j in 1:2) {
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
    for (i in length(rownames(cell_metadata))) {
      rownames(cell_metadata)[i] = paste0(rownames(cell_metadata)[i], "_CD", i, "col", j)
      colnames(expression_matrix)[i] = paste0(colnames(expression_matrix)[i], "_CD", i, "col", j)
    }

    if (i == 1 and j == 1) {
      expression_master <- expression_matrix
      genes_master <- genes_matrix
      cells_master <- cell_metadata
    } else {
      expression_master, genes_master, cells_master <- update_masters(expression_master,
        genes_master, cells_master, expression_matrix, genes_matrix, cell_metadata)
    }
    rm(expression_matrix, cell_metadata, genes_matrix)
  }
}
expression_sparse <- Matrix(as.matrix(expression_master), sparse=TRUE)

save.image("monocle_construction_cds_all.RData")

# create cell_data_set object
cds <- new_cell_data_set(expression_master,
  cell_metadata = cells_master,
  gene_metadata = genes_master)

saveRDS(cds, "./rds_objects/sc_fulldata_monocle.RDS")
