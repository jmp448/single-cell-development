library(stringr)
library(Seurat)

### The following is, until the end, copied directly from Reem Elorbany's work
### found at
### https://relorbany.github.io/Exp1_prelim_analysis/scCM_round1_lowpass_combined.html

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


## My goal is to combine all 18 collections into a single Seurat object.  First, I
## need to read in the raw data for all collections.
# # Load previously saved raw data
# load("day1raw.RData")
# load("day2raw.RData")
# load("day3raw.RData")
# #OR
# load('fullrawdata.RData')

## Next, I convert Ensembl gene numbers to gene names/symbols to make it more
## readable.  To do that, first I need to get the matched Ensembl IDs and gene
## symbols through biomart.  get data from biomart
geneinfo <- readRDS("geneinfo.rds")

for (i in 1:3) {
  for (j in 1:6) {
    # read in every raw data matrix
    # raw_dat_temp <- eval(as.name(paste0("rawdata_C",i,"c",j)))
    # rm(list=(paste0("rawdata_C",i,"c",j)))
    
    raw_dat_temp <- read.table(rawdata[6*(i-1)+j],header = T, stringsAsFactors = F, row.names = 1)
    
    # remove the version numbers from those gene IDs 
    gene_id_temp <-
      str_replace(rownames(raw_dat_temp),
                  pattern = ".[0-9]+$", replacement="")
    # and subset to only those genes from the full list of genes
    mygeneinfo_temp <- geneinfo[geneinfo$ensembl_gene_id %in% gene_id_temp,]
    
    # make a temporary raw data matrix with the non-version gene IDs
    bm_rawdat_temp <- raw_dat_temp
    rownames(bm_rawdat_temp) <- gene_id_temp
    rm(gene_id_temp, raw_dat_temp)
    # subset temp raw data matrix to only the ones with info on biomart
    bm_rawdat_temp <- bm_rawdat_temp[rownames(bm_rawdat_temp) %in% mygeneinfo_temp$ensembl_gene_id,]
    
    # get all the index numbers for a gene list sorted by number
    mygeneinfo_idsort_temp <- order(mygeneinfo_temp$ensembl_gene_id,decreasing=F)
    # and sort the genes with those index numbers
    mygeneinfo_sort_temp <- mygeneinfo_temp[mygeneinfo_idsort_temp,]
    rm(mygeneinfo_temp, mygeneinfo_idsort_temp)
    # get rid of duplicate ensembl IDs
    mygeneinfo_s_uni_temp <- mygeneinfo_sort_temp[order(mygeneinfo_sort_temp$ensembl_gene_id)[!duplicated(mygeneinfo_sort_temp$ensembl_gene_id)],]
    rm(mygeneinfo_sort_temp)
    # deal with duplicate gene name/symbols (mostly ""s) by identifying the duplicates and then creating a new gene name for them that is genesymbol.ensemblID
    dupl_temp <- unique(mygeneinfo_s_uni_temp$hgnc_symbol[duplicated(mygeneinfo_s_uni_temp$hgnc_symbol)])
    for (k in dupl_temp)
    {
      rows_temp <- which(mygeneinfo_s_uni_temp$hgnc_symbol==k)
      for (m in 1:length(rows_temp))
      {
        mygeneinfo_s_uni_temp$hgnc_symbol[rows_temp[m]] <- paste(mygeneinfo_s_uni_temp$hgnc_symbol[rows_temp[m]], '.', mygeneinfo_s_uni_temp$ensembl_gene_id[rows_temp[m]],sep='')
      }
    }
    # add those new gene names/symbols as rownames in the raw data matrix
    rownames(bm_rawdat_temp) <- mygeneinfo_s_uni_temp$hgnc_symbol
    rm(mygeneinfo_s_uni_temp)
    # output the raw data matrix with the proper collection name
    # assign(paste0("bm_rawdat_C",i,"c",j),bm_rawdat_temp)
    
    saveRDS(bm_rawdat_temp, paste0("bm_rawdat_C",i,"c",j, ".RDS"))
    rm(bm_rawdat_temp)
  }
}

rm(geneinfo)