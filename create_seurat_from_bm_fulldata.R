library(stringr)
library(Seurat)

# Now let's create the individual Seurat objects for each collection. I do this
# first so I can assign individual and diffday labels more easily, before merging
# and subsetting the object.  create individual seurat objects

for (i in 1:3) {
  for (j in 1:6) {
    
    bmrawdat_temp <- readRDS(paste0("bm_rawdat_C", i, "c", j, ".RDS"))
    
    
    
  }
}

rm(i, j)

# these are Seurat objects for each collection with gene symbols!  filtered by
# only genes in at least 3 cells and only cells with at least 200 genes

# this filter destroyed a few samples (in terms of cell numbers): CD2col1,
# CD2col2, CD2col6 especially maybe i should rethink how strict it is?

# ``` The next step is assigning labels for 'individual' and 'diffday' to each
# cell barcode. I do this using the results from demuxlet, which was run
# previously on the data from each collection.  ```{r demux for the combined one}
# this one runs p fast
for (i in 1:3) {
  for (j in 1:6) {
    
    
  }
}

rm(i, j, m, tmp_dday, tmp_ind, SObject, demux_temp)
rm(list = ls(pattern = "SObj"))

# if we want to explore the data, can do: table(CD1col1_lbld$individual)
# table(CD1col1_lbld$diffday)

# Now that I've added metadata tags for 'individual' and 'diffday' for each
# collection, I can combine data all collections into a single Seurat object. I
# will subset to remove all cells that were labeled doublets, or cells that had
# the wrong individual label given the collection (which I check separately).

# r merge all Seurat objects and subset out useless cells

# merge AFTER giving them individual and diffday labels
all_cols_S <- merge(CD1col1_lbld, y = c(CD1col2_lbld, CD1col3_lbld, CD1col4_lbld, 
  CD1col5_lbld, CD1col6_lbld, CD2col1_lbld, CD2col2_lbld, CD2col3_lbld, CD2col4_lbld, 
  CD2col5_lbld, CD2col6_lbld, CD3col1_lbld, CD3col2_lbld, CD3col3_lbld, CD3col4_lbld, 
  CD3col5_lbld, CD3col6_lbld), add.cell.ids = c("CD1col1", "CD1col2", "CD1col3", 
  "CD1col4", "CD1col5", "CD1col6", "CD2col1", "CD2col2", "CD2col3", "CD2col4", 
  "CD2col5", "CD2col6", "CD3col1", "CD3col2", "CD3col3", "CD3col4", "CD3col5", 
  "CD3col6"), project = "scCM_round1_fulldata")

# this object has 26381 features across 195218 samples

table(all_cols_S$orig.ident)
# a range of cell nos per col from 3045 in C2c2 to 19725 in C1c1

# now i want to subset so only the cells with an individual and a diffday are
# being used (i.e. no mislabeled individuals on each collection and no doublets)
Idents(all_cols_S) <- "diffday"
# make diffday the relevant identity marker
table(Idents(all_cols_S))
# a range of cell nos per diffday from 10026 in D15 to 33179 in D1 (lots of NAs:
# 44468)

all_cols_noNA_S = subset(all_cols_S, idents = c("Day 1", "Day 3", "Day 5", "Day 7", 
  "Day 11", "Day 15", "Day 0"))

# now this object has 26381 features across 150750 samples there should be no
# doublets or mis-IDed cells here

# finally, i'll remove the individual Seurat objects for each col
rm(list = ls(pattern = "lbld"))
# and the raw data which i shouldn't need anymore
rm(list = ls(pattern = "rawdata"))

# and save this RObj so I can get to it again without running all the previous
# code saveRDS(all_cols_noNA_S,file='Sobj_allhp_noNAdays.rds') and to load
# all_cols_noNA_S <- readRDS('Sobj_allhp_noNAdays.rds')


# I will make a 'sample' label in the metadata tag, which has the combined
# individual and diffday information in the form: NAInd.Day# make a sample label

# let's make a column that's a combo of individual and diffday
all_cols_noNA_S$sample <- "IND.DAY"

all_cols_noNA_S$sample <- paste0(all_cols_noNA_S$individual, ".", all_cols_noNA_S$diffday)

# i will also add metadata info for percent mitochondrial genes
all_cols_noNA_S[["percent.mito"]] <- PercentageFeatureSet(all_cols_noNA_S, pattern = "^MT-")

# first, let's order the levels in a more intuitive way for diffday
mydaylevels <- c("Day 0", "Day 1", "Day 3", "Day 5", "Day 7", "Day 11", "Day 15")
sc@meta.data$diffday <- factor(x = sc@meta.data$diffday, levels = mydaylevels)
# for individual
myindlevels <- c("NA19093", "NA18858", "NA18912", "NA18520", "NA18508", "NA18511")
sc@meta.data$individual <- factor(x = sc@meta.data$individual, levels = myindlevels)
# for sample
mysamplelevels <- as.character(sc$sample)
mysamplelevels <- mixedsort(mysamplelevels, decreasing = F)
mysamplelevels <- unique(mysamplelevels)
sc@meta.data$sample <- factor(sc@meta.data$sample, levels = mysamplelevels)

# Assign colday to the cells
sc$colday <- "colday"
sc$colday <- substr(sc$orig.ident, 3, 3)
table(sc$colday)

# get rid of high mito
sc <- subset(sc, subset = percent.mito < 30)

saveRDS(sc, "sc_fulldata.RDS")
