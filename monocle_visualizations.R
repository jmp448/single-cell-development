library(monocle3)

sc <- readRDS("./rds_objects/sc_fulldata_monocle.RDS")

sc <- preprocess_cds(sc, norm_method="log", scaling=T, verbose=T)