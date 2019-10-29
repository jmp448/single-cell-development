if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager", repos = "http://cran.us.r-project.org")
BiocManager::install()

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor'))

install.packages("devtools", repos = "http://cran.us.r-project.org")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
