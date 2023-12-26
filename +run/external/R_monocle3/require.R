if (!requireNamespace("monocle3", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
}
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))
}

if (!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools", repo="http://cran.rstudio.com/")
}

if (!requireNamespace("monocle3", quietly = TRUE)){
    devtools::install_github('cole-trapnell-lab/monocle3')
}


if (!requireNamespace("rhdf5", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
}
    BiocManager::install("rhdf5")
}
if (!(packageVersion("Matrix")=='1.6.1.1')){
    remotes::install_version("Matrix", version = "1.6-1.1")
}