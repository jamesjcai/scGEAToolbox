if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
}
if (!requireNamespace("BUSseq", quietly = TRUE)){
BiocManager::install("BUSseq")
}
if (!requireNamespace("rhdf5", quietly = TRUE)){
BiocManager::install("rhdf5")
}
