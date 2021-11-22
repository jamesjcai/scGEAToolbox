if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
}

if (!requireNamespace("celda", quietly = TRUE)){
BiocManager::install("celda")
}

if (!requireNamespace("rhdf5", quietly = TRUE)){
BiocManager::install("celda")
}
