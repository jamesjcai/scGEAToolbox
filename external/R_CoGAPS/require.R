if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repo = "http://cran.rstudio.com/")
}

if (!requireNamespace("CoGAPS", quietly = TRUE)) {
    BiocManager::install("CoGAPS")
}

if (!requireNamespace("rhdf5", quietly = TRUE)) {
    BiocManager::install("rhdf5")
}
