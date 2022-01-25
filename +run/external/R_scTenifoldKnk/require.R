if (!requireNamespace("scTenifoldKnk", quietly = TRUE)){
    install.packages("scTenifoldKnk", repo="http://cran.rstudio.com/")
}

if (!requireNamespace("rhdf5", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
}
    BiocManager::install("rhdf5")
}
