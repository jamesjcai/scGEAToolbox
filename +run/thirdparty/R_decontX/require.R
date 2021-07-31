if (!requireNamespace("celda", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
}
BiocManager::install("celda")
}
if (!requireNamespace("R.matlab", quietly = TRUE)){
    install.packages("R.matlab", repo="http://cran.rstudio.com/")
}