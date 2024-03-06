if (!requireNamespace("Matrix", quietly = TRUE)){
    install.packages("Matrix", repo="http://cran.rstudio.com/")
}
if (!requireNamespace("SCORPION", quietly = TRUE)){
    install.packages("SCORPION", repo="http://cran.rstudio.com/")
}
if (!requireNamespace("rhdf5", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
}
    BiocManager::install("rhdf5")
}