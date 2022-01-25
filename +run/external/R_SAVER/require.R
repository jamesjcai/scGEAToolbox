if (!requireNamespace("SAVER", quietly = TRUE)){
    install.packages("SAVER", repo="http://cran.rstudio.com/")
}
if (!requireNamespace("R.matlab", quietly = TRUE)){
    install.packages("R.matlab", repo="http://cran.rstudio.com/")
}

if (!requireNamespace("rhdf5", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
}
BiocManager::install("rhdf5")
}
