if (!requireNamespace("monocle", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
}
BiocManager::install("monocle")
}

if (!requireNamespace("rhdf5", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
}
    BiocManager::install("rhdf5")
}

#if (!requireNamespace("R.matlab", quietly = TRUE)){
#    install.packages("R.matlab", repo="http://cran.rstudio.com/")
#}
