if (!requireNamespace("iDESC", quietly = TRUE)){
if (!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools", repo="http://cran.rstudio.com/")
}
devtools::install_github("yl883/iDESC")
}


if (!requireNamespace("rhdf5", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
}
BiocManager::install("rhdf5")
}

