#setwd(getSrcDirectory()[1])
#setwd("e:\\GitHub\\scGEAToolbox\\thirdparty\\R_MAST")

if (!requireNamespace("MAST", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
}
    BiocManager::install("MAST")
}

if (!requireNamespace("rhdf5", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
}
    BiocManager::install("rhdf5")
}


if (!requireNamespace("Seurat", quietly = TRUE)){
    install.packages("Seurat", repo="http://cran.rstudio.com/")
}

#if (!requireNamespace("R.matlab", quietly = TRUE)){
#    install.packages("R.matlab", repo="http://cran.rstudio.com/")
#}