if (!requireNamespace("sctransform", quietly = TRUE)){
    install.packages("sctransform", repo="http://cran.rstudio.com/")
}
if (!requireNamespace("Matrix", quietly = TRUE)){
    install.packages("Matrix", repo="http://cran.rstudio.com/")
}
if (!requireNamespace("Seurat", quietly = TRUE)){
    install.packages("Seurat", repo="http://cran.rstudio.com/")
}
if (!requireNamespace("rhdf5", quietly = TRUE)){
    install.packages("rhdf5", repo="http://cran.rstudio.com/")
}

if (!requireNamespace("rhdf5", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
}
BiocManager::install("rhdf5")
}

if (!requireNamespace("glmGamPoi", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
}
BiocManager::install('glmGamPoi')
}