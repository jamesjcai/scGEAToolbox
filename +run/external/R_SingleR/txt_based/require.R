if (!requireNamespace("SingleR", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
BiocManager::install("SingleR")
}

if (!requireNamespace("scater", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
BiocManager::install("scater")
}

if (!requireNamespace("celldex", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
BiocManager::install("celldex")
}
