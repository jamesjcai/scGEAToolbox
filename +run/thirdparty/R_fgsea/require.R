# library(fgsea)

if (!requireNamespace("fgsea", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}
    BiocManager::install("fgsea")
}
