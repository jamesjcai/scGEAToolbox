#setwd(getSrcDirectory()[1])
#setwd("e:\\GitHub\\scGEAToolbox\\thirdparty\\R_MAST")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

#BiocManager::install("MAST")
#install.packages("Seurat")

if (!requireNamespace("MAST", quietly = TRUE))
    BiocManager::install("MAST")

if (!requireNamespace("Seurat", quietly = TRUE))
    install.packages("Seurat")
    
#library(MAST)
#library(Seurat)