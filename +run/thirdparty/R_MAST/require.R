#setwd(getSrcDirectory()[1])
#setwd("e:\\GitHub\\scGEAToolbox\\thirdparty\\R_MAST")

if (!requireNamespace("MAST", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}
    BiocManager::install("MAST")
}

if (!requireNamespace("Seurat", quietly = TRUE)){
    install.packages("Seurat")
}
