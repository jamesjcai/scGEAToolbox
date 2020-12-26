#setwd(getSrcDirectory()[1])
#setwd("e:\\GitHub\\scGEAToolbox\\thirdparty\\R_MAST")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("MAST")
#install.packages("Seurat")

library(MAST)
library(Seurat)
X <- as.matrix(read.table("input1.txt", sep=","))
Y <- as.matrix(read.table("input2.txt", sep=","))

rownames(X) <- rownames(Y) <- paste0('', seq_len(nrow(X)))

X <- CreateSeuratObject(X, project = 'X')
Y <- CreateSeuratObject(Y, project = 'Y')

ALL <- merge(X,Y)
ALL <- NormalizeData(ALL)
ALL <- ScaleData(ALL)
DE <- FindMarkers(ALL, ident.1 = 'X', ident.2 = 'Y', logfc.threshold = 0, test.use = 'MAST')
write.csv(DE, 'output.csv')