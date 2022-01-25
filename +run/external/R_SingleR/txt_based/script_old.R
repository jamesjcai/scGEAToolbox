# https://nbisweden.github.io/excelerate-scRNAseq/session-celltypeid/celltypeid.html

#setwd(getSrcDirectory()[1])
#setwd("e:\\GitHub\\scGEAToolbox\\thirdparty\\R_SingleR")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("SingleR")
#install.packages("Seurat")

library(SingleR)
hpca.se <- HumanPrimaryCellAtlasData()

library(Seurat)
X <- as.matrix(read.table("input.txt", sep="\t", row.names = 1, header=TRUE))
Xx <- CreateSeuratObject(X, project = 'Xx')
Xx <- NormalizeData(Xx, normalization.method = "LogNormalize", scale.factor = 10000)

# genelist<-rownames(X)
# genelist <- read.csv('input.txt', row.names = 1)
# rownames(X) <- paste0('', seq_len(nrow(X)))
pred <- SingleR(test = Xx, ref = hpca.se, labels = hpca.se$label.main)
write.csv(pred$labels, 'output.csv')
