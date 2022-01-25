library(Seurat)
library(Matrix)

countMatrix <- read.table('input.txt', sep = '\t', stringsAsFactors = FALSE)
geneList <- make.unique(countMatrix[,1])[-1]
bcList <- countMatrix[1,][-1]
countMatrix <- Matrix(as.matrix(countMatrix[-1,-1]))
rownames(countMatrix) <- geneList
colnames(countMatrix) <- bcList
countMatrix <- CreateSeuratObject(countMatrix)
sce <- NormalizeData(countMatrix)
saveRDS(sce, file = "output.Rds")

#countMatrix <- ScaleData(countMatrix, features = unlist(cc.genes.updated.2019))
#countMatrix <- CellCycleScoring(countMatrix, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
#countMatrix <- data.frame(S=countMatrix$S.Score, G2M=countMatrix$G2M.Score, Phase = countMatrix$Phase)
#write.csv(countMatrix, file = 'output.csv')

