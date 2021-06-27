library(Seurat)
library(Matrix)

countMatrix <- read.table('input.txt', sep = '\t', stringsAsFactors = FALSE)
geneList <- make.unique(countMatrix[,1])[-1]
bcList <- countMatrix[1,][-1]
countMatrix <- Matrix(as.matrix(countMatrix[-1,-1]))
rownames(countMatrix) <- geneList
colnames(countMatrix) <- bcList
countMatrix <- CreateSeuratObject(countMatrix)
countMatrix <- NormalizeData(countMatrix)
countMatrix <- FindVariableFeatures(countMatrix, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(countMatrix)
countMatrix <- ScaleData(countMatrix, features = all.genes)
countMatrix <- RunPCA(countMatrix, features = VariableFeatures(object = countMatrix))
#countMatrix <- RunUMAP(countMatrix, reduction = "pca", dims = 1:20)
# countMatrix <- RunUMAP(countMatrix, dims = 1:10, n.components = 3L)
countMatrix <- RunUMAP(countMatrix, dims = 1:10)
write.csv(countMatrix@reductions$umap@cell.embeddings, file = 'output.csv')
#countMatrix <- CellCycleScoring(countMatrix, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
#countMatrix <- data.frame(S=countMatrix$S.Score, G2M=countMatrix$G2M.Score, Phase = countMatrix$Phase)
#write.csv(countMatrix, file = 'output.csv')
