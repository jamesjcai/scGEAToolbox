library(Seurat)
library(Matrix)

A <- read.table('input.txt', sep = '\t', stringsAsFactors = FALSE)
geneList <- make.unique(A[,1])[-1]
bcList <- A[1,][-1]
A <- Matrix(as.matrix(A[-1,-1]))
rownames(A) <- geneList
colnames(A) <- bcList
A <- CreateSeuratObject(A)

A <- NormalizeData(A)
A <- FindVariableFeatures(A,nfeatures = 2000)
A <- ScaleData(A)
A <- RunPCA(A)
A <- RunUMAP(A,dims = 1:50)
write.csv(A@reductions$umap@cell.embeddings, file = 'output.csv')



A <- NormalizeData(A)
A <- FindVariableFeatures(A, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(A)
A <- ScaleData(A, features = all.genes)
A <- RunPCA(A, features = VariableFeatures(object = A))
#A <- RunUMAP(A, reduction = "pca", dims = 1:20)
# A <- RunUMAP(A, dims = 1:10, n.components = 3L)
A <- RunUMAP(A, dims = 1:10)
write.csv(A@reductions$umap@cell.embeddings, file = 'output.csv')
#A <- CellCycleScoring(A, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
#A <- data.frame(S=A$S.Score, G2M=A$G2M.Score, Phase = A$Phase)
#write.csv(A, file = 'output.csv')



