library(Seurat)
library(Matrix)
library(scTenifoldNet)

pbmc.counts <- read.table('input.txt', sep = '\t', stringsAsFactors = FALSE)
pbmc.counts <- scQC(pbmc.counts)
geneList <- make.unique(pbmc.counts[,1])[-1]
bcList <- pbmc.counts[1,][-1]
pbmc.counts <- Matrix(as.matrix(pbmc.counts[-1,-1]))
rownames(pbmc.counts) <- geneList
colnames(pbmc.counts) <- bcList

# https://satijalab.org/seurat/articles/essential_commands.html
pbmc <- CreateSeuratObject(counts = pbmc.counts)

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc,nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)
pbmc <- RunUMAP(pbmc,dims = 1:50)


#Apply sctransform normalization
#Note that this single command replaces NormalizeData, ScaleData, and FindVariableFeatures.
#https://satijalab.org/seurat/archive/v3.0/sctransform_vignette.html

#pbmc <- SCTransform(pbmc)
#pbmc <- RunPCA(object = pbmc)
pbmc <- RunTSNE(object = pbmc, dims = 1:50)
#pbmc <- RunUMAP(object = pbmc, dims = 1:30)
pbmc <- FindNeighbors(object = pbmc, dims = 1:50)
pbmc <- FindClusters(object = pbmc)

write.csv(pbmc@reductions$tsne@cell.embeddings, file = 'tsneoutput.csv')
write.csv(pbmc@reductions$umap@cell.embeddings, file = 'umapoutput.csv')
write.csv(pbmc@active.ident, file = 'activeidentoutput.csv')
