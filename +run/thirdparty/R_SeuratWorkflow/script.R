# setwd("U:\\GitHub\\scGEAToolbox\\+run\\thirdparty\\R_SeuratWorkflow")
library(Seurat)
library(Matrix)

pbmc.counts <- read.table('input.txt', sep = '\t', stringsAsFactors = FALSE)
geneList <- make.unique(pbmc.counts[,1])[-1]
bcList <- pbmc.counts[1,][-1]
pbmc.counts <- Matrix(as.matrix(pbmc.counts[-1,-1]))
rownames(pbmc.counts) <- geneList
colnames(pbmc.counts) <- bcList

# https://satijalab.org/seurat/articles/essential_commands.html
pbmc <- CreateSeuratObject(counts = pbmc.counts)
pbmc <- NormalizeData(object = pbmc)
pbmc <- FindVariableFeatures(object = pbmc)
#all.genes <- rownames(pbmc)
#pbmc <- ScaleData(object = pbmc, features = all.genes)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc)
pbmc <- FindNeighbors(object = pbmc)
pbmc <- FindClusters(object = pbmc)
pbmc <- RunTSNE(object = pbmc)
pbmc <- RunUMAP(object = pbmc, dims = 1:10)
write.csv(pbmc@reductions$tsne@cell.embeddings, file = 'tsneoutput.csv')
write.csv(pbmc@reductions$umap@cell.embeddings, file = 'umapoutput.csv')
write.csv(pbmc@active.ident, file = 'activeidentoutput.csv')


#pbmc <- CellCycleScoring(pbmc, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
#pbmc <- data.frame(S=pbmc$S.Score, G2M=pbmc$G2M.Score, Phase = pbmc$Phase)
# write.csv(pbmc, file = 'cellcycleoutput.csv')

# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
# pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(pbmc, file = "output.rds")


# DimPlot(object = pbmc, reduction = "tsne")
# pbmc <- CreateSeuratObject(pbmc)
# pbmc <- NormalizeData(pbmc)
# pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(pbmc)
# pbmc <- ScaleData(pbmc, features = all.genes)
# pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:20)
# pbmc <- RunUMAP(pbmc, dims = 1:10, n.components = 3L)
# pbmc <- RunUMAP(pbmc, dims = 1:10)
# write.csv(pbmc@reductions$umap@cell.embeddings, file = 'output.csv')
#pbmc <- CellCycleScoring(pbmc, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
#pbmc <- data.frame(S=pbmc$S.Score, G2M=pbmc$G2M.Score, Phase = pbmc$Phase)
#write.csv(pbmc, file = 'output.csv')
