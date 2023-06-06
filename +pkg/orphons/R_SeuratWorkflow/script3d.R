suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(R.matlab))


if (file.exists("input.mat")){
    mat<-readMat("input.mat")
    pbmc.counts<-Matrix(mat$X)
    rownames(pbmc.counts) <- make.unique(unlist(mat$genelist))
    colnames(pbmc.counts) <- paste0(colnames(pbmc.counts), 1:ncol(pbmc.counts))
    pbmc <- CreateSeuratObject(counts = pbmc.counts)
} else {
    countMatrix <- read.table('input.txt', sep = '\t', stringsAsFactors = FALSE)
    geneList <- make.unique(countMatrix[,1])[-1]
    bcList <- countMatrix[1,][-1]
    countMatrix <- Matrix(as.matrix(countMatrix[-1,-1]))
    rownames(countMatrix) <- geneList
    colnames(countMatrix) <- bcList
    pbmc <- CreateSeuratObject(countMatrix)
}

pbmc <- SCTransform(pbmc)

pbmc <- RunPCA(object = pbmc)
pbmc <- RunTSNE(object = pbmc, dims = 1:30, dim.embed = 3)
pbmc <- RunUMAP(object = pbmc, dims = 1:30, n.components = 3L)
pbmc <- FindNeighbors(object = pbmc, dims = 1:30)
pbmc <- FindClusters(object = pbmc)

s_tsne<-as.matrix(pbmc@reductions$tsne@cell.embeddings)
s_umap<-as.matrix(pbmc@reductions$umap@cell.embeddings)
c_ident<-as.matrix(as.numeric(pbmc@active.ident))
writeMat("output.mat",s_tsne=s_tsne,s_umap=s_umap,c_ident=c_ident)

