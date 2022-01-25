suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(R.matlab))
# setwd("U:\\GitHub\\scGEAToolbox\\+run\\thirdparty\\R_SeuratWorkflow")
#pbmc.counts <- read.table('input.txt', sep = '\t', stringsAsFactors = FALSE)
#geneList <- make.unique(pbmc.counts[,1])[-1]
#bcList <- pbmc.counts[1,][-1]
#pbmc.counts <- Matrix(as.matrix(pbmc.counts[-1,-1]))
#rownames(pbmc.counts) <- geneList
#colnames(pbmc.counts) <- bcList


if (file.exists("input.mat")){
    mat<-readMat("input.mat")
    pbmc.counts<-Matrix(mat$X)
    rownames(pbmc.counts) <- make.unique(unlist(mat$genelist))
    colnames(pbmc.counts) <- paste0(colnames(pbmc.counts), 1:ncol(pbmc.counts))
    # https://satijalab.org/seurat/articles/essential_commands.html
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


#Apply sctransform normalization
#Note that this single command replaces NormalizeData, ScaleData, and FindVariableFeatures.
#https://satijalab.org/seurat/archive/v3.0/sctransform_vignette.html
pbmc <- SCTransform(pbmc)

pbmc <- RunPCA(object = pbmc)
pbmc <- RunTSNE(object = pbmc, dims = 1:30)
pbmc <- RunUMAP(object = pbmc, dims = 1:30)
pbmc <- FindNeighbors(object = pbmc, dims = 1:30)
pbmc <- FindClusters(object = pbmc)

#write.csv(pbmc@reductions$tsne@cell.embeddings, file = 'tsneoutput.csv')
#write.csv(pbmc@reductions$umap@cell.embeddings, file = 'umapoutput.csv')
#write.csv(pbmc@active.ident, file = 'activeidentoutput.csv')
s_tsne<-as.matrix(pbmc@reductions$tsne@cell.embeddings)
s_umap<-as.matrix(pbmc@reductions$umap@cell.embeddings)
c_ident<-as.matrix(as.numeric(pbmc@active.ident))
writeMat("output.mat",s_tsne=s_tsne,s_umap=s_umap,c_ident=c_ident)

