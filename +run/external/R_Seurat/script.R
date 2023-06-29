suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(rhdf5))


cMatrix <- h5read(file = "input.mat", name = "/X")
ndim <- h5read(file = "input.mat", name = "/ndim")

g <- read.table('g.txt')
rownames(cMatrix) <- t(g)
colnames(cMatrix) <- paste0("C", seq_len(ncol(cMatrix)))
pbmc <- CreateSeuratObject(counts = cMatrix)

#Apply sctransform normalization
#Note that this single command replaces NormalizeData, ScaleData, and FindVariableFeatures.
#https://satijalab.org/seurat/archive/v3.0/sctransform_vignette.html
pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(object = pbmc)
pbmc <- RunTSNE(object = pbmc, dims = 1:30, dim.embed = ndim[1])

if (ndim[1] == 3) {
    pbmc <- RunUMAP(object = pbmc, dims = 1:30, dn.components = 3L)
} else {
    pbmc <- RunUMAP(object = pbmc, dims = 1:30)
}


pbmc <- FindNeighbors(object = pbmc, dims = 1:30)
pbmc <- FindClusters(object = pbmc)
s_tsne<-as.matrix(pbmc@reductions$tsne@cell.embeddings)
s_umap<-as.matrix(pbmc@reductions$umap@cell.embeddings)
c_ident<-as.matrix(as.numeric(pbmc@active.ident))

h5write(s_tsne, "output.h5", "s_tsne")
h5write(s_umap, "output.h5", "s_umap")
h5write(c_ident, "output.h5", "c_ident")