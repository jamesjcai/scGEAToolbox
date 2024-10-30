library(Seurat)
library(Matrix)
library(rhdf5)

X <- h5read(file = "input.h5", name = "/X")
g <- h5read(file = "input.h5", name = "/g")
celltype <- h5read(file = "input.h5", name = "/celltype")


# all_datasets <- h5ls(file_path)
#if (dataset_name %in% all_datasets$name) {
#  data <- h5read(file_path, dataset_name)
#} else {  
#  cat("Dataset not found in the HDF5 file.\n")
#}

#countMatrix <- read.table('input.txt', sep = '\t', stringsAsFactors = FALSE)
#geneList <- make.unique(countMatrix[,1])[-1]
#bcList <- countMatrix[1,][-1]

countMatrix <- Matrix(as.matrix(X))
rownames(countMatrix) <- g
colnames(countMatrix) <- paste0("C", seq_len(ncol(countMatrix)))

countMatrix <- CreateSeuratObject(countMatrix)
sce <- NormalizeData(countMatrix)
#sce <- SCTransform(sce)
#sce <- RunPCA(object = sce)
#sce <- RunUMAP(object = sce, dims = 1:30)


#celltype <- as.data.frame(celltype)
#rownames(celltype) <- colnames(x = sce)
#celltype <- setNames(celltype, colnames(x = sce))


celltype <- data.frame(strings = celltype)
rownames(celltype) <- colnames(x = sce)

sce <- AddMetaData(
  object = sce,
  metadata = celltype,
  col.name = 'CellType'
)
saveRDS(sce, file = "output.Rds")

#countMatrix <- ScaleData(countMatrix, features = unlist(cc.genes.updated.2019))
#countMatrix <- CellCycleScoring(countMatrix, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
#countMatrix <- data.frame(S=countMatrix$S.Score, G2M=countMatrix$G2M.Score, Phase = countMatrix$Phase)
#write.csv(countMatrix, file = 'output.csv')
