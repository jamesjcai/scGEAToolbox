library(Seurat)
library(Matrix)
library(rhdf5)


# Validate datasets in the HDF5 file
datasets <- h5ls("input.h5")
required_datasets <- c("X", "g", "celltype", "batchid")
if (!all(required_datasets %in% datasets$name)) {
    stop("One or more required datasets are missing in the HDF5 file.")
}


X <- h5read(file = "input.h5", name = "/X")
g <- h5read(file = "input.h5", name = "/g")
celltype <- h5read(file = "input.h5", name = "/celltype")
batchid <- h5read(file = "input.h5", name = "/batchid")

# Validate dimensions
if (length(g) != nrow(X)) {
    stop("Number of genes (g) does not match the number of rows in X.")
}

if (length(celltype) != ncol(X) || length(batchid) != ncol(X)) {
    stop("Lengths of celltype or batchid do not match the number of columns in X.")
}

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
sce <- SCTransform(sce)
sce <- RunPCA(object = sce)

#sce <- RunUMAP(object = sce, dims = 1:30)


#celltype <- as.data.frame(celltype)
#rownames(celltype) <- colnames(x = sce)
#celltype <- setNames(celltype, colnames(x = sce))


#sce <- AddMetaData(
#  object = sce,
#  metadata = celltype,
#  col.name = 'CellType'
#)

# Add metadata for CellType
celltype <- data.frame(CellType = celltype)
rownames(celltype) <- colnames(sce)
sce <- AddMetaData(sce, metadata = celltype$CellType, col.name = 'CellType')

sce <- AddMetaData(
  object = sce,
  metadata = as.factor(unlist(batchid)),
  col.name = 'BatchID'
)

saveRDS(sce, file = "output.Rds")

#countMatrix <- ScaleData(countMatrix, features = unlist(cc.genes.updated.2019))
#countMatrix <- CellCycleScoring(countMatrix, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
#countMatrix <- data.frame(S=countMatrix$S.Score, G2M=countMatrix$G2M.Score, Phase = countMatrix$Phase)
#write.csv(countMatrix, file = 'output.csv')
