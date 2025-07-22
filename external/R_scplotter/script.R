library(Seurat)
library(Matrix)
library(rhdf5)
library(scplotter)
setwd("C:\\Users\\jcai\\Documents\\GitHub\\scGEAToolbox\\+run\\external\\R_scplotter")
X <- h5read(file = "input.h5", name = "/X")
g <- h5read(file = "input.h5", name = "/g")
s <- h5read(file = "input.h5", name = "/s")
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
sce <- SCTransform(sce)
sce <- RunPCA(object = sce)
# sce <- RunUMAP(object = sce, dims = 1:30)

#celltype <- as.data.frame(celltype)
#rownames(celltype) <- colnames(x = sce)
#celltype <- setNames(celltype, colnames(x = sce))


celltype <- data.frame(strings = celltype)
rownames(celltype) <- colnames(x = sce)

s<-as.data.frame(s)
rownames(s) <- colnames(x = sce)
sce@reductions[['umap']] <- CreateDimReducObject(embeddings = s, key = 'Umap_', assay = 'RNA')


sce <- AddMetaData(
  object = sce,
  metadata = celltype,
  col.name = 'CellType'
)


a<-CellDimPlot(sce, group_by = "CellType", reduction = "umap",
            label = TRUE)

png(filename="output.png", width = 800, height = 600)
plot(a)
dev.off()

#b<-FeatureStatPlot(sce, plot_type = "dim", features = "Rbp4", reduction = "umap",
#   add_density = TRUE, density_filled = TRUE)
#png(filename="output2.png")
#plot(b)
#dev.off()

saveRDS(sce, file = "output.Rds")

