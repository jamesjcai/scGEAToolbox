library(Seurat)
library(Matrix)
library(sctransform)
library(rhdf5) 

# setwd("U:\\GitHub\\scGEAToolbox\\+run\\external\\R_SeuratSctransform")       
if (file.exists("input.mat")){
    sce.counts <- h5read(file = "input.mat", name = "/X")
    sce.counts <- as(sce.counts, "dgCMatrix")
    genelist<-read.table('g.txt', sep = '\t', stringsAsFactors = FALSE)
    rownames(sce.counts) <- make.unique(unlist(genelist))
    colnames(sce.counts) <- paste0(colnames(sce.counts), 1:ncol(sce.counts))
    sce <- CreateSeuratObject(counts = sce.counts)
} else {
    countMatrix <- read.table('input.txt', sep = '\t', stringsAsFactors = FALSE)
    geneList <- make.unique(countMatrix[,1])[-1]
    bcList <- countMatrix[1,][-1]
    countMatrix <- Matrix(as.matrix(countMatrix[-1,-1]))
    rownames(countMatrix) <- geneList
    colnames(countMatrix) <- bcList
    sce <- CreateSeuratObject(countMatrix)
}

# https://satijalab.org/seurat/articles/sctransform_vignette.html
# store mitochondrial percentage in object meta data
sce <- PercentageFeatureSet(sce, pattern = "^MT-", col.name = "percent.mt")
sce <- SCTransform(sce, vars.to.regress = "percent.mt", verbose = FALSE, vst.flavor = "v2")

#Apply sctransform normalization
#Note that this single command replaces NormalizeData, ScaleData, and FindVariableFeatures.
#https://satijalab.org/seurat/archive/v3.0/sctransform_vignette.html

# sce <- SCTransform(sce)

#counts<-as.matrix(sce@assays$SCT@counts)
#data<-as.matrix(sce@assays$SCT@data)
#scale_data<-as.matrix(sce@assays$SCT@scale.data)
#Access normalized data: 
data <- GetAssayData(sce, assay = "SCT", slot = "data")
#Access Pearson residuals: 
scale_data <- GetAssayData(sce, assay = "SCT", slot = "scale.data")

#Access counts (may be empty): 
#counts <- GetAssayData(sce, assay = "SCT", slot = "counts")

tryCatch({
    if (file.exists("output.h5")) {
      file.remove("output.h5")
    }
    h5write(as.matrix(data),"output.h5","data")
    h5write(as.matrix(scale_data),"output.h5","scale_data")
},
error = function(err){
    write.table(as.matrix(data),file="output_data.txt",row.names=FALSE,col.names = FALSE) # drops the rownames
    write.table(as.matrix(scale_data),file="output_scale_data.txt",row.names=FALSE,col.names = FALSE) # drops the rownames
})
