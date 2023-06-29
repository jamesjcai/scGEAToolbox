library(Seurat)
library(Matrix)
library(sctransform)
library(rhdf5) 

# setwd("U:\\GitHub\\scGEAToolbox\\+run\\external\\R_SeuratSctransform")       
if (file.exists("input.mat")){
    sce.counts <- h5read(file = "input.mat", name = "/X")
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
sce <- SCTransform(sce, vars.to.regress = "percent.mt", verbose = FALSE)

#Apply sctransform normalization
#Note that this single command replaces NormalizeData, ScaleData, and FindVariableFeatures.
#https://satijalab.org/seurat/archive/v3.0/sctransform_vignette.html

# sce <- SCTransform(sce)

X<-as.matrix(sce@assays$SCT@counts)


tryCatch({    
    h5write(as.matrix(X),"output.h5","X")
},
error = function(err){
    write.table(X,file="output.txt",row.names=FALSE,col.names = FALSE) # drops the rownames
})
