suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(R.matlab))
suppressMessages(library(sctransform))
library(rhdf5)


# setwd("U:\\GitHub\\scGEAToolbox\\+run\\external\\R_SeuratSctransform")
setwd("C:\\Users\\jcai\\Documents\\GitHub\\scGEAToolbox\\+run\\external\\R_SeuratSctransform")       
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


# https://satijalab.org/seurat/articles/sctransform_vignette.html
# store mitochondrial percentage in object meta data
# pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
# pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)


#Apply sctransform normalization
#Note that this single command replaces NormalizeData, ScaleData, and FindVariableFeatures.
#https://satijalab.org/seurat/archive/v3.0/sctransform_vignette.html

pbmc <- SCTransform(pbmc)

X2<-as.matrix(pbmc@assays$SCT@counts)

library(rhdf5)

tryCatch({
# writeMat("output.mat",X2=X2)
h5write(as.matrix(X2),"output.h5","X2")
},
error = function(err){
write.table(X2,file="output.txt",row.names=FALSE,col.names = FALSE) # drops the rownames
})

