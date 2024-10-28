library(copykat)
#library(Seurat)
library(Matrix)
library(rhdf5)

X <- h5read(file = "input.h5", name = "/X")
g <- h5read(file = "input.h5", name = "/g")

countMatrix <- Matrix(as.matrix(X))
rownames(countMatrix) <- g
colnames(countMatrix) <- paste0("C", seq_len(ncol(countMatrix)))

exp.rawdata <- countMatrix
copykat.test <- copykat(rawmat=exp.rawdata, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, sam.name="test", distance="euclidean", norm.cell.names="",output.seg="FLASE", plot.genes="TRUE", genome="hg20",n.cores=1)
