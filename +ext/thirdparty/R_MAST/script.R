#setwd("u:\\GitHub\\scGEAToolbox\\+run\\thirdparty\\R_MAST")
suppressMessages(library(MAST))
suppressMessages(library(Seurat))
suppressMessages(library(rhdf5))

#mat<-readMat('input.mat')
#X<-as.matrix(mat$X)
#Y<-as.matrix(mat$Y)

X <- h5read(file = "input.mat", name = "/X")
Y <- h5read(file = "input.mat", name = "/Y")

X<-as.matrix(X)
Y<-as.matrix(Y)

#X <- as.matrix(read.table("input1.txt", sep=","))
#Y <- as.matrix(read.table("input2.txt", sep=","))
rownames(X) <- rownames(Y) <- paste0('', seq_len(nrow(X)))
colnames(X) <- paste0(colnames(X), 1:ncol(X))
colnames(Y) <- paste0(colnames(Y), 1:ncol(Y))

X <- CreateSeuratObject(X, project = 'X')
Y <- CreateSeuratObject(Y, project = 'Y')

ALL <- merge(X,Y)
ALL <- NormalizeData(ALL)
ALL <- ScaleData(ALL)
DE <- FindMarkers(ALL, ident.1 = 'X', ident.2 = 'Y', logfc.threshold = 0, test.use = 'MAST')
write.csv(DE, 'output.csv')
