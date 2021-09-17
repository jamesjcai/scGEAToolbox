# setwd("U:\\GitHub\\scGEAToolbox\\+run\\thirdparty\\R_SuperCell")
suppressMessages(library(R.matlab))
suppressMessages(library(SuperCell))
mat<-readMat("input.mat")
X<-mat$X
rownames(X) <- paste0(rownames(X), 1:nrow(X))
colnames(X) <- paste0(colnames(X), 1:ncol(X))
gammavalue <- mat$gammavalue # graining level    default = 5
kvalue<-mat$kvalue 

SC <- SCimplify(X,  # gene expression matrix 
                k.knn = kvalue, # number of nearest neighbors to build kNN network
                gamma = gammavalue, # graining level
                n.var.genes = 1000) # number of the top variable genes to use for dimentionality reduction 

X2 <- supercell_GE(X, SC$membership)
writeMat("output.mat",X2=X2)
