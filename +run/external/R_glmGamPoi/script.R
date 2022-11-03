#setwd("E:\\GitHub\\scGEAToolbox\\+run\\external\\R_DESeq2")
suppressMessages(library(glmGamPoi))
suppressMessages(library(rhdf5))

#mat<-readMat('input.mat')
#X<-as.matrix(mat$X)
#Y<-as.matrix(mat$Y)

X <- h5read(file = "input.mat", name = "/X")
Y <- h5read(file = "input.mat", name = "/Y")

X<-as.matrix(X)
Y<-as.matrix(Y)
# rownames(X) <- rownames(Y) <- paste0('', seq_len(nrow(X)))

CT1 = t(replicate(ncol(X),1))
CT2 = t(replicate(ncol(Y),2))
metadata = t(cbind(CT1,CT2))
colnames(metadata) = ('clusters')
data.counts = cbind(X,Y)

dds <- DESeqDataSetFromMatrix(countData=data.counts,
                              colData=metadata,
                              design=~clusters)

dds <- DESeq(dds)
res <- results(dds, contrast = c(1,2))
write.csv(res,"output.csv")
