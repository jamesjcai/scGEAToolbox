#setwd("E:\\GitHub\\scGEAToolbox\\+run\\external\\R_DESeq2")
suppressMessages(library(DESeq2))
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




d <- DGEList(counts=data.counts,group=factor(metadata))
d
dim(d)
d.full <- d # keep the old one in case we mess up
head(d$counts)
head(cpm(d))
apply(d$counts, 2, sum)
keep <- rowSums(cpm(d)>100) >= 2
d <- d[keep,]
dim(d)
d$samples$lib.size <- colSums(d$counts)
d$samples
d <- calcNormFactors(d)
d
plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)

#Estimating Dispersion
d1 <- estimateCommonDisp(d, verbose=T)
names(d1)
d1 <- estimateTagwiseDisp(d1)
names(d1)
plotBCV(d1)

