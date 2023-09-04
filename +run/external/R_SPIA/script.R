library(SPIA)
library(rhdf5)


library(SPIA)
data(colorectalcancer)
options(digits=3)
head(top)


###################################################
### code chunk number 2: SPIA.Rnw:99-108
###################################################
library(hgu133plus2.db)
x <- hgu133plus2ENTREZID 
top$ENTREZ<-unlist(as.list(x[top$ID]))
top<-top[!is.na(top$ENTREZ),]
top<-top[!duplicated(top$ENTREZ),]
tg1<-top[top$adj.P.Val<0.1,]
DE_Colorectal=tg1$logFC
names(DE_Colorectal)<-as.vector(tg1$ENTREZ)
ALL_Colorectal=top$ENTREZ


res=spia(de=DE_Colorectal,all=ALL_Colorectal,organism="hsa",nB=2000,plots=FALSE,beta=NULL,combine="fisher",verbose=FALSE)
#make the output fit this screen
res$Name=substr(res$Name,1,10)
#show first 15 pathways, omit KEGG links
res[1:20,-12]


require(node2vec)
gene_edges<-read.csv("input.txt")
#emb<-node2vecR(gene_edges,p=2,q=1,num_walks=5,walk_length=5,dim=10)
emb<-node2vecR(gene_edges,directed=TRUE)
write.csv(emb,file="output.txt")


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

