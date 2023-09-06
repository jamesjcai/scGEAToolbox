#setwd('D:\\GitHub\\scGEAToolbox\\+run\\external\\R_SPIA')
library(SPIA)
options(digits=3)
#gene_all<-as.list(read.csv("input1.txt", header = FALSE, colClasses=c('character')))
#gene_deg<-as.list(read.csv("input2.txt", header = FALSE, colClasses=c('character')))
#gene_fc<-as.matrix(read.csv("input3.txt", header = FALSE))
#names(gene_fc)<-as.vector(gene_deg);
#res=spia(de=gene_fc,all=gene_all,organism="hsa",nB=2000,plots=FALSE,beta=NULL,combine="fisher",verbose=FALSE)

 
allx<-scan("input1.txt",sep=',',what="", quiet = TRUE)
degx<-scan("input2.txt",sep=',',what="", quiet = TRUE)
dex<-scan("input3.txt",sep=',',what=numeric(), quiet = TRUE)
names(dex)<-degx
res=spia(de=dex,all=allx,organism="hsa",nB=2000,plots=FALSE,beta=NULL,combine="fisher",verbose=FALSE)
write.csv(res,'output.csv')

#data(colorectalcancer)
#library(hgu133plus2.db)
#x <- hgu133plus2ENTREZID 
#top$ENTREZ<-unlist(as.list(x[top$ID]))
#top<-top[!is.na(top$ENTREZ),]
#top<-top[!duplicated(top$ENTREZ),]
#tg1<-top[top$adj.P.Val<0.1,]
#DE_Colorectal=tg1$logFC
#names(DE_Colorectal)<-as.vector(tg1$ENTREZ)
#ALL_Colorectal=top$ENTREZ
#res=spia(de=DE_Colorectal,all=ALL_Colorectal,organism="hsa",nB=2000,plots=FALSE,beta=NULL,combine="fisher",verbose=FALSE)
#res$Name=substr(res$Name,1,10)
#res[1:20,-12]res