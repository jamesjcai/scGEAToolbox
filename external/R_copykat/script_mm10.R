library(copykat)
#library(Seurat)
#library(Matrix)
library(rhdf5)

X <- h5read(file = "input.h5", name = "/X")
g <- h5read(file = "input.h5", name = "/g")

countMatrix <- as.matrix(X)
rownames(countMatrix) <- g
colnames(countMatrix) <- paste0("C", seq_len(ncol(countMatrix)))

exp.rawdata <- countMatrix
copykat.test <- copykat(rawmat=exp.rawdata, id.type="S", 
ngene.chr=5, win.size=25, KS.cut=0.1, 
sam.name="test", distance="euclidean", 
norm.cell.names="",output.seg=FALSE,
plot.genes=FALSE, genome="mm10", n.cores=1)

#pred.test <- data.frame(copykat.test$prediction)
#pred1.test <- pred.test[which(pred.test$copykat.pred %in% c("aneuploid","diploid")),]  ##keep defined cells
#CNA.test <- data.frame(copykat.test$CNAmat)

#write.csv(pred.test,'output1.txt')
#write.csv(pred1.test,'output2.txt')
#write.csv(CNA.test,'output3.txt')

