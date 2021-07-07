#setwd("U:\\GitHub\\scGEAToolbox\\+run\\thirdparty\\R_decontX")
library(celda)
library(R.matlab)

#sc <- read.table('input.csv', sep = ',', stringsAsFactors = FALSE)
inputd <-readMat('input.mat')
sc <- inputd$X
pbmc4k <- decontX(x = sc)

dsc<-data.matrix(pbmc4k$decontXcounts)
# write.table(dsc,file="output.csv", sep=",",col.names=FALSE,row.names = FALSE)
writeMat("output.mat", X = dsc)