#setwd("U:\\GitHub\\scGEAToolbox\\+run\\thirdparty\\R_decontX")
library(celda)
library(R.matlab)

#sc <- read.table('input.csv', sep = ',', stringsAsFactors = FALSE)
inputd <-readMat('input.mat')
sc <- inputd$X
res <- decontX(x = sc)
X=res$decontXcounts
contamination=res$contamination
writeMat("output.mat", X = X, contamination=contamination)
# write.table(as.matrix(X),file="output.csv", sep=",",col.names=FALSE,row.names = FALSE)
