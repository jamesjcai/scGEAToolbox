# setwd("U:\\GitHub\\scGEAToolbox\\+run\\thirdparty\\R_decontX")
library(celda)
library(rhdf5)
# library(R.matlab)

#sc <- read.table('input.csv', sep = ',', stringsAsFactors = FALSE)
#inputd <-readMat('input6.mat')
#sc2 <- inputd$X

sc <- h5read(file = "input.mat", name = "/X")
res <- decontX(x = sc)
X=res$decontXcounts
contamination=res$contamination
h5createFile("output.h5")
h5write(as.matrix(X), "output.h5","/X")
h5write(contamination, "output.h5","/contamination")

#writeMat("output.mat", X = X, contamination=contamination)
#write.table(as.matrix(X),file="output.csv", sep=",",col.names=FALSE,row.names = FALSE)
