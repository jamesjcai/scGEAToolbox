# setwd("U:\\GitHub\\scGEAToolbox\\+run\\thirdparty\\R_SAVER")
suppressMessages(library(R.matlab))
suppressMessages(library(SAVER))
mat<-readMat("input.mat")
X<-mat$X

X2 <- saver(X, ncores = 12, estimates.only = TRUE)

writeMat("output.mat",X2=X2)
