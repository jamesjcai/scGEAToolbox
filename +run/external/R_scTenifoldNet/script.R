#setwd('U:\\GitHub\\scGEAToolbox\\+run\\thirdparty\\R_scTenifoldNet')
#setwd('C:\\Users\\jcai.AUTH\\Documents\\GitHub\\scGEAToolbox\\+run\\thirdparty\\R_scTenifoldNet')
#library(rhdf5)
library(scTenifoldNet)

A0<-readRDS("input1.rds")
A1<-readRDS("input2.rds")
 
X0=A0@assays$RNA@counts
X1=A1@assays$RNA@counts

Output <- scTenifoldNet(X = X0, Y = X1, qc_minLibSize = 0)
write.csv(Output$diffRegulation,'output.txt')