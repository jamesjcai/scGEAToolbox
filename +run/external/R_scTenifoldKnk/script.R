#setwd('U:\\GitHub\\scGEAToolbox\\+run\\thirdparty\\R_scTenifoldKnk')
#setwd('C:\\Users\\jcai.AUTH\\Documents\\GitHub\\scGEAToolbox\\+run\\thirdparty\\R_scTenifoldKnk')
library(rhdf5)
library(scTenifoldKnk)

X <- h5read(file = "input.h5", name = "/X")
g <- h5read(file = "input.h5", name = "/g")
targetg <- h5read(file = "input.h5", name = "/targetg")

rownames(X) <- g
#rownames(X) <- read.csv('genelist.txt', header=FALSE)[,1]
colnames(X) <- paste0("C", seq_len(ncol(X)))

KO_res <- scTenifoldKnk(X, gKO=targetg, qc_minLSize = 0, nc_nNet = 1)
# save(KO_res, file = 'output.RData')
write.csv(KO_res$diffRegulation,'output.txt')