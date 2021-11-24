#computeUmap <- function(X, outputFile){
#  cMatrix <- read.csv(X, header = FALSE)
#  require(uwot)
#  s_umap <- umap(cMatrix, n_components = 3, n_neighbors = 50, learning_rate = 0.5, init = "random")
#  write.csv(s_umap, file = outputFile)
#}
#computeUmap(X = "input.csv", outputFile = "output.csv")
# setwd("E:\\GitHub\\scGEAToolbox\\thirdparty\\R_BUSseq")

library(BUSseq)
library(rhdf5)

X0 <- as.matrix(read.csv("input0.csv", header = FALSE))
X1 <- as.matrix(read.csv("input1.csv", header = FALSE))
CountData<-list(X0,X1)

BUSseqfits_res <- BUSseq_MCMC(ObservedData = CountData, n.celltypes = 4, n.iterations = 500, working_dir = ".", showIteration = TRUE, seed = 123)
corrected_countdata <- corrected_read_counts(BUSseqfits_res)
write.csv(corrected_countdata[[1]], file = "output0.csv")
write.csv(corrected_countdata[[2]], file = "output1.csv")

