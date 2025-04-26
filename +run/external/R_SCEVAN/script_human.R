library(SCEVAN)
#library(Seurat)
#library(Matrix)
library(rhdf5)

X <- h5read(file = "input.h5", name = "/X")
g <- h5read(file = "input.h5", name = "/g")

countMatrix <- as.matrix(X)
rownames(countMatrix) <- g
colnames(countMatrix) <- paste0("C", seq_len(ncol(countMatrix)))

results <- pipelineCNA(countMatrix, SUBCLONES = FALSE, plotTree = FALSE, organism = "human")
write.csv(results,'output.csv')

#pipelineCNA(raw_count_mtx, par_cores = nproc, SUBCLONES = TRUE, plotTree = FALSE, organism = "mouse")

