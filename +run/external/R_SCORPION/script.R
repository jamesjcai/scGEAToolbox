library(Matrix)
library(SCORPION)
library(rhdf5)

load("hg38_PPI.RData")
load("hg38_TF.RData")

# hsaPPI hsaTF
X <- h5read(file = "input.h5", name = "/X")
g <- h5read(file = "input.h5", name = "/g")

countMatrix <- Matrix(as.matrix(X))
rownames(countMatrix) <- g
colnames(countMatrix) <- paste0("C", seq_len(ncol(countMatrix)))

scorpionOutput <- scorpion(tfMotifs = hsaTF,
                           gexMatrix = as(countMatrix, "dgCMatrix"),
                           ppiNet = hsaPPI,
                           alphaValue = 0.1)

if (file.exists("output.h5")) {
  file.remove("output.h5")
}

h5_file <- "output.h5"
h5createFile(h5_file)
h5write(as.matrix(scorpionOutput$regNet), "regNet", file = h5_file)
h5write(colnames(scorpionOutput$regNet),"colnames", file = h5_file)
h5write(rownames(scorpionOutput$regNet),"rownames", file = h5_file)
h5closeAll()

# h5write(as.matrix(scorpionOutput$regNet), 'filename.h5', "regNet")
# b = h5read('filename.h5','/regNet');

