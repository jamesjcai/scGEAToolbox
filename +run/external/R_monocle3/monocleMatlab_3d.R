computePseudoTime <- function(X, outputFile){
  cMatrix <- read.csv(X, header = FALSE)
  rownames(cMatrix) <- paste0("G", seq_len(nrow(cMatrix)))
  colnames(cMatrix) <- paste0("C", seq_len(ncol(cMatrix)))
  cMatrix <- cMatrix[rowSums(cMatrix) > 0,]
  require(monocle)
  fd <- data.frame('gene_short_name' = rownames(cMatrix))
  rownames(fd) <- rownames(cMatrix)
  fd <- new("AnnotatedDataFrame", data = fd)
  cds <- newCellDataSet(as.matrix(cMatrix), featureData = fd, expressionFamily = negbinomial.size())
  cds <- estimateSizeFactors(cds)
  cds <- reduceDimension(cds, reduction_method = "DDRTree", verbose = TRUE, max_components = 3)
  cds <- orderCells(cds)
  pseudoTiveV <- pData(cds)
  pseudoTiveV <- cbind(pseudoTiveV,t(cds@reducedDimS[,rownames(pseudoTiveV)]))
  pseudoTiveV <- pseudoTiveV[,c(2,4,5,6)]
  colnames(pseudoTiveV) <- c("pseudoTime", "DDRTree1", "DDRTree2", "DDRTree3")
  write.csv(pseudoTiveV, file = outputFile)
}

computePseudoTime(X = "input.csv", outputFile = "output.csv")
