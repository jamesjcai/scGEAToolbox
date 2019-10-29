fw_glmPCA <- function(X, outputFile){
  require(glmpca)
  y <- read.csv(X, header = FALSE)
  temp <- glmpca(Y = y,L = 2, fam = 'nb')
  betaCoefficient <- (temp$coefX[,1])
  pValue <- pnorm(scale(betaCoefficient))
  FDR <- p.adjust(pValue, method = 'fdr')
  out <- data.frame(featureCoefficient = betaCoefficient, p.value = pValue, p.adj = FDR)
  rownames(out) <- rownames(X)
  write.csv(out, file = outputFile)
}
fw_glmPCA(X = "input.csv", outputFile = "output.csv")
