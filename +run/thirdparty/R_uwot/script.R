computeUmap <- function(X, outputFile){
  cMatrix <- read.csv(X, header = FALSE)
  require(uwot)
  s_umap <- umap(cMatrix, n_components = 2, n_neighbors = 50, learning_rate = 0.5, init = "random")
  write.csv(s_umap, file = outputFile)
}
computeUmap(X = "input.csv", outputFile = "output.csv")
