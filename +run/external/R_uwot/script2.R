computeUmap <- function(X, outputFile){
  cMatrix <- read.csv(X, header = FALSE)
  require(umap)
  custom.config = umap.defaults
  custom.config$n_components=3
  custom.config$n_neighbors=50
  s_umap <- umap(cMatrix, custom.config)
  write.csv(s_umap, file = outputFile)
}
computeUmap(X = "input2.csv", outputFile = "output2.csv")
