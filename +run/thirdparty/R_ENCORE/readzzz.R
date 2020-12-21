readzzz <- function(countMatrix, geneList){
  X <- read.csv(countMatrix, header = FALSE)
  rownames(X) <- readLines('genelist.txt')
  colnames(X) <- paste0('C', seq_len(ncol(X)))
  X <- t(X)
  X <- as.data.frame(X)
  return(X)
} 
