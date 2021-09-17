library(RJcluster)

X <- as.matrix(read.csv("input.csv", header = FALSE))
#X <- as.matrix(read.table("input.csv", sep=","))
#res = RJclust(data = X, scale = TRUE)
res = RJclust(data = X, n_bins = 20, scaleRJ = TRUE, C_max = 20)
write.csv(res$class,'output.csv')
