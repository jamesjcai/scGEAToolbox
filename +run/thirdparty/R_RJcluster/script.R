library(RJcluster)

X <- as.matrix(read.csv("input.csv", header = FALSE))
#X <- as.matrix(read.table("input.csv", sep=","))
#res = RJclust(data = X, scale = TRUE)
res = RJclust(data = X)
write.csv(res$class,'output.csv')
