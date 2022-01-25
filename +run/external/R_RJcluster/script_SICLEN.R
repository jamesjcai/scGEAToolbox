library(SICLEN)

X <- as.matrix(read.csv("input.csv", header = FALSE))
#X <- as.matrix(read.table("input.csv", sep=","))
res <- siclen(t(X))
write.csv(res$clusterid,'output.csv')
