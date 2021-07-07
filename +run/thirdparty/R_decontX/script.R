library(celda)

sc <- read.table('input.csv', sep = ',', stringsAsFactors = FALSE)
sc <- as.matrix(sc)
pbmc4k <- decontX(x = sc)

dsc<-data.matrix(pbmc4k$decontXcounts)
write.table(dsc,file="output.csv", sep=",",col.names=FALSE,row.names = FALSE)
