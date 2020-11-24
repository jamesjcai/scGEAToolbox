library(Seurat)

countMatrix <- read.csv('counts.txt', header = FALSE)
rownames(countMatrix) <- readLines('genelist.txt')
colnames(countMatrix) <- paste0('c', seq_len(ncol(countMatrix)))

countMatrix <- CreateSeuratObject(countMatrix)
Idents(countMatrix) <- readLines('clusterid.txt')

DE <- FindAllMarkers(countMatrix)
write.csv(DE, 'output.txt')
