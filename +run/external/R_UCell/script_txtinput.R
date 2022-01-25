setwd("U:\\GitHub\\scGEAToolbox\\+run\\thirdparty\\R_UCell")
library(UCell)
# library(Matrix)
# library(Seurat)

pbmc.counts <- read.table('input.txt', sep = '\t', stringsAsFactors = FALSE)
geneList <- make.unique(pbmc.counts[,1])[-1]
bcList <- pbmc.counts[1,][-1]
pbmc.counts <- Matrix(as.matrix(pbmc.counts[-1,-1]))
rownames(pbmc.counts) <- geneList
colnames(pbmc.counts) <- bcList

markers <- list()
# markers$Tcell_NK <- c("FGFBP2", "SPON2", "KLRF1", "FCGR3A", "KLRD1", "TRDC")
markers$Tcell_exhaustion <-c("CD69","PDCD1","TGFB1","CTLA4","SPN","LAG3")


#pbmc <- CreateSeuratObject(counts = pbmc.counts)
#pbmc <- AddModuleScore_UCell(pbmc, features = markers)
#write.csv(pbmc@meta.data$Tcell_exhaustion_UCell, file = 'output.csv')

scores <- ScoreSignatures_UCell(pbmc.counts, features=markers)
write.csv(scores,file="output.csv")
