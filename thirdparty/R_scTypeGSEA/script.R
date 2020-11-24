library(Seurat)

countMatrix <- read.csv('counts.txt', header = FALSE)
rownames(countMatrix) <- readLines('genelist.txt')
colnames(countMatrix) <- paste0('c', seq_len(ncol(countMatrix)))

countMatrix <- CreateSeuratObject(countMatrix)
Idents(countMatrix) <- readLines('clusterid.txt')

getLog2FC <- function(X){
  require(Matrix)
  allIdent <- levels(Idents(X))
  countData <- X@assays$RNA@counts
  countData <- (t(t(countData)/colSums(countData)))*1e4
  lapply(allIdent, function(i){
    iMean <- rowMeans(countData[,Idents(X) %in% i])
    oMean <- rowMeans(countData[,!(Idents(X) %in% i)])
    sGenes <- ((iMean > 0) & (oMean > 0))
    iMean <- iMean[sGenes]
    oMean <- oMean[sGenes]
    FC <- iMean/oMean
    log2FC <- log2(FC)
    return(log2FC)
  })
}

FC <- getLog2FC(countMatrix)

library(clustermole)
CT <- clustermole_markers(species = 'mm')
CT <- CT[grepl('Intestin|GI',CT$organ, ignore.case = TRUE),]
CT$celltype <- (gsub('s$','',CT$celltype))
ctList <- unique(CT$celltype)
CT <- lapply(ctList, function(cellType){
  unique(CT$gene[CT$celltype %in% cellType])
})
names(CT) <- ctList

library(fgsea)
E <- lapply(FC, function(X){
  # E <- fgseaMultilevel(CT, X, eps = 0)
  E <- fgseaMultilevel(CT, X)
  E <- E[E$NES > 0,]
  E <- E[order(E$NES, decreasing = TRUE),]
  E[1,, drop= FALSE]
})
O <- do.call(rbind.data.frame, E)

O <- data.frame(cluster = levels(Idents(countMatrix)), O)
# O$leadingEdge <- unlist(lapply(O$leadingEdge, function(X){paste0(X, collapse = ';')}))
# write.csv(O, file = 'cellAssignation.csv')

levels(Idents(countMatrix)) <- O$pathway
writeLines(as.vector(Idents(countMatrix)), 'output.txt')
