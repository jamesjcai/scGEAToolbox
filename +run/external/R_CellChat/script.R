library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Matrix)
library(rhdf5)


X <- h5read(file = "input.h5", name = "/X")
g <- h5read(file = "input.h5", name = "/g")
countMatrix <- Matrix(as.matrix(X))
rownames(countMatrix) <- g
colnames(countMatrix) <- paste0("C", seq_len(ncol(countMatrix)))
countMatrix <- CreateSeuratObject(countMatrix)
sce <- NormalizeData(countMatrix)
saveRDS(sce, file = "output.Rds")



cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
# showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use



X <- h5read(file = "input.h5", name = "/X")
g <- h5read(file = "input.h5", name = "/g")

#countMatrix <- read.table('input.txt', sep = '\t', stringsAsFactors = FALSE)
#geneList <- make.unique(countMatrix[,1])[-1]
#bcList <- countMatrix[1,][-1]

countMatrix <- Matrix(as.matrix(X))
rownames(countMatrix) <- g
colnames(countMatrix) <- paste0("C", seq_len(ncol(countMatrix)))

countMatrix <- CreateSeuratObject(countMatrix)
sce <- NormalizeData(countMatrix)
saveRDS(sce, file = "output.Rds")

#countMatrix <- ScaleData(countMatrix, features = unlist(cc.genes.updated.2019))
#countMatrix <- CellCycleScoring(countMatrix, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
#countMatrix <- data.frame(S=countMatrix$S.Score, G2M=countMatrix$G2M.Score, Phase = countMatrix$Phase)
#write.csv(countMatrix, file = 'output.csv')

https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html

