suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
# suppressMessages(require(R.matlab))

#setwd('U:\\GitHub\\scGEAToolbox\\+run\\thirdparty\\R_SeuratReadRds')
#setwd("C:\\Users\\jcai.AUTH\\Documents\\GitHub\\scGEAToolbox\\+run\\thirdparty\\R_SeuratReadRds")
# countMatrix <- read.table('input.txt', sep = '\t', stringsAsFactors = FALSE)
# fileName <- 'input.txt'
# f <- readChar(fileName, file.info(fileName)$size-2)

filename<-readLines("inputrdsfile.txt")
A<-readRDS(filename)

# counts_matrix = GetAssayData(seurat_obj, slot="counts")
# see: https://github.com/broadinstitute/inferCNV/wiki/infercnv-10x

tryCatch(
{writeMat("output.mat",X=A@assays$RNA@counts)},
error = function(msg){
	# write.csv(A@assays$RNA@counts, file = 'X.csv', col.names=F)  #
	write.table(A@assays$RNA@counts, file = 'X.csv', sep=",", col.names=FALSE, row.names=FALSE)
	}
	
)
write.csv(rownames(A@assays$RNA),file='g.csv')
#write.csv(A@reductions$umap@cell.embeddings, file = 'umap.csv')
#write.csv(A$orig.ident, file = 'batchid.csv');

#geneList <- make.unique(countMatrix[,1])[-1]
#bcList <- countMatrix[1,][-1]
#countMatrix <- Matrix(as.matrix(countMatrix[-1,-1]))
#rownames(countMatrix) <- geneList
#colnames(countMatrix) <- bcList
#countMatrix <- CreateSeuratObject(countMatrix)
#sce <- NormalizeData(countMatrix)
#saveRDS(sce, file = "output.Rds")

#countMatrix <- ScaleData(countMatrix, features = unlist(cc.genes.updated.2019))
#countMatrix <- CellCycleScoring(countMatrix, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
#countMatrix <- data.frame(S=countMatrix$S.Score, G2M=countMatrix$G2M.Score, Phase = countMatrix$Phase)
#write.csv(countMatrix, file = 'output.csv')

