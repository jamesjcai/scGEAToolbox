doit <- function(inputFile, outputFile){
#BiocManager::install('GSVA')
#install.packages('clustermole', dependencies = TRUE)
library(clustermole)
library(fgsea)

# Building the database for Mouse
CT <- clustermole_markers()
CT <- CT[CT$species %in% 'Mouse',]
CT <- CT[CT$db %in% c('CellMarker', 'PanglaoDB'),]
ctNames <- unique(CT$celltype)
CT <-lapply(ctNames, function(X){
  unique(CT$gene[CT$celltype %in% X])
})
names(CT) <- ctNames

# Reading List
geneList <- toupper(readLines(inputFile))

# Creating a rank (is prefered to use the FC here)
gList <- sort(rnorm(length(geneList)), decreasing = TRUE)

# Assigning the ID
names(gList) <- geneList

# Running the cell type assignation
E <- fgseaMultilevel(CT, gList)

# Sorting the result
E <- E[order(E$NES,1/E$padj, decreasing = TRUE),]

# Result
head(E)
write.csv(E$pathway, file = outputFile,row.names = FALSE)

}
doit(inputFile = "input.txt", outputFile = "output.txt")
