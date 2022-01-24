#dir <- dirname(parent.frame(2)$ofile)
#setwd(dir)
#BiocManager::install(c('SingleR', 'scater'))
library(SingleR)
library(scater)
# source('https://raw.githubusercontent.com/dosorio/utilities/master/idConvert/hsa2mmu_SYMBOL.R')
library(celldex)

hsa2mmu_SYMBOL <- function(Z){
  Z <- sapply(Z, function(I){
    I <- unlist(strsplit(I,''))
    I <- c(I[1],tolower(I[-1]))
    I <- paste0(I, collapse = "")
    return(I)
  })
  Z <- as.character(Z)
  return(Z)
}

# Reference dataset
refdata <- celldex::MouseRNAseqData()
# refdata <-ImmGenData()
# HumanPrimaryCellAtlasData ImmGenData BlueprintEncodeData MonacoImmuneData MouseRNAseqData NovershternHematopoieticData

# New dataset
X <- read.csv('input.txt', sep = '\t', row.names = 1, na.string=".")
rownames(X) <- hsa2mmu_SYMBOL(rownames(X))
X <- SummarizedExperiment(list(counts = X))
X <- logNormCounts(X)

# Cell types
CT <- SingleR(test = X, ref = refdata, labels = refdata$label.main)
writeLines(CT$labels, 'output.csv')
