#BiocManager::install(c('SingleR', 'scater'))
library(SingleR)
library(scater)

# Reference dataset
refdata <- MouseRNAseqData()
# HumanPrimaryCellAtlasData ImmGenData BlueprintEncodeData MonacoImmuneData MouseRNAseqData NovershternHematopoieticData

# New dataset
X <- read.csv('input.txt', sep = '\t', row.names = 1, na.string=".")
X <- SummarizedExperiment(list(counts = X))
X <- logNormCounts(X)

# Cell types
CT <- SingleR(test = X, ref = refdata, labels = refdata$label.main)
writeLines(CT$labels, 'output.csv')
