#BiocManager::install(c('SingleR', 'scater'))
library(SingleR)
library(scater)

# Reference dataset
hpca.se <- HumanPrimaryCellAtlasData()

# New dataset
X <- read.csv('input.txt', sep = '\t', row.names = 1)
X <- SummarizedExperiment(list(counts = X))
X <- logNormCounts(X)

# Cell types
CT <- SingleR(test = X, ref = hpca.se, labels = hpca.se$label.main)
writeLines(CT$labels, 'output.csv')
