#BiocManager::install(c('SingleR', 'scater'))
library(SingleR)
library(scater)
library(celldex)

# Reference dataset
hpca.se <- celldex::HumanPrimaryCellAtlasData()
# ImmGenData BlueprintEncodeData MonacoImmuneData MouseRNAseqData NovershternHematopoieticData

# New dataset
X <- read.csv('input.txt', sep = '\t', row.names = 1, na.string=".")
X <- SummarizedExperiment(list(counts = X))
X <- logNormCounts(X)

# Cell types
CT <- SingleR(test = X, ref = hpca.se, labels = hpca.se$label.main)
writeLines(CT$labels, 'output.csv')
