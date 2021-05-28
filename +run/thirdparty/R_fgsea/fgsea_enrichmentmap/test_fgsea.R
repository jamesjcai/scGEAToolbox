library(fgsea)
library(data.table)
library(ggplot2)

setwd("C:\\Users\\jcai\\Desktop\\fgsea_enrichmentmap")
GOiea <- gmtPathways('Supplementary_Table3_Human_GOBP_AllPathways_no_GO_iea_July_01_2017_symbol.gmt')
df <- read.csv('Supplementary_Table2_MesenvsImmuno_RNASeq_ranks.rnk', sep="\t", header=T)
Z=df$rank
names(Z) <- df$GeneName
fgseaRes <- fgsea(pathways = GOiea, stats = Z, minSize  = 15, maxSize  = 500)
#head(fgseaRes[order(pval), ])
#write.csv(fgseaRes, 'output.txt')

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(GOiea[topPathways], Z, fgseaRes, gseaParam=0.5)


plotEnrichment(GOiea[["ALZHEIMER DISEASE-PRESENILIN PATHWAY%PANTHER PATHWAY%P00004"]],
               Z) + labs(title="ALZHEIMER DISEASE-PRESENILIN PATHWAY%PANTHER PATHWAY%P00004")