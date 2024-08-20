library(enrichR)
MS <- read.csv('input.txt')
#write.csv(MS, 'output.txt')
#listEnrichrSites()

setEnrichrSite("Enrichr") # Human genes
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2023", "GO_Biological_Process_2023")
#enriched <- enrichr(c("Runx1", "Gfi1", "Gfi1b", "Spi1", "Gata1", "Kdr"), dbs)
enriched <- enrichr(c(MS$genelist), dbs)

write.csv(enriched[["GO_Biological_Process_2023"]],'output1.txt')
write.csv(enriched[["GO_Molecular_Function_2023"]],'output2.txt')

# plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")


