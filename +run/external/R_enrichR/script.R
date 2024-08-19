library(enrichR)
#listEnrichrSites()
#> Enrichr ... Connection is Live!
#> FlyEnrichr ... Connection is available!
#> WormEnrichr ... Connection is available!
#> YeastEnrichr ... Connection is available!
#> FishEnrichr ... Connection is available!
#> OxEnrichr ... Connection is available!
setEnrichrSite("Enrichr") # Human genes
#> Connection changed to https://maayanlab.cloud/Enrichr/
#> Connection is Live!
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")
enriched <- enrichr(c("Runx1", "Gfi1", "Gfi1b", "Spi1", "Gata1", "Kdr"), dbs)
#> Uploading data to Enrichr... Done.
#>   Querying GO_Molecular_Function_2015... Done.
#>   Querying GO_Cellular_Component_2015... Done.
#>   Querying GO_Biological_Process_2015... Done.
#> Parsing results... Done.
enriched[["GO_Biological_Process_2015"]]
plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")


