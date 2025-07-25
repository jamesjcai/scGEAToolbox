library(scCancer)

dataPath <- "./data/KC-example"     # The path of cell ranger processed data
savePath <- "./results/KC-example"  # A path to save the results
sampleName <- "KC-example"          # The sample name
authorName <- "G-Lab@THU"           # The author name to mark the report

# Run scStatistics
stat.results <- runScStatistics(
    dataPath = dataPath,
    savePath = savePath,
    sampleName = sampleName,
    authorName = authorName
)

# Run scAnnotation
anno.results <- runScAnnotation(
    dataPath = dataPath,
    statPath = statPath,
    savePath = savePath,
    authorName = authorName,
    sampleName = sampleName,
    geneSet.method = "average"       # or "GSVA"
)