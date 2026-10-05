if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repo = "http://cran.rstudio.com/")
}
for (p in c("rhdf5", "SingleCellExperiment", "scuttle", "tricycle")) {
    if (!requireNamespace(p, quietly = TRUE)) {
        BiocManager::install(p, update = FALSE, ask = FALSE)
    }
}
if (!requireNamespace("Matrix", quietly = TRUE)) {
    install.packages("Matrix", repo = "http://cran.rstudio.com/")
}
