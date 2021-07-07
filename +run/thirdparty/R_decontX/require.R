if (!requireNamespace("celda", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}
BiocManager::install("celda")
}
if (!requireNamespace("R.matlab", quietly = TRUE)){
    install.packages("R.matlab")
}