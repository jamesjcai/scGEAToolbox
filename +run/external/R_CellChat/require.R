if (!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools", repo="http://cran.rstudio.com/")
}
if (!requireNamespace("NMF", quietly = TRUE)){
    install.packages("NMF", repo="http://cran.rstudio.com/")
}
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)){
    devtools::install_github("jokergoo/ComplexHeatmap")
}
if (!requireNamespace("circlize", quietly = TRUE)){
    devtools::install_github("jokergoo/circlize")
}
if (!requireNamespace("CellChat", quietly = TRUE)){
    devtools::install_github("sqjin/CellChat")
}

if (!requireNamespace("rhdf5", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
}
BiocManager::install("rhdf5")
}
if (!requireNamespace("Matrix", quietly = TRUE)){
    install.packages("Matrix", repo="http://cran.rstudio.com/")
}
library(CellChat)
library(patchwork)
library(Matrix)
library(rhdf5)

