if (!requireNamespace("Matrix", quietly = TRUE)){
    install.packages("Matrix", repo="http://cran.rstudio.com/")
}
if (!requireNamespace("Seurat", quietly = TRUE)){
    install.packages("Seurat", repo="http://cran.rstudio.com/")
}

if (!requireNamespace("rhdf5", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
}
BiocManager::install("rhdf5")
}

checkPkg <- function(pkg){
    return(requireNamespace(pkg, quietly = TRUE))
}
if(!checkPkg("devtools")) install.packages("devtools")
library(devtools)
if(!checkPkg("scplotter")) devtools::install_github("pwwang/scplotter")
