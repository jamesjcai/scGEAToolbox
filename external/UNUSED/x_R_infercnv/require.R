#setwd(getSrcDirectory()[1])
checkPkg <- function(pkg){
    return(requireNamespace(pkg, quietly = TRUE))
}
if(!checkPkg("BiocManager")) install.packages("BiocManager")
if(!checkPkg("infercnv")) BiocManager::install("infercnv")


if (!requireNamespace("infercnv", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
}
    BiocManager::install("infercnv")
}

if (!requireNamespace("rhdf5", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
}
    BiocManager::install("rhdf5")
}


if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
BiocManager::install("infercnv")




checkPkg <- function(pkg){
    return(requireNamespace(pkg, quietly = TRUE))
}
if(!checkPkg("BiocManager")) install.packages("BiocManager")
if(!checkPkg("devtools")) install.packages("devtools")

library(devtools)
if(!checkPkg("harmony")) install_github("immunogenomics/harmony")

if(!checkPkg("RcppArmadillo")) install.packages("RcppArmadillo")
if(!checkPkg("RcppProgress")) install.packages("RcppProgress")
if(!checkPkg("NNLM")) install_github("linxihui/NNLM")

if(!checkPkg("liger")) install_github("MacoskoLab/liger")

install_github("wguo-research/scCancer")
