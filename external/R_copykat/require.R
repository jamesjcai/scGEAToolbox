#setwd(getSrcDirectory()[1])
checkPkg <- function(pkg){
    return(requireNamespace(pkg, quietly = TRUE))
}
if(!checkPkg("devtools")) install.packages("devtools")
library(devtools)
if(!checkPkg("copykat")) install_github("navinlabcode/copykat")

if (!requireNamespace("rhdf5", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
}
BiocManager::install("rhdf5")
}
