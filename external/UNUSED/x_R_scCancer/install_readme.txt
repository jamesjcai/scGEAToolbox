For R (version>=4.0) under Windows system, you need to update the Rtools to version 4.0 from
https://cran.r-project.org/bin/windows/Rtools/
 
Installing Rtools44
https://cran.r-project.org/bin/windows/base/howto-R-devel.html
install.packages(devtools)
devtools::install_github("PolMine/RcppCWB")
 

https://github.com/wguo-research/scCancer/wiki/2.-Installation

checkPkg <- function(pkg){
    return(requireNamespace(pkg, quietly = TRUE))
}
if(!checkPkg("BiocManager")) install.packages("BiocManager")
if(!checkPkg("devtools")) install.packages("devtools")

library(devtools)
if(!checkPkg("harmony")) install_github("immunogenomics/harmony")
if(!checkPkg("RcppArmadillo")) install.packages("RcppArmadillo")
if(!checkPkg("RcppProgress")) install.packages("RcppProgress")
if(!checkPkg("NNLM")) install_github("linxihui/NNLM")
install.packages("rliger")


rstudio searching all files
replace 'liger' with 'rliger'
install.packages("D:/downloads/scCancer/", repos = NULL, type = "source")
