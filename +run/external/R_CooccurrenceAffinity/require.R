if (!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools", repo="http://cran.rstudio.com/")
}
if (!requireNamespace("CooccurrenceAffinity", quietly = TRUE)){
    require(devtools)
    install_github("kpmainali/CooccurrenceAffinity")
}

if (!requireNamespace("rhdf5", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
}
BiocManager::install("rhdf5")
}