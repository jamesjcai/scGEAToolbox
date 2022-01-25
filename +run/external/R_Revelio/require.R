if (!requireNamespace("R.matlab", quietly = TRUE)){
    install.packages("R.matlab", repo="http://cran.rstudio.com/")
}

if (!requireNamespace("Revelio", quietly = TRUE)){
devtools::install_github('danielschw188/Revelio')
}

if (!requireNamespace("rhdf5", quietly = TRUE)){
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repo="http://cran.rstudio.com/")
}
BiocManager::install("rhdf5")
}
