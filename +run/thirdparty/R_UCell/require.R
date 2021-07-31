if (!requireNamespace("UCell", quietly = TRUE)){
    if (!requireNamespace("remotes", quietly = TRUE)){
        install.packages("remotes", repo="http://cran.rstudio.com/")
    }
    library(remotes)
    remotes::install_github("carmonalab/UCell")
}