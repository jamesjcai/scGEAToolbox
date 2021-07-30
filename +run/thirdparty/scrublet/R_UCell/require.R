if (!requireNamespace("UCell", quietly = TRUE)){
    if (!requireNamespace("remotes", quietly = TRUE)){
        install.packages("remotes")
    }
    library(remotes)
    remotes::install_github("carmonalab/UCell")
}