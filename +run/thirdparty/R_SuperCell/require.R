if (!requireNamespace("RANN", quietly = TRUE)){
    install.packages("RANN", repo="http://cran.rstudio.com/")
}

if (!requireNamespace("igraph", quietly = TRUE)){
    install.packages("igraph", repo="http://cran.rstudio.com/")
}

if (!requireNamespace("WeightedCluster", quietly = TRUE)){
    install.packages("WeightedCluster", repo="http://cran.rstudio.com/")
}

if (!requireNamespace("corpcor", quietly = TRUE)){
    install.packages("corpcor", repo="http://cran.rstudio.com/")
}

if (!requireNamespace("weights", quietly = TRUE)){
    install.packages("weights", repo="http://cran.rstudio.com/")
}

if (!requireNamespace("Hmisc", quietly = TRUE)){
    install.packages("Hmisc", repo="http://cran.rstudio.com/")
}

if (!requireNamespace("Matrix", quietly = TRUE)){
    install.packages("Matrix", repo="http://cran.rstudio.com/")
}

if (!requireNamespace("patchwork", quietly = TRUE)){
    install.packages("patchwork", repo="http://cran.rstudio.com/")
}

if (!requireNamespace("plyr", quietly = TRUE)){
    install.packages("plyr", repo="http://cran.rstudio.com/")
}

if (!requireNamespace("irlba", quietly = TRUE)){
    install.packages("irlba", repo="http://cran.rstudio.com/")
}

if (!requireNamespace("remotes")) install.packages("remotes", repo="http://cran.rstudio.com/")
if (!requireNamespace("SuperCell")) remotes::install_github("GfellerLab/SuperCell")

if (!requireNamespace("R.matlab", quietly = TRUE)){
    install.packages("R.matlab", repo="http://cran.rstudio.com/")
}

