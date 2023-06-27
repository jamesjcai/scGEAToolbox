rm(list=ls())
library(amap)

homedir = '/data/hoan/'
datasets = c('camp1','grun','li','patel', 'pollen', 'wang', 'xin', 'zeisel')
for (dataset in datasets) {
    data_dir = paste0(homedir, 'spectral_clustering_matlab/data/',dataset ,'-prepare-log_count.csv')
    data <- read.table(data_dir, sep=',')
    data <- data[,2:ncol(data)]
    distance_matrix <- as.matrix(Dist(data, method = "euclidean", nbproc = 32, diag = FALSE, upper = T))
    write.table(distance_matrix, file = paste0(homedir, 'spectral_clustering_matlab/data/',dataset ,'_distance_matrix.csv'),
                sep=',', col.names = F, row.names = F)
}


