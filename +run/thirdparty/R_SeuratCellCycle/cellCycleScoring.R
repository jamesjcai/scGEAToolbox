library(Seurat)

features <- list(c(
    'CD79B',
    'CD79A',
    'CD19',
    'CD180',
    'CD200',
    'CD3D',
    'CD2',
    'CD3E',
    'CD7',
    'CD8A',
    'CD14',
    'CD1C',
    'CD68',
    'CD9',
    'CD247'
  ))

object <- pbmc_small
ctrl = 5
pool = NULL
nbin = 24
k = FALSE
assay = NULL
name = 'Cluster'
seed = 1
search = FALSE

set.seed(seed = seed)
cluster.length <- length(x = features)
object <- object@assays$RNA@data
pool <- rownames(x = object)
data.avg <- Matrix::rowMeans(x = object)
data.avg <- data.avg[order(data.avg)]
data_cut <- function(data.avg, nbin){
  n.points <- length(data.avg)
  assigned.bin <- rep(FALSE, n.points)
  bin.size <- n.points/nbin
  for(i in seq_len(nbin)){
    bin.match <- data.avg <= data.avg[round(bin.size*(i))]
    pos.avail <- assigned.bin == FALSE
    assigned.bin[(bin.match & pos.avail)] <- i
  }
  return(assigned.bin)
}

data.cut <- data_cut(data.avg, nbin)

names(x = data.cut) <- names(x = data.avg)
ctrl.use <- vector(mode = "list", length = cluster.length)
for (i in 1:cluster.length) {
  features.use <- features[[i]]
  for (j in 1:length(x = features.use)) {
    ctrl.use[[i]] <- c(
      ctrl.use[[i]],
      names(x = sample(
        x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
        size = ctrl,
        replace = FALSE
      ))
    )
  }
}
ctrl.use <- lapply(X = ctrl.use, FUN = unique)

ctrl.scores <- matrix(
  data = numeric(length = 1L),
  nrow = length(x = ctrl.use),
  ncol = ncol(x = object)
)
for (i in 1:length(ctrl.use)) {
  features.use <- ctrl.use[[i]]
  ctrl.scores[i, ] <- Matrix::colMeans(x = object)
}
features.scores <- matrix(
  data = numeric(length = 1L),
  nrow = cluster.length,
  ncol = ncol(x = object)
)
for (i in 1:cluster.length) {
  features.use <- features[[i]]
  data.use <- object[features.use, , drop = FALSE]
  features.scores[i, ] <- Matrix::colMeans(x = data.use)
}
features.scores.use <- features.scores - ctrl.scores
rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
features.scores.use <- as.data.frame(x = t(x = features.scores.use))
rownames(x = features.scores.use) <- colnames(x = object)
