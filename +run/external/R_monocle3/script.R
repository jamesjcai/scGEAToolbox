require(monocle3)
require(rhdf5)

#X <- h5read(file = "input.h5", name = "/X")

X <- h5read(file = "input.mat", name = "/X")
idx <- h5read(file = "input.mat", name = "/idx")

cMatrix <-as.matrix(X)
rownames(cMatrix) <- paste0("G", seq_len(nrow(cMatrix)))
colnames(cMatrix) <- paste0("C", seq_len(ncol(cMatrix)))
#cMatrix <- cMatrix[rowSums(cMatrix) > 0,]

g <- data.frame('gene_short_name' = rownames(cMatrix))
rownames(g) <- rownames(cMatrix)
c <- data.frame('cell_short_name' = colnames(cMatrix))
rownames(c) <- colnames(cMatrix)

cds <- new_cell_data_set(cMatrix,
                         cell_metadata = c,
                         gene_metadata = g)
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

cds <- order_cells(cds, root_cells = paste0("C", idx))
# cds <- order_cells(cds, root_cells = c("C20", "C21"))
ps_tim <- pseudotime(cds)
# https://rdrr.io/github/cole-trapnell-lab/monocle3/man/pseudotime.html
h5write(ps_tim, "output.h5", "t")

# X <- h5read(file = "input.h5", name = "/X")
