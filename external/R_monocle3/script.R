require(monocle3)
require(rhdf5)

#X <- h5read(file = "input.h5", name = "/X")

X <- h5read(file = "input.mat", name = "/X")
idx <- h5read(file = "input.mat", name = "/idx")
ndim <- h5read(file = "input.mat", name = "/ndim")

cMatrix <-as.matrix(X)
rownames(cMatrix) <- paste0("G", seq_len(nrow(cMatrix)))
colnames(cMatrix) <- paste0("C", seq_len(ncol(cMatrix)))
#cMatrix <- cMatrix[rowSums(cMatrix) > 0,]

g <- data.frame('gene_short_name' = rownames(cMatrix))
rownames(g) <- rownames(cMatrix)
cell_meta <- data.frame('cell_short_name' = colnames(cMatrix))
rownames(cell_meta) <- colnames(cMatrix)

cds <- new_cell_data_set(cMatrix,
                         cell_metadata = cell_meta,
                         gene_metadata = g)
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, max_components = ndim)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

cds <- order_cells(cds, root_cells = paste0("C", idx))
# cds <- order_cells(cds, root_cells = c("C20", "C21"))
ps_tim <- pseudotime(cds)
# https://rdrr.io/github/cole-trapnell-lab/monocle3/man/pseudotime.html
ss_mat <- SingleCellExperiment::reducedDims(cds)[["UMAP"]]
dp_mst <- cds@principal_graph_aux[["UMAP"]]$dp_mst

# Check if the file exists and remove it if it does
if (file.exists("output.h5")) {
  file.remove("output.h5")
}

h5write(ps_tim, "output.h5", "t")
h5write(ss_mat, "output.h5", "s")
h5write(as.matrix(igraph::as_adjacency_matrix(dp_mst)), "output.h5", "m")

ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
q <- ciliated_cds_pr_test_res$q_value
h5write(q, "output.h5", "q")

# X <- h5read(file = "input.h5", name = "/X")
