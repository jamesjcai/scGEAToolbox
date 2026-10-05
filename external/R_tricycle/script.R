# Cell cycle position with tricycle (Zheng et al. 2022, Genome Biology).
# Called by run.r_tricycle, which writes input.h5, species.txt and, when a
# reference is given, reference.h5 into the working folder.
suppressMessages({
    library(rhdf5)
    library(Matrix)
    library(SingleCellExperiment)
    library(scuttle)
    library(tricycle)
})

# e_writeh5 stores a sparse matrix as 1-based CSC: /data, /indices (rows),
# /indptr (column starts) and /shape.
read_counts <- function(f) {
    shape <- as.integer(h5read(f, "/shape"))
    m <- sparseMatrix(i = as.integer(h5read(f, "/indices")),
                      p = as.integer(h5read(f, "/indptr")) - 1L,
                      x = as.numeric(h5read(f, "/data")),
                      dims = shape)
    rownames(m) <- as.character(h5read(f, "/g"))
    m
}

species <- readLines("species.txt", warn = FALSE)[1]
X <- read_counts("input.h5")
nCells <- ncol(X)
useRef <- file.exists("reference.h5")
counts <- X
if (useRef) {
    # Normalised together, so sample and reference share one size-factor scale.
    counts <- cbind(X, read_counts("reference.h5"))
}
colnames(counts) <- paste0("C", seq_len(ncol(counts)))

sce <- SingleCellExperiment(list(counts = counts))
sce <- logNormCounts(sce)
sce <- project_cycle_space(sce, gname.type = "SYMBOL", species = species)
emb <- reducedDim(sce, "tricycleEmbedding")

# project_cycle_space centres each gene on the mean of all the cells it is
# given, so the angle depends on the sample's make-up. With a reference, the
# angle is taken around the reference cells' mean instead: shifting the
# embedding by their mean is the same as centring the genes on it.
centre <- c(0, 0)
if (useRef) {
    centre <- colMeans(emb[-seq_len(nCells), , drop = FALSE])
}
sce <- estimate_cycle_position(sce, center.pc1 = centre[1], center.pc2 = centre[2])

keep <- seq_len(nCells)
out <- data.frame(position = sce$tricyclePosition[keep],
                  pc1 = emb[keep, 1] - centre[1],
                  pc2 = emb[keep, 2] - centre[2])
write.csv(out, "output.csv", row.names = FALSE)
