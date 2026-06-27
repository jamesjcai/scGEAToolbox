require(CoGAPS)
require(rhdf5)

# ---- Read inputs written by MATLAB (save -v7.3 => HDF5) --------------------
X         <- h5read("input.mat", "/X")            # genes x cells (numeric)
nPatterns <- as.integer(h5read("input.mat", "/nPatterns"))
nIter     <- as.integer(h5read("input.mat", "/nIterations"))
sparseOpt <- as.logical(h5read("input.mat", "/sparseOpt"))
genes     <- readLines("genes.txt")               # real gene symbols (one per line)

M <- as.matrix(X)
rownames(M) <- genes
colnames(M) <- paste0("C", seq_len(ncol(M)))

# ---- Run CoGAPS (Bayesian NMF) --------------------------------------------
params <- CogapsParams(nPatterns = nPatterns,
                       nIterations = nIter,
                       seed = 42,
                       sparseOptimization = sparseOpt)
res <- CoGAPS(M, params)

A <- res@featureLoadings      # genes x patterns  (A matrix)
P <- res@sampleFactors        # cells x patterns  (P matrix)

# ---- PatternMarker statistic ----------------------------------------------
# patternMarkers() assigns each gene to the pattern it is most specific to.
pm <- patternMarkers(res, threshold = "all")
pmList <- pm$PatternMarkers
pnames <- names(pmList)
if (is.null(pnames)) {
    pnames <- paste0("Pattern_", seq_along(pmList))
}
markerTbl <- do.call(rbind, lapply(seq_along(pmList), function(i) {
    genesI <- as.character(pmList[[i]])
    if (length(genesI) == 0) {
        return(NULL)
    }
    data.frame(gene = genesI,
               pattern = pnames[i],
               rank = seq_along(genesI),
               stringsAsFactors = FALSE)
}))
write.csv(markerTbl, "output_markers.csv", row.names = FALSE)

# ---- Write factorization matrices -----------------------------------------
if (file.exists("output.h5")) {
    file.remove("output.h5")
}
h5write(as.matrix(A), "output.h5", "A")
h5write(as.matrix(P), "output.h5", "P")
