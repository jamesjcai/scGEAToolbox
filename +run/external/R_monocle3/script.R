require(monocle3)
require(rhdf5)
#setwd('D:\\GitHub\\scGEAToolbox\\+run\\external\\R_monocle3')
X <- h5read(file = "input.h5", name = "/X")
cMatrix <-as.matrix(X)
rownames(cMatrix) <- paste0("G", seq_len(nrow(cMatrix)))
colnames(cMatrix) <- paste0("C", seq_len(ncol(cMatrix)))
cMatrix <- cMatrix[rowSums(cMatrix) > 0,]

fd <- data.frame('gene_short_name' = rownames(cMatrix))
rownames(fd) <- rownames(cMatrix)
fd <- new("AnnotatedDataFrame", data = fd)
cds <- newCellDataSet(as.matrix(cMatrix), featureData = fd, expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- reduceDimension(cds, reduction_method = "DDRTree", verbose = TRUE, max_components = 3)
cds <- orderCells(cds)
pseudoTiveV <- pData(cds)
pseudoTiveV <- cbind(pseudoTiveV,t(cds@reducedDimS[,rownames(pseudoTiveV)]))
pseudoTiveV <- pseudoTiveV[,c(2,4,5,6)]
colnames(pseudoTiveV) <- c("pseudoTime", "DDRTree1", "DDRTree2", "DDRTree3")
write.csv(pseudoTiveV, file = "output.csv")



# https://cole-trapnell-lab.github.io/monocle3/docs/getting_started/

#The expression value matrix must:
#have the same number of columns as the cell_metadata has rows.
#have the same number of rows as the gene_metadata has rows.
#Additionally:
#row names of the cell_metadata object should match the column names of the expression matrix.
#row names of the gene_metadata object should match row names of the expression matrix.
#one of the columns of the gene_metadata should be named "gene_short_name", which represents the gene symbol or simple name (generally used for plotting) for each gene.

# Load the data
expression_matrix <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_expression.rds"))
cell_metadata <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_colData.rds"))
gene_annotation <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_rowData.rds"))


cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)
## Step 2: Remove batch effects with cell alignment
#cds <- align_cds(cds, alignment_group = "batch")
## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)

## Step 4: Cluster the cells
cds <- cluster_cells(cds)
## Step 5: Learn a graph
cds <- learn_graph(cds)

## Step 6: Order cells
cds <- order_cells(cds)

# With regression:
gene_fits <- fit_models(cds, model_formula_str = "~embryo.time")
fit_coefs <- coefficient_table(gene_fits)
emb_time_terms <- fit_coefs %>% filter(term == "embryo.time")
emb_time_terms <- emb_time_terms %>% mutate(q_value = p.adjust(p_value))
sig_genes <- emb_time_terms %>% filter (q_value < 0.05) %>% pull(gene_short_name)

# With graph autocorrelation:
pr_test_res <- graph_test(cds,  neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(pr_test_res, q_value < 0.05))




