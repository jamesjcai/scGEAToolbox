library(rhdf5)
library(harmony)

#data(cell_lines)
#V <- cell_lines$scaled_pcs
#meta_data <- cell_lines$meta_data
#harmony_embeddings <- harmony::RunHarmony(
#  V, meta_data, 'dataset', verbose=FALSE
#)

s <- h5read(file = "input.h5", name = "/s")
batchid <- h5read(file = "input.h5", name = "/batchid")


meta_data <- data.frame(
  batchid = as.character(batchid)
)


harmony_embeddings <- RunHarmony(
  s, meta_data, 'batchid', verbose = FALSE
)


#h5create("output.h5", "/harmony_embeddings", dim(harmony_embeddings))
#h5write(harmony_embeddings, "output.h5", "/harmony_embeddings")

h5createFile("output.h5")
h5write(as.matrix(harmony_embeddings), "output.h5", "/harmony_embeddings")

