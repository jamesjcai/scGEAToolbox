library(Seurat)
library(Matrix)
library(rhdf5)

filename<-readLines("inputrdsfile.txt")

A<-readRDS(filename)

# counts_matrix = GetAssayData(seurat_obj, slot="counts")
# see: https://github.com/broadinstitute/inferCNV/wiki/infercnv-10x

#gname<-rownames(A@assays$RNA)
#if (is.null(gname)){
#    gname<-rownames(A@assays$RNA@counts)
#}

#gname <- rownames(GetAssayData(A, layer = "counts"))
#write.csv(gname,file='g.csv')


#tryCatch(
#{
#    X=A@assays$RNA@counts
#    #X is a dgCMatrix
#    X<-as.matrix(X)   # this requires large memory
#    h5write(X, "output.h5", "X")
#},
#error = function(msg){
#	# write.csv(A@assays$RNA@counts, file = 'X.csv', col.names=F)  #
#	write.table(A@assays$RNA@counts, file = 'X.csv', sep=",", col.names=FALSE, row.names=FALSE)
#	}
#)

tryCatch(
{
    #X=A@assays$RNA@counts
    #X=A@assays$RNA$counts
    X<-GetAssayData(A, layer="counts")

if (file.exists("output.h5")) {
  file.remove("output.h5")
}

#    h5createDataset(file = ex_file, dataset = "counts_chunked", 
#                    dims = dim(m1), storage.mode = "integer", 
#                    chunk = c(100,100), level = 6)
#    h5write(obj = m1, file = ex_file, name = "counts_chunked")

    h5createFile("output.h5")
    # X is a dgCMatrix
    h5write(X@x, "output.h5", "data")
    h5write(X@i, "output.h5", "indices")
    h5write(X@p, "output.h5", "indptr")
    h5write(X@Dim, "output.h5", "shape")

    # Define chunk size
    #chunk_size_x <- c(1000)  # Adjust the chunk size as needed
    #h5write(X_matrix, "output.h5", "data", chunk = chunk_size_x)

    h5closeAll()

    gname <- rownames(X)
    write.csv(gname,file='g.csv')

},
error = function(msg){

	}
)




# bcode<-colnames(A@assays$RNA@counts)

tryCatch(
    {
        if (!is.null(rownames(A@meta.data))) {
         write.csv(rownames(A@meta.data), file='barcodes.csv')
        }
    },
    error = function(e){ 
        # (Optional)
        # Do this if an error is caught...
    }
)

tryCatch(
    {
        if (!is.null(A@reductions$umap@cell.embeddings)) {
            write.csv(A@reductions$umap@cell.embeddings, file = 'umap.csv')
        }
    },
    error = function(e){ 
        # Do this if an error is caught...
    }
)


#tryCatch(
#    {
#         write.csv(A@meta.data$annotation, file = 'annotation.csv')
#    },
#    error = function(e){ 
        # (Optional)
        # Do this if an error is caught...
#    }
#)


tryCatch(
    {
        if (!is.null(A@meta.data$BatchID)) {
            write.csv(A@meta.data$BatchID, file = 'batch.csv')
        }else if (!is.null(A@meta.data$orig.ident)) {
          write.csv(A@meta.data$orig.ident, file = 'batch.csv')
        }
    },
    error = function(e){ 
          # Do this if an error is caught...
    }
)

tryCatch(
    {
        if (!is.null(A@meta.data$CellType)) {
            write.csv(A@meta.data$CellType, file = 'celltype.csv')
        }
    },
    error = function(e){ 
        # Do this if an error is caught...
    }
)

tryCatch(
  {
    if (!is.null(A@meta.data$celltype)) {
      write.csv(A@meta.data$celltype, file = 'celltype.csv')
    }
  },
  error = function(e){ 
    # Do this if an error is caught...
  }
)


tryCatch(
    {
        if (!is.null(A@meta.data)) {
            write.csv(A@meta.data, file = 'metadata.csv')
        }
    },
    error = function(e){ 
        # Do this if an error is caught...
    }
)