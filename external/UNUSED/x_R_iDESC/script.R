setwd('D:\\GitHub\\scGEAToolbox\\+run\\external\\R_iDESC')
library(Matrix)
library(rhdf5)
library(iDESC)


cMatrix <- h5read(file = "input.mat", name = "/X")
# ndim <- h5read(file = "input.mat", name = "/ndim")

a <- read.csv('a.txt')
g <- read.table('g.txt')
rownames(cMatrix) <- t(g)
# colnames(cMatrix) <- paste0("C", seq_len(ncol(cMatrix)))
colnames(cMatrix) <- a$CellID

#metadata <- cbind(a$BatchID, a$CellType, a$sequencing_depth)
#metadata <- as.data.frame(metadata)
#colnames(metadata) <- c('BatchID', 'CellType', 'sequencing_depth')
#rownames(metadata) <- a$CellID

metadata <- data.frame(
  BatchID = a$BatchID,
  CellType = a$CellType,
  sequencing_depth = a$sequencing_depth,
  row.names = a$CellID,
  stringsAsFactors = FALSE
)


mat=cMatrix
sequencing_depth=a$sequencing_depth

result=iDESC(mat,metadata,subject_var="CellType",group_var="BatchID",
             norm_opt="User",user_sf = sequencing_depth,span = 0.7)

write.csv(result, "output.csv", row.names = FALSE)