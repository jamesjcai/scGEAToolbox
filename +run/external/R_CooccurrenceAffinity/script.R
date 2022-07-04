suppressMessages(library(CooccurrenceAffinity))
library(rhdf5)

# data(finches)


X <- h5read(file = "input.h5", name = "/X")
# g <- h5read(file = "input.h5", name = "/g")
myout <- affinity(data = X, row.or.col = "row", squarematrix = c("all"))
