suppressMessages(library(CooccurrenceAffinity))
library(rhdf5)

X <- h5read(file = "input.h5", name = "/X")
myout <- affinity(data = X, row.or.col = "row", squarematrix = c("p_value"))
h5createFile("output.h5")
h5write(as.matrix(myout$p_value), "output.h5","/p")
