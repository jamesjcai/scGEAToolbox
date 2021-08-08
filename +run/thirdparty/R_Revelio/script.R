setwd("C:\\Users\\jcai.AUTH\\Documents\\GitHub\\scGEAToolbox\\+run\\thirdparty\\R_Revelio")
# C:\Users\jcai.AUTH\Documents\GitHub\scGEAToolbox\+run\thirdparty
library(Revelio)
library(R.matlab)

#sc <- read.table('input.csv', sep = ',', stringsAsFactors = FALSE)
mat <-readMat('input.mat')
pbmc.counts <- mat$X
rownames(pbmc.counts) <- make.unique(unlist(mat$genelist))


# revelioTestData_rawDataMatrix   

myData <- createRevelioObject(rawData = pbmc.counts,
                              cyclicGenes = revelioTestData_cyclicGenes)
myData <- getCellCyclePhaseAssignInformation(dataList = myData)
myData <- getPCAData(dataList = myData)
myData <- getOptimalRotation(dataList = myData)
normalizedDataWithoutCCEffects <- removeCCEffects(dataList = myData)


#res <- decontX(x = sc)
#X=res$decontXcounts
#contamination=res$contamination
# ccPhase=myData@cellInfo$ccPhase
# writeMat("output.mat", X=as.matrix(normalizedDataWithoutCCEffects))
# write.table(as.matrix(X),file="output.csv", sep=",",col.names=FALSE,row.names = FALSE)

write.table(myData@cellInfo$ccPhase,file="output1.csv", sep=",",col.names=FALSE,row.names = FALSE)
write.table(myData@cellInfo$cellID,file="output2.csv", sep=",",col.names=FALSE,row.names = FALSE)