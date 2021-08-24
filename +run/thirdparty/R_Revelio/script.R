#setwd("C:\\Users\\jcai.AUTH\\Documents\\GitHub\\scGEAToolbox\\+run\\thirdparty\\R_Revelio")
#setwd("u:\\GitHub\\scGEAToolbox\\+run\\thirdparty\\R_Revelio")
library(Revelio)
library(R.matlab)

#sc <- read.table('input.csv', sep = ',', stringsAsFactors = FALSE)
mat <-readMat('input.mat')
pbmc.counts <- mat$X
rownames(pbmc.counts) <- make.unique(unlist(mat$genelist))


# revelioTestData_rawDataMatrix   pbmc.counts

myData <- createRevelioObject(rawData = pbmc.counts,
                              cyclicGenes = revelioTestData_cyclicGenes)

myData <- getCellCyclePhaseAssignInformation(dataList = myData)
#myData <- getPCAData(dataList = myData)
#myData <- getOptimalRotation(dataList = myData)

myData <- getPCAData(dataList = myData, boolPlotResults = FALSE)
myData <- getOptimalRotation(dataList = myData, boolPlotResults = FALSE)

# normalizedDataWithoutCCEffects <- removeCCEffects(dataList = myData)


#res <- decontX(x = sc)
#X=res$decontXcounts
#contamination=res$contamination
# ccPhase=myData@cellInfo$ccPhase
dcdata=as.matrix(myData@transformedData$dc$data[c('DC1','DC2'),])
writeMat("output.mat", dc=dcdata)
# write.table(as.matrix(X),file="output.csv", sep=",",col.names=FALSE,row.names = FALSE)
a<-myData@cellInfo$ccPhase
b<-myData@cellInfo$cellID
write.table(list(a,b),file="output.csv", sep=",",col.names=FALSE,row.names = FALSE)
# write.table(,file="output2.csv", sep=",",col.names=FALSE,row.names = FALSE)