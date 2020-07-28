#setwd(getSrcDirectory()[1])
#setwd("e:\\GitHub\\scGEAToolbox\\thirdparty\\R_MAST")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("MAST")

library(MAST)
X <- as.matrix(read.table("input1.txt", sep="\t"))
Y <- as.matrix(read.table("input2.txt", sep="\t"))


# write.table(A, paste(dir,"/A.txt",sep=""), row.names=F, col.names=F, sep="\t")

cond<-factor(colData(sca)$condition)
cond<-relevel(cond,"Unstim")
colData(sca)$condition<-cond
zlmCond <- zlm(~condition + cngeneson, sca)
# The following are equivalent
## lrt <- lrTest(zlm, "condition")
## lrt <- lrTest(zlm, CoefficientHypothesis('conditionStim'))

# This would test if 2*cngeneson=conditionStim
#  This is sheer nonsense biologically and statistically, but gives an example of the flexibility.
## lrt <- lrTest(zlm, Hypothesis('2*cngeneson-conditionStim'))
We could run a likelihood ratio test here, testing for differences when we drop the condition factor. Note that any arbitrary contrast matrix can be tested here, and specified either using a matrix or syntactically. See Hypothesis for details.

#only test the condition coefficient.
summaryCond <- summary(zlmCond, doLRT='conditionStim') 
#print the top 4 genes by contrast using the logFC
print(summaryCond, n=4)
