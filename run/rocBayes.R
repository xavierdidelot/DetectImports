#This is a ROC curve for testBayes
rm(list=ls())
library(DetectImports)
library(ape)

ntip=500
pvalres=1000
roc=roc_test1=matrix(0,pvalres,2)
repeats=10
for (i in 1:repeats) {
  print(i)
  set.seed(i)
  tree=simImports(localPopStart=2020,importDates=c(2020.25,2020.5),
                samplingStartDate=2020,samplingEndDate=2021,samplingNumber=ntip)
  real=tree$imports#real imports
  unreal=setdiff(1:ntip,real)
  pvals=suppressMessages(testBayes(tree))
  for (j in 1:pvalres) {
    p=j/pvalres
    roc[j,1]=roc[j,1]+length(which(pvals[  real]<p))/length(real)
    roc[j,2]=roc[j,2]+length(which(pvals[unreal]<p))/length(unreal)
  }
  pvals=suppressMessages(test1(tree,showPlot = F))
  for (j in 1:pvalres) {
    p=j/pvalres
    roc_test1[j,1]=roc_test1[j,1]+length(which(pvals[  real]<p))/length(real)
    roc_test1[j,2]=roc_test1[j,2]+length(which(pvals[unreal]<p))/length(unreal)
  }
}
roc=roc/repeats
plot(roc[,2],roc[,1],type='l',xlab ='FPR',ylab='TPR')
points(roc[pvalres*c(0.01,0.05),2],roc[pvalres*c(0.01,0.05),1])
roc_test1=roc_test1/repeats
lines(roc_test1[,2],roc_test1[,1],type='l',xlab ='FPR',ylab='TPR',col='red')
points(roc_test1[pvalres*c(0.01,0.05),2],roc_test1[pvalres*c(0.01,0.05),1],col='red')
legend('bottomright',legend=c("testBayes", "test1"),col=c("black", "red"), lty=1, cex=0.8)
