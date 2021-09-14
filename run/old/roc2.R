#This is a ROC curve for test2
rm(list=ls())
library(DetectImports)
library(ape)

ntip=500
pvalres=100
roc=matrix(0,pvalres,2)
repeats=10
for (i in 1:repeats) {
  print(i)
  set.seed(i)
  tree=simImports(localPopStart=2020,importDates=c(2020.25,2020.5),
                samplingStartDate=2020,samplingEndDate=2021,samplingNumber=ntip)
  real=tree$imports#real imports
  unreal=setdiff(1:ntip,real)
  for (j in 1:pvalres) {
    p=j/pvalres
    pvals=suppressMessages(test2(tree,showPlot = F,alpha=p))
    roc[j,1]=roc[j,1]+length(which(pvals[  real]<=p))/length(real)
    roc[j,2]=roc[j,2]+length(which(pvals[unreal]<=p))/length(unreal)
  }
}
roc=roc/repeats
plot(roc[,2],roc[,1],type='l',xlab ='FPR',ylab='TPR')
points(roc[pvalres*c(0.01,0.05),2],roc[pvalres*c(0.01,0.05),1])
