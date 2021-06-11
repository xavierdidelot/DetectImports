rm(list=ls())
library(DetectImports)
library(ape)

ntip=1000
fpr=tpr=0
repeats=1
alpha=1
for (i in 1:repeats) {
  print(i)
  set.seed(i)
  tree=simImports(localPopStart=2020,importDates=c(2020.25,2020.5),
                samplingStartDate=2020,samplingEndDate=2021,samplingNumber=ntip)
  real=tree$imports#real imports
  unreal=setdiff(1:ntip,real)
  pvals=test2(tree,showPlot = T,alpha=alpha,maxi=100)
  tpr=tpr+length(which(pvals[  real]==alpha))/length(real)
  fpr=fpr+length(which(pvals[unreal]==alpha))/length(unreal)
}
fpr=fpr/repeats
tpr=tpr/repeats
