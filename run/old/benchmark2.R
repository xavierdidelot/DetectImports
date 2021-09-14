#Computes and compares FPR and TPR for all tests
rm(list=ls())
library(DetectImports)
library(ape)

ntip=500
fpr=tpr=c(0,0,0)
repeats=10
alpha=0.05
for (i in 1:repeats) {
  print(i)
  set.seed(i)
  tree=simImports(localPopStart=2020,importDates=c(2020.25,2020.5),
                samplingStartDate=2020,samplingEndDate=2021,samplingNumber=ntip)
  real=tree$imports#real imports
  unreal=setdiff(1:ntip,real)
  pvals=suppressMessages(test1(tree)$pvals)
  tpr[1]=tpr[1]+length(which(pvals[  real]<=alpha))/length(real)
  fpr[1]=fpr[1]+length(which(pvals[unreal]<=alpha))/length(unreal)
  pvals=suppressMessages(test1(tree,online=T)$pvals)
  tpr[2]=tpr[2]+length(which(pvals[  real]<=alpha))/length(real)
  fpr[2]=fpr[2]+length(which(pvals[unreal]<=alpha))/length(unreal)
  pvals=suppressMessages(test2(tree,alpha=alpha,maxi=100)$pvals)
  tpr[3]=tpr[3]+length(which(pvals[  real]<=alpha))/length(real)
  fpr[3]=fpr[3]+length(which(pvals[unreal]<=alpha))/length(unreal)
}
fpr=fpr/repeats
tpr=tpr/repeats
results=rbind(fpr,tpr)
colnames(results)<-c('test1','test1Online','test2')
print(results)
