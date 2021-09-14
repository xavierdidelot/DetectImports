#This is a ROC curve for test1
rm(list=ls())
library(DetectImports)
library(ape)

ntip=500
pvalres=1000
roc=rocOnline=matrix(0,pvalres,2)
repeats=10
for (i in 1:repeats) {
  print(i)
  set.seed(i)
  tree=simImports(localPopStart=2020,importDates=c(2020.25,2020.5),
                samplingStartDate=2020,samplingEndDate=2021,samplingNumber=ntip)
  real=tree$imports#real imports
  unreal=setdiff(1:ntip,real)
  pvals=suppressMessages(test1(tree)$pvals)
  for (j in 1:pvalres) {
    p=j/pvalres
    roc[j,1]=roc[j,1]+length(which(pvals[  real]<p))/length(real)
    roc[j,2]=roc[j,2]+length(which(pvals[unreal]<p))/length(unreal)
  }
  pvals=suppressMessages(test1(tree,online=T)$pvals)
  for (j in 1:pvalres) {
    p=j/pvalres
    rocOnline[j,1]=rocOnline[j,1]+length(which(pvals[  real]<p))/length(real)
    rocOnline[j,2]=rocOnline[j,2]+length(which(pvals[unreal]<p))/length(unreal)
  }
}
roc=roc/repeats
plot(roc[,2],roc[,1],type='l',xlab ='FPR',ylab='TPR')
points(roc[pvalres*c(0.01,0.05),2],roc[pvalres*c(0.01,0.05),1])
rocOnline=rocOnline/repeats
lines(rocOnline[,2],rocOnline[,1],type='l',xlab ='FPR',ylab='TPR',col='red')
points(rocOnline[pvalres*c(0.01,0.05),2],rocOnline[pvalres*c(0.01,0.05),1],col='red')
legend('bottomright',legend=c("Offline", "Online"),col=c("black", "red"), lty=1, cex=0.8)
