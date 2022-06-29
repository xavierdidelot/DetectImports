#This is a ROC curve comparison between methods
rm(list=ls())
library(DetectImports)
library(ape)

methods=c('full','constant','fast')
set.seed(0)
ntip=500
pvalres=10000
roc=array(0,dim=c(pvalres,2,length(methods)))
repeats=100
for (i in 1:repeats) {
  print(i)
  set.seed(i)
  tree=NULL
  while(is.null(tree)) {try(
    tree<-simImports(localPopStart=2020,importDates=c(2020.25,2020.5),
                    samplingStartDate=2020,samplingEndDate=2021,samplingNumber=ntip)
  ,silent=T)
  if (is.null(tree)) next
  d=dist.nodes(tree)[Ntip(tree)+1,1:Ntip(tree)]
  if (is.element(which(d==min(d)),tree$imports)) tree=NULL#force first leaf to be from non-imported local pop
  }
  real=tree$imports#real imports
  unreal=setdiff(1:ntip,real)
  for (m in 1:length(methods)) {
    if (methods[m]=='full') pvals=detectImports(tree,verbose=F,seed=0)$pvals
    if (methods[m]=='constant') pvals=detectImports(tree,verbose=F,seed=0,constant=T)$pvals
    if (methods[m]=='fast') pvals=detectImportsFAST(tree)$pvals
    for (j in 1:pvalres) {
      p=j/pvalres
      roc[j,1,m]=roc[j,1,m]+length(which(pvals[  real]<p))/length(real)
      roc[j,2,m]=roc[j,2,m]+length(which(pvals[unreal]<p))/length(unreal)
    }
  }
}
roc=roc/repeats
save.image('roc.RData')

for (p in c(0.05,0.01,0.001)) for (m in 1:length(methods))
  print(sprintf('With cutoff p=%.3f, %s method has sensitivity=%.3f and specificity=%.3f',p,methods[m],roc[pvalres*p,1,m],1-roc[pvalres*p,2,m]))

#Plotting
methods=c('Model with variable population size','Model with constant population size')
cols=c('red','blue')
pdf('/tmp/roc.pdf',7,7)
plot(c(0,1),c(0,1),type='l',lty=2,xlab ='False Positive Rate (1-Specificity)',ylab='True Positive Rate (Sensitivity)')
for (m in 1:length(methods)) {
  lines(c(0,roc[,2,m],1),c(0,roc[,1,m],1),col=cols[m])
  points(roc[pvalres*0.01,2,m],roc[pvalres*0.01,1,m],col=cols[m])
}
legend('bottomright',legend=methods,col=cols, lty=1, cex=0.8)
dev.off()
system('open /tmp/roc.pdf')

m=1
sum(roc[,1,m]*diff(c(0,roc[,2,m])))

m=2
sum(roc[,1,m]*diff(c(0,roc[,2,m])))

