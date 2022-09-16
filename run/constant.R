rm(list=ls())
library(DetectImports)
library(ape)

methods=c('full','constant','fast')
set.seed(0)
repeats=100
falsepos=array(0,dim=c(repeats,length(methods)))
for (i in 1:repeats) {
  print(i)
  set.seed(i)
  tree<-simCoal(dates=runif(500,2020,2021),NeFun = function(x){return(1)})
  for (m in 1:length(methods)) {
    ptm <- proc.time()
    if (methods[m]=='full') pvals=detectImports(tree,verbose=F,seed=0)$pvals
    if (methods[m]=='constant') pvals=detectImports(tree,verbose=F,seed=0,constant=T)$pvals
    if (methods[m]=='fast') pvals=detectImportsFAST(tree)$pvals
    ptm = proc.time()-ptm
    falsepos[i,m]=length(which(pvals<0.01))
  }
}

pdf('/tmp/falsepos.pdf',7,7)
boxplot(falsepos[,c(1,2)]/500,names=c('Variable population size model','Constant population size model'),notch=T,ylim=c(0,0.02),ylab='Proportion of false positives')
dev.off()
system('open /tmp/falsepos.pdf')
