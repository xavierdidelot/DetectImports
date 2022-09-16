rm(list=ls())
library(DetectImports)
library(ape)

methods=c('full','constant','fast')
set.seed(0)
repeats=90
tim=array(0,dim=c(repeats,1+length(methods)))
for (i in 1:repeats) {
  ntip=100+10*i
  tim[i,1]=ntip
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
  for (m in 1:length(methods)) {
    ptm <- proc.time()
    if (methods[m]=='full') pvals=detectImports(tree,verbose=F,seed=0)$pvals
    if (methods[m]=='constant') pvals=detectImports(tree,verbose=F,seed=0,constant=T)$pvals
    if (methods[m]=='fast') pvals=detectImportsFAST(tree)$pvals
    ptm = proc.time()-ptm
      tim[i,m+1]=ptm[3]
  }
}

pdf('/tmp/timings.pdf',7,7)
plot(tim[,1],tim[,2],ylim=c(0,60),xlab='Number of genomes',ylab='Time (seconds)')
dev.off()
system('open /tmp/timings.pdf')
