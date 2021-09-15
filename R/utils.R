#' Plot the coalescent intervals
#' @param tree Tree
#' @param ... Additional parameters are passed on
#' @return Plot
#' @export
plotCoalInt=function(tree,...)
{
  if (is.null(tree$stats)) m=keyStats(tree)$stats else m=tree$stats
  dates=m[1:Ntip(tree),'dates']
  coalints=m[1:Ntip(tree),'coalint']
  plot(dates,coalints,
       xlab='Sampling dates',ylab='Coalescent intervals',...)
}

#' Plot the tree and imports
#' @param tree Tree
#' @param imports Imports
#' @param ... Additional parameters are passed on
#' @return Plot
#' @export
plotImports=function(tree,imports,...)
{
  if (is.null(tree$stats)) m=keyStats(tree)$stats else m=tree$stats
  cols=rep('black',Nedge(tree))
  if (length(imports>0)) for (i in 1:length(imports)) {
    a=imports[i]
    ci=m[a,'coalintdiffdate']
    while (ci>0) {
      w=which(tree$edge[,2]==a)
      cols[w]='red'
      ci=ci-tree$edge.length[w]
      a=tree$edge[w,1]
    }
  }
  plot(tree,show.tip.label = F,edge.color = cols,...)
  axisPhylo(1,backward = F)

}

#' Plotting method for DetectImports results
#' @param x Output from DetectImports test methods
#' @param ... Additional parameters are passed on
#' @return Plot of DetectImports results
#' @export
plot.resDetectImports=function(x,...)
{
  stopifnot(inherits(x, "resDetectImports"))
  n=Ntip(x$tree)
  stats=x$tree$stats
  dates=stats[1:n,'dates']
  symbols=rep(1,length(dates))
  cols=rep(1,length(dates))
  cols[which(x$pvals<=0.01 )]=2
  cols[which(x$pvals<=0.001)]=3
  cols=c('black','red3','red')[cols]
  if (!is.null(x$tree$imports)) symbols[x$tree$imports]=2
  args=list(x=dates,y=stats[1:n,'coalint'],pch=symbols,col=cols)
  args=modifyList(args,list(...))
  if (!hasArg(xlab)) args=modifyList(args,list(xlab='Sampling dates'))
  if (!hasArg(ylab)) args=modifyList(args,list(ylab='Coalescent intervals'))
  do.call("plot",args)
  ix=sort(dates,index.return=T)$ix
  if (!is.null(x$mus)) lines(dates[ix],x$mus[ix],col='blue')
  if (!is.null(x$mus_low )) lines(dates[ix],x$mus_low [ix],col='blue',lty=2)
  if (!is.null(x$mus_high)) lines(dates[ix],x$mus_high[ix],col='blue',lty=2)
  legend('topleft',legend=c("p>0.01","0.001<p<=0.01","p<=0.001"),col=c("black","red3","red"), lty=1, cex=0.8)
}

#' Print function for DetectImports results
#' @param x output from DetectImports test methods
#' @param ... Passed on
#' @return Print out details of DetectImports results
#' @export
print.resDetectImports <- function(x, ...)
{
  stopifnot(inherits(x, "resDetectImports"))
  cat('Results of DetectImports test\n')
  print(sprintf('Tree with %d leaves. %d imports were found with p<=0.05. Lowest p-value was %.2e.',Ntip(x$tree),length(which(x$pvals<=0.05)),min(x$pvals,na.rm=T)))
  invisible(x)
}

#' Summary function for DetectImports objects
#' @param object output from DetectImports test methods
#' @param ... Passed on
#' @return Print out details of DetectImports results
#' @export
summary.resDetectImports <- function(object, ...){
  print(object,...)
}
