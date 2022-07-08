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
#' @param tree Phylogenetic tree, object of class phylo
#' @param imports Vector containing the indexes of imported tips
#' @param show.axis Whether or not to show the time axis
#' @param colorBase Base branch color
#' @param colorImports Color for imports
#' @param colorDescendants Color for the local descendants of imports
#' @param ... Additional parameters are passed on to the plot.phylo function
#' @return Plot of the tree with colored imports
#' @export
plotImports=function(tree,imports=c(),show.axis=T,colorBase="black",colorImports="red",colorDescendants="blue",...)
{
  if (is.null(tree$stats)) m=keyStats(tree)$stats else m=tree$stats
  ecols=rep(colorBase,Nedge(tree))
  tcols=rep(colorBase,Ntip(tree))
  if (length(imports>0)) {
    wave=c()
    reds=c()
    for (i in 1:length(imports)) {
      a=imports[i]
      ci=m[a,"coalintdiffdate"]
      if (!is.na(ci)) {
        while (ci>0) {
          w=which(tree$edge[,2]==a)
          reds=c(reds,w)
          ci=ci-tree$edge.length[w]
          a=tree$edge[w,1]
          if (length(ci)==0) ci=0
        }
        wave=c(wave,w)
      }
    }
    while (length(wave)>0) {
      new_wave=c()
      for (i in wave) {
        ecols[i]=colorDescendants
        w=which(tree$edge[,1]==tree$edge[i,2])
        for (j in w)
          if (ecols[j]!=colorDescendants) new_wave=c(new_wave,j)
      }
      wave=new_wave
    }
    for (i in reds) ecols[i]=colorImports
    for (i in 1:Nedge(tree)) {
      j=tree$edge[i,2]
      if (j<=Ntip(tree)) tcols[j]=ecols[i]
    }
  }
  args=list(x=tree,show.tip.label=F,edge.color=ecols,tip.color=tcols)
  args=modifyList(args,list(...))
  do.call(plot,args)
  if (show.axis) axisPhylo(1,backward = F)
}

#' Plotting method for DetectImports results
#' @param x Output from DetectImports test methods
#' @param type Type of plot to produce. Can be 'scatter' (default) or 'tree'
#' @param ... Additional parameters are passed on
#' @return Plot of DetectImports results
#' @export
plot.resDetectImports=function(x,type='scatter',...)
{
  stopifnot(inherits(x, "resDetectImports"))
  if (type=='scatter') {
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
    if (!hasArg('xlab')) args=modifyList(args,list(xlab='Sampling dates'))
    if (!hasArg('ylab')) args=modifyList(args,list(ylab='Coalescent intervals'))
    do.call("plot",args)
    ix=sort(dates,index.return=T)$ix
    if (!is.null(x$mus)) lines(dates[ix],x$mus[ix],col='blue')
    if (!is.null(x$mus_low )) lines(dates[ix],x$mus_low [ix],col='blue',lty=2)
    if (!is.null(x$mus_high)) lines(dates[ix],x$mus_high[ix],col='blue',lty=2)
    legend('topleft',legend=c("p>0.01","0.001<p<=0.01","p<=0.001"),col=c("black","red3","red"), lty=1, cex=0.8)
  }
  if (type=='tree') plotImports(x$tree,which(x$pvals<0.01),...)
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
  print(sprintf('Tree with %d leaves. %d imports were found with p<=0.01. Lowest p-value was %.2e.',Ntip(x$tree),length(which(x$pvals<=0.01)),min(x$pvals,na.rm=T)))
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
