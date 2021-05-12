#' Test for the presence of imports assuming that the local population has a constant, known population size
#' @param tree Tree
#' @param Ne Population size
#' @param adjust Method for adjusting p-values (default is fdr)
#' @param showPlot Whether to show a plot of the test
#' @return p-values for importation
#' @export
test0=function(tree,Ne=NULL,adjust='fdr',showPlot=T)
{
  if (is.null(tree$stats)) m=keyStats(tree)$stats else m=tree$stats
  coalints=m[1:Ntip(tree),'coalint']
  if (is.null(Ne)) Ne=median(coalints,na.rm = T)/log(2)
  if (showPlot) {
    h=hist(coalints,main='',breaks=20,xlab='Coalescent intervals',ylab='Frequency')
    br=h$breaks
    lines(br[-1]-diff(br)[1]/2,length(coalints)*(pexp(br[-1],1/Ne)-pexp(br[-length(br)],1/Ne)))
  }
  pvals=1-pexp(coalints,1/Ne)
  pvals=p.adjust(pvals,adjust)
  message(sprintf('%d imports were found with p<0.01. Lowest p-value was %.2e',length(which(pvals<0.01)),min(pvals,na.rm=T)))
  return(pvals)
}

#' Semi-parametric test for the presence of imports
#' @param tree Tree
#' @param epsilon Precision parameter
#' @param adjust Method for adjusting p-values (default is fdr)
#' @param showPlot Whether to show a plot of the test
#' @return p-values for importation
#' @export
test1=function(tree,epsilon,adjust='fdr',showPlot=T)
{
  if (is.null(tree$stats)) m=keyStats(tree)$stats else m=tree$stats
  dates=m[1:Ntip(tree),'dates']
  coalints=m[1:Ntip(tree),'coalint']
  if (missing(epsilon)) epsilon=(max(dates)-min(dates))/20
  l=length(dates)
  NeHat=rep(NA,l)
  for (i in 1:l) {
    w=setdiff(which(abs(dates-dates[i])<epsilon),i)
    NeHat[i]=median(coalints[w],na.rm=T)/log(2)
  }
  if (showPlot) {
    plot(dates,coalints,xlab='',ylab='')
    lines(dates,NeHat,col='red')
  }
  pvals=1-pexp(coalints,1/NeHat)
  pvals=p.adjust(pvals,adjust)
  message(sprintf('%d imports were found with p<0.05. Lowest p-value was %.2e',length(which(pvals<0.05)),min(pvals,na.rm=T)))
  return(pvals)
}
