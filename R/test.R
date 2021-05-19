#' Test for the presence of imports assuming that the local population has a constant population size
#' @param tree Tree
#' @param Ne Population size (default is estimated from data)
#' @param adjust Method for adjusting p-values (default is fdr)
#' @param showPlot Whether to show a plot of the test
#' @return p-values for importation
#' @export
test0=function(tree,Ne,adjust='fdr',showPlot=T)
{
  if (is.null(tree$stats)) m=keyStats(tree)$stats else m=tree$stats
  coalints=m[1:Ntip(tree),'coalint']
  if (missing(Ne)) Ne=median(coalints,na.rm = T)/log(2)
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
    k=0
    w=c()
    while (length(w)<5) {#just in case there are not enough leaves within epsilon
      k=k+1
      w=setdiff(which(abs(dates-dates[i])<epsilon*k),i)
    }
    NeHat[i]=median(coalints[w],na.rm=T)/log(2)
  }
  if (showPlot) {
    plot(dates,coalints,xlab='',ylab='')
    ix=sort(dates,index.return=T)$ix
    lines(dates[ix],NeHat[ix],col='red')
  }
  pvals=1-pexp(coalints,1/NeHat)
  pvals=p.adjust(pvals,adjust)
  message(sprintf('%d imports were found with p<0.05. Lowest p-value was %.2e',length(which(pvals<0.05)),min(pvals,na.rm=T)))
  return(pvals)
}

#' ESD test for the presence of imports
#' @param tree Tree
#' @return p-values for importation
#' @export
test2=function(tree)
{
  if (is.null(tree$stats)) m=keyStats(tree)$stats else m=tree$stats
  dates=m[1:Ntip(tree),'dates']
  coalints=m[1:Ntip(tree),'coalint']
  Ne=median(coalints,na.rm = T)/log(2)
  #residuals=(coalints-Ne)/sqrt(Ne)#bad normalization
  residuals=qnorm(pexp(coalints,1/Ne))#transform Exp(1/Ne) into N(0,1)
  plot(dates,residuals)
  qqnorm(residuals)
  pvals=1-pnorm(coalints)
  pvals=p.adjust(pvals,'fdr')
  message(sprintf('%d imports were found with p<0.05. Lowest p-value was %.2e',length(which(pvals<0.05)),min(pvals,na.rm=T)))
  return(pvals)
}
