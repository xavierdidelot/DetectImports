#' Test for the presence of imports assuming that the local population has a constant, known population size
#' @param tree Tree
#' @param Ne Population size
#' @return p-values for importation
#' @export
test0=function(tree,Ne=NULL)
{
  if (is.null(tree$stats)) m=keyStats(tree) else m=tree$stats
  coalints=m[2:Ntip(tree),'coalint']
  if (is.null(Ne)) Ne=median(coalints)/log(2)
  h=hist(coalints,main='',breaks=20,xlab='Coalescent intervals',ylab='Frequency')
  br=h$breaks
  lines(br[-1]-diff(br)[1]/2,length(coalints)*(pexp(br[-1],1/Ne)-pexp(br[-length(br)],1/Ne)))
  pvals=1-pexp(coalints,1/Ne)
  pvals=p.adjust(pvals,"fdr")
  return(pvals)
}

#' Semi-parametric test for the presence of imports
#' @param tree Tree
#' @param epsilon Precision parameter
#' @return p-values for importation
#' @export
test1=function(tree,epsilon)
{
  if (is.null(tree$stats)) m=keyStats(tree) else m=tree$stats
  dates=m[2:Ntip(tree),'dates']
  coalints=m[2:Ntip(tree),'coalint']
  l=length(dates)
  NeHat=rep(NA,l)
  for (i in 1:l) {
    w=setdiff(which(abs(dates-dates[i])<epsilon),i)
    NeHat[i]=median(coalints[w])/log(2)
  }
  plot(dates,coalints,xlab='',ylab='')
  lines(dates,NeHat,col='red')
  pvals=1-pexp(coalints,1/NeHat)
  pvals=p.adjust(pvals,"fdr")
  return(pvals)
}
