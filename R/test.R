#' Test for the presence of imports assuming that the local population has a constant, known population size
#' @param dates Dates of samples
#' @param coalints Coalescent intervals
#' @param Ne Population size
#' @return p-values for importation
#' @export
test0=function(dates,coalints,Ne)
{
  h=hist(coalints,main='',breaks=20)
  br=h$breaks
  lines(br[-1]-diff(br)[1]/2,length(coalints)*(pexp(br[-1],1/Ne)-pexp(br[-length(br)],1/Ne)))
  pvals=pexp(coalints,1/Ne)
  pvals=p.adjust(pvals,"fdr")
  return(pvals)
}
