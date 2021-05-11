#' Plot the coalescent intervals
#' @param tree Tree
#' @return Plot
#' @export
plotCoalInt=function(tree)
{
if (is.null(tree$stats)) m=keyStats(tree)$stats else m=tree$stats
dates=m[1:Ntip(tree),'dates']
coalints=m[1:Ntip(tree),'coalint']
plot(dates,coalints,
     xlab='Sampling dates',ylab='Coalescent intervals')
}
