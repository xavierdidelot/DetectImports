#' Plot the coalescent intervals
#' @param tree Tree
#' @return Plot
#' @export
plotCoalInt=function(tree)
{
if (is.null(tree$stats)) m=keyStats(tree) else m=tree$stats
dates=m[2:Ntip(tree),'dates']
coalints=m[2:Ntip(tree),'coalint']
plot(dates,coalints,
     xlab='Sampling dates',ylab='Coalescent intervals')
}
