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

#' Plot the tree and imports
#' @param tree Tree
#' @param imports Imports
#' @return Plot
#' @export
plotImports=function(tree,imports)
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
  plot(tree,show.tip.label = F,edge.color = cols)
  axisPhylo(1,backward = F)

}
