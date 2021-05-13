#Test stats computations
context("Test stats computations")

test_that("Coalescent intervals are as expected on small example.", {
  n=10
  t=list()
  t$Nnode=n-1
  t$tip.label=as.character(1:n)
  t$edge=cbind(c((n+1):(2*n-1),2*n-1,(n+1):(2*n-2)),c(1:n,(n+2):(2*n-1)))
  t$edge.length=c(rep(0.5,n),rep(1,n-2))
  class(t)<-'phylo'
  coalint=DetectImports:::coalIntervals(t)[[1]]
  expect_equal(coalint,c(NA,rep(0.5,n-1)))
})

test_that("Changing order of leaves does not affect stats.", {
  n=10
  t=rtree(n)
  s=keyStats(t)$stats
  reorder=sample(1:n,n,replace = F)
  w=which(t$edge[,2]<=n)
  t$edge[w,2]=reorder[t$edge[w,2]]
  s2=keyStats(t)$stats
  expect_equal(s,s2[c(reorder,(n+1):nrow(s2)),])
})
