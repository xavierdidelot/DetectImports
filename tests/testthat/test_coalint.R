#Test coalescent interval function
context("Test coalescent interval function")

test_that("Coalescent intervals are as expected on small example.", {
  n=10
  t=list()
  t$Nnode=n-1
  t$tip.label=as.character(1:n)
  t$edge=cbind(c((n+1):(2*n-1),2*n-1,(n+1):(2*n-2)),c(1:n,(n+2):(2*n-1)))
  t$edge.length=c(rep(0.5,n),rep(1,n-2))
  class(t)<-'phylo'
  coalint=DetectImports:::coalIntervals(t)
  expect_equal(coalint,c(NA,rep(0.5,n-1)))
})
