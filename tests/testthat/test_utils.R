#Test utility functions
context("Test utility functions")

test_that("Utility functions work as expected.", {
  set.seed(0)
  NeFun=function(t){return(5)}
  tree=simCoal(seq(2000,2020,1),NeFun)
  expect_silent(plotImports(tree,c(2,3)))
  expect_silent(plotCoalInt(tree))
  expect_silent(plotBoth(tree,NeFun))
  suppressWarnings(res<-detectImports(tree,verbose=F,iter=100,nchains=1))
  expect_is(capture_output(print(res)),'character')
  expect_is(capture_output(summary(res)),'character')
  expect_silent(plot(res))
})

