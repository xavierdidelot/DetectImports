#Test all is running without error
context("Test running without error")

test_that("Basic functions are running without error.", {
  set.seed(0)
  expect_silent(tree<-simCoal(2000:2020,function(t){return(5)}))
  expect_is(tree,'phylo')
  expect_silent(suppressMessages(pvals<-test0(tree,showPlot=F)))
  expect_gt(min(pvals,na.rm = T),0.05)
  expect_silent(suppressMessages(pvals<-test1(tree,epsilon=10,showPlot=F)))
  expect_gt(min(pvals,na.rm = T),0.05)
})
