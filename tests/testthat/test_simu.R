#Test simulation functions
context("Test simulation functions")

test_that("Coalescent simulation is as expected.", {
  texp=1999.9
  expect_silent(tree<-simCoal(2000:2020,function(t){return(pmax(0,t-texp))},1e-2))
  expect_is(tree,'phylo')
  expect_equal(Ntip(tree),length(2000:2020))
  expect_gt(tree$root.time,texp)
  dates=tree$root.time+unname(dist.nodes(tree)[Ntip(tree)+1,1:Ntip(tree)])
  expect_equal(dates,2000:2020)
})

test_that("Simulation with imports is as expected.",{
  tree=simImports(localPopStart=2020,importDates=c(2020.5),
                  samplingStartDate=2020,samplingEndDate=2021,samplingNumber=100)
  expect_is(tree,'phylo')
  expect_equal(Ntip(tree),100)
  dates=tree$root.time+unname(dist.nodes(tree)[Ntip(tree)+1,1:Ntip(tree)])
  expect_gt(min(dates),2020)
  expect_lt(max(dates),2021)
})
