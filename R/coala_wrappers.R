# Note: MK made these functions simply for learning about and testing applications of the coala package.
# I wouldn't intend them to be distributed as part of this package, although perhaps they could be useful if we generalise them a bit, add tests etc.

#' Simulate and plot trees using the coala package
#' @param pop1_sample_size The size of the sample of the local population (population 1)
#' @param pop2_sample_size The size of the sample of the non-local population (population 2)
#' @param rate Migration rate from population 2 to 1 (default = 0.5)
#' @param divergence.time Time of divergence of population 2 from population 1 (default = 5)
#' @param loci_number The number of loci to add - see `coala::coal_model` (default = 100)
#' @param loci_length The length of the loci to add - see `coala::coal_model` (default = 1000)
#' @param ploidy The number of chromosomes that will be simulated per individual - see `coala::coal_model` (default = 1)
#' @param Ne1 For specifying the effective population size of population 1. At the moment it's implemented quite crudely: if (as default) it is NA then the effective population size is constant over time, if not, it will change abruptly to the proportion specified here at the time given by `Ne1.time.change`
#' @param Ne2 as above but for population 2
#' @param Ne1.time.change as above
#' @param Ne2.time.change as above
#' @param seed Option to set the seed (default NA)
#' @param epsilon Precision parameter to pass to DetectImports::test1 (default 1)
#' @param alpha The threshold for statistical significance (default 0.05)
#' @param ... Extra options to pass to the ape plot.phylo function
#' @return Plot
#' @export
coalaSimThenDetectAndPlot <- function(pop1_sample_size, pop2_sample_size, rate=0.5, divergence.time=5, loci_number=100, loci_length=1000, ploidy=1, Ne1=NA, Ne2=NA, Ne1.time.change, Ne2.time.change, seed=NA, epsilon=1, alpha=0.05, ...) {

ifelse(is.na(seed), sampleno <- sample(1:999, 1), sampleno <- seed) # use supplied random seed or create random sample number here

# prepare models
modelGlobalSample <- coal_model(sample_size = c(pop1_sample_size, pop2_sample_size), loci_number = loci_number, loci_length=loci_length, ploidy=ploidy) +
  feat_migration(rate = rate, symmetric = FALSE, pop_from = 2, pop_to = 1) + # symmetric migration rate of 0.5 between populations
  #feat_pop_merge(0.5, 3, 2) + # population 3 merges into population 2 0.5 coalescent time units in the past
  feat_pop_merge(divergence.time, 2, 1) + # population 2 merges into / diverges from population 1 0.8 coalescent time units in the past (0.3 further into the past)
  sumstat_trees() # + output trees

modelLocalSample <- coal_model(sample_size = c(pop1_sample_size, 0), loci_number = loci_number, loci_length=loci_length, ploidy=ploidy) +
  feat_migration(rate = rate, symmetric = FALSE, pop_from = 2, pop_to = 1) + # symmetric migration rate of 0.5 between populations
  #feat_pop_merge(0.5, 3, 2) + # population 3 merges into population 2 0.5 coalescent time units in the past
  feat_pop_merge(divergence.time, 2, 1) + # population 2 merges into / diverges from population 1 0.8 coalescent time units in the past (0.3 further into the past)
  sumstat_trees() # + output trees

if (!is.na(Ne1) && !is.na(Ne2)) {
  modelGlobalSample <- modelGlobalSample + feat_size_change(new_size=Ne1, time=Ne1.time.change, population=1)
  modelLocalSample <- modelLocalSample + feat_size_change(new_size=Ne2, time=Ne2.time.change, population=1)
}

# simulate
sumstatsGlobal <- simulate(modelGlobalSample, seed = sampleno)
sumstatsLocal <- simulate(modelLocalSample, seed = sampleno)

# extract trees
GlobalTree <- read.tree(text=sumstatsGlobal$trees[[1]])
LocalTree <- read.tree(text=sumstatsLocal$trees[[1]])

pvals=test1(LocalTree,epsilon=1,showPlot=F)
putative_imports <- which(pvals < alpha) + 1

layout(matrix(1:2, 1, 2))
plot(GlobalTree, label.offset=0.1, main="Globally sampled", ...)
tiplabels(tip=which(GlobalTree$tip.label %in% as.character(setdiff(1:length(GlobalTree$tip.label),1:pop1_sample_size))), pch=24)
plot(LocalTree, label.offset=0.1, main="Locally sampled", ...)
tiplabels(tip=putative_imports, pch=21, bg="green")
}
