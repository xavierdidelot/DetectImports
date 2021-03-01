
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DetectImports

The goal of DetectImports is to detect importations in a dated phylogeny
of local samples only, ie no samples from external locations (otherwise
it becomes more like a phylogeography type of inference). I would like
to make as few assumptions as possible about what is going on outside of
the local population. If there were only local transmission, we can
describe what phylogenies are expected to look like under a coalescent
model with varying population size. We can think about this as the null
distribution, and calculate for each sample the probability of
importation based on how unlikely the observation would be under a model
of local transmission only.

## Installation

You can install DetectImports from github with:

``` r
devtools::install_github("xavierdidelot/DetectImports")
```

The package can then be loaded using:

``` r
library(DetectImports)
```

## Example

Cf vignette.
