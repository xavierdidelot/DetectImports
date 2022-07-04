library(DetectImports)
library(ape)
set.seed(0)
mod <- cmdstanr::cmdstan_model(system.file('stan','gpmodel.stan',package='DetectImports',mustWork = T),compile=F)
mod$compile()
Ne=function(t){if (t>2020.5) return(1) else return(0.2)}
tree=simCoal(seq(2020,2021,length.out = 100),Ne)
start_time <- Sys.time()
res2=detectImports(tree,seed=0)
end_time <- Sys.time()
print(end_time - start_time)
sum(res2$pvals)#49.32972

