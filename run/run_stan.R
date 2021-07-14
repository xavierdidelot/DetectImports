library(DetectImports)
library(posterior)
library(cmdstanr)
library(ggplot2)
compute_ci <- function(x, conf=0.95) {
  ci <- c()
  x_ord <- order(x)
  if(length(x)%%2==0) {
    l<-length(x)/2
    p1 <- x[x_ord][(l+1):length(x)]
    p2 <- x[x_ord][1:l]
  } else {
    l <-floor(length(x)/2)
    p1 <- x[x_ord][(l+2):length(x)]
    p2 <- x[x_ord][1:l]
  }
  ci[1] <- p2[l*(1-conf)]
  ci[2] <- p1[l*conf]
  return(ci)
}

ntip<-500
set.seed(123)
tree<-simImports(localPopStart=2020,importDates=c(2020.25,2020.5),
                samplingStartDate=2020,samplingEndDate=2021,samplingNumber=ntip)

real=tree$imports#real imports

if (is.null(tree$stats)) m<-keyStats(tree)$stats else m<-tree$stats

m<-cbind(m,"is_imp"=sapply(c(1:nrow(m)), function(i) i %in% real))

toana<-which(!is.na(m[,'coalint']))
dates<-m[toana,'dates']
coalints<-m[toana,'coalint']
imp<-m[toana, "is_imp"]
imp <- as.logical(imp)

print(dates)
print(coalints)

cnames <- sapply(c(1:length(coalints)),function(i) paste0(
  "f[",i,"]"))

mod <- cmdstan_model(paste0("../stan/gpmodel.stan"))
data_list <- list(N = length(coalints),intervals=coalints, T_s=dates, l=0.7, sigma=1)
fit <- mod$sample(
  data = data_list,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 500, 
  iter_sampling = 1e4
)

draws_array <- fit$draws()
draws_df <- as_draws_df(draws_array)

c_means <- draws_df[cnames]

ci_lo <- apply(c_means, 2, function(x) compute_ci(x, conf=0.98)[1])
ci_hi <- apply(c_means, 2, function(x) compute_ci(x, conf=0.98)[2])

data_df <- data.frame(x=dates, hi=ci_hi, lo=ci_lo, int=coalints, imp=imp)

p <- ggplot(data_df, aes(x)) + geom_ribbon(aes(ymin=lo, ymax=hi)) + geom_point(aes(y=int,color=imp))
plot(p)
dev.off()
