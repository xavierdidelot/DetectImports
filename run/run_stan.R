library(DetectImports)
library(posterior)
library(cmdstanr)
library(ggplot2)
library(evd)
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

compute_quantile <- function(x, quant=0.95) {
  x_ord <- order(x)
  return(x[x_ord][1-quant])
}

ntip<-1000
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

cnames <- sapply(c(1:length(coalints)),function(i) paste0(
  "f[",i,"]"))

coalnames <- sapply(c(1:length(coalints)),function(i) paste0(
  "coal_means[",i,"]"))

enames <- sapply(c(1:length(coalints)),function(i) paste0(
  "e_tilde[",i,"]"))


mod <- cmdstan_model(paste0("../stan/gpmodel.stan"))
data_list <- list(N = length(coalints), intervals=coalints, T_s=dates, shape=5, scale=5, M=30, c=2.0)
fit <- mod$sample(
  data = data_list,
  seed = 123,
  adapt_delta=0.9,
  chains = 4,
  parallel_chains = 4,
  refresh = 500, 
  iter_sampling = 3e3,
  iter_warmup = 1e3
)

draws_array <- fit$draws()
draws_df <- as_draws_df(draws_array)

#probs <- draws_df[cnames]
#e_probs <- apply(probs,2,mean)

c_means <- draws_df[cnames]
coal_m <- draws_df[coalnames]

ci_hi <- apply(c_means, 2, function(x) compute_quantile(x, quant=0.99))
g_med <- apply(coal_m, 2, median)

data_df <- data.frame(x=dates, hi=ci_hi, int=coalints, imp=imp, g=g_med)

e_draws <- draws_df[enames]
exp_e <- apply(e_draws,2,median)

edf <- data.frame(x=dates, y=exp_e)


#data_df <- data.frame(x=dates, prob=e_probs, int=coalints, imp=imp)

p <- ggplot(data_df, aes(x)) +  geom_ribbon(aes(ymin=0, ymax=ci_hi)) + geom_point(aes(y=int,shape=imp)) + geom_line(aes(y=g))
plot(p)
p <- ggplot(edf, aes(x)) + geom_point(aes(y=exp_e,shape=imp)) + geom_hline(yintercept=qexp(0.99,1))
plot(p)
dev.off()
