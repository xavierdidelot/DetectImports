#' Test for the presence of imports assuming that the local population has a constant population size
#' @param tree Tree
#' @param Ne Population size (default is estimated from data)
#' @param adjust Method for adjusting p-values (default is fdr)
#' @param showPlot Whether to show a plot of the test
#' @return p-values for importation
#' @export
test0=function(tree,Ne,adjust='fdr',showPlot=T)
{
  if (is.null(tree$stats)) m=keyStats(tree)$stats else m=tree$stats
  coalints=m[1:Ntip(tree),'coalint']
  if (missing(Ne)) Ne=median(coalints,na.rm = T)/log(2)
  if (showPlot) {
    h=hist(coalints,main='',breaks=20,xlab='Coalescent intervals',ylab='Frequency')
    br=h$breaks
    lines(br[-1]-diff(br)[1]/2,length(coalints)*(pexp(br[-1],1/Ne)-pexp(br[-length(br)],1/Ne)))
  }
  pvals=1-pexp(coalints,1/Ne)
  pvals=p.adjust(pvals,adjust)
  message(sprintf('%d imports were found with p<0.01. Lowest p-value was %.2e',length(which(pvals<0.01)),min(pvals,na.rm=T)))
  return(pvals)
}

#' Semi-parametric test for the presence of imports
#' @param tree Tree
#' @param epsilon Smoothing precision parameter
#' @param adjust Method for adjusting p-values (default is fdr)
#' @param online Whether to perform the online test
#' @param showPlot Whether to show a plot of the test
#' @return p-values for importation
#' @export
test1=function(tree,epsilon,adjust='fdr',online=F,showPlot=T)
{
  if (is.null(tree$stats)) m=keyStats(tree)$stats else m=tree$stats
  dates=m[1:Ntip(tree),'dates']
  coalints=m[1:Ntip(tree),'coalint']
  if (missing(epsilon)) epsilon=(max(dates)-min(dates))/20
  l=length(dates)
  NeHat=rep(NA,l)
  for (i in 1:l) {
    if (online) preexisting=length(which(dates[i]-dates>0))
    k=0
    w=c()
    while (length(w)<5) {#increase window size in case there are not enough leaves within epsilon
      k=k+1
      if (online) {
        w=which(dates[i]-dates>0 & dates[i]-dates<epsilon*k)
        if (length(w)==preexisting) break
      } else {
        w=setdiff(which(abs(dates-dates[i])<epsilon*k),i)
        if (length(w)==l-1) break
      }
    }
    if (length(w)<3) w=c()
    NeHat[i]=median(coalints[w],na.rm=T)/log(2)
  }
  if (showPlot) {
    plot(dates,coalints,xlab='',ylab='')
    ix=sort(dates,index.return=T)$ix
    lines(dates[ix],NeHat[ix],col='red')
  }
  pvals=1-pexp(coalints,1/NeHat)
  pvals=p.adjust(pvals,adjust)
  message(sprintf('%d imports were found with p<0.05. Lowest p-value was %.2e',length(which(pvals<0.05)),min(pvals,na.rm=T)))
  return(pvals)
}

#' ESD test for the presence of imports
#' @param tree Tree
#' @param epsilon Smoothing precision parameter
#' @param alpha Threshold for p-value significance, default is 0.05
#' @param maxi Maximum number of outliers in dataset
#' @param showPlot Whether to show a plot of the test
#' @return p-values for importation
#' @export
test2=function(tree,epsilon,alpha=0.05,maxi=round(Ntip(tree)/10),showPlot=F)
{
  if (is.null(tree$stats)) m=keyStats(tree)$stats else m=tree$stats
  toana=which(!is.na(m[,'coalint']))
  dates=m[toana,'dates']
  coalints=m[toana,'coalint']
  if (missing(epsilon)) epsilon=(max(dates)-min(dates))/20

  l=length(dates)
  NeHat=rep(NA,l)
  for (i in 1:l) {
    k=0
    w=c()
    while (length(w)<5) {#just in case there are not enough leaves within epsilon
      k=k+1
      w=setdiff(which(abs(dates-dates[i])<epsilon*k),i)
    }
    NeHat[i]=median(coalints[w],na.rm=T)/log(2)
  }

  #Ne=median(coalints,na.rm = T)/log(2)
  y=qnorm(pexp(coalints,1/NeHat))#transform Exp(1/NeHat(t)) into N(0,1)

  if (showPlot) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow=c(2,1),mar=c(4,4,1,4))
    plot(dates,y)
    qqnorm(y)
  }

  n = length(y)
  lam=R=out=rep(NA,maxi)
  inds=1:n
  for (i in 1:maxi){
    #Compute test statistics
    v=y[inds]
    v2 = abs(v-mean(v))/sd(v)
    R[i]=max(v2)
    out[i]=inds[which(v2==R[i])[1]]
    inds=setdiff(inds,out[i])
    #Compute critical values
    p = 1 - alpha/(2*(n-i+1))
    t = qt(p,(n-i-1))
    lam[i] = t*(n-i) / sqrt((n-i-1+t**2)*(n-i+1))
  }

  w=which(R>lam)
  pvals=rep(1,Ntip(tree))
  if (length(w)>0) {
    w=w[length(w)]
    pvals[toana[out[1:w]]]=alpha
  }
  message(sprintf('%d imports were found with p<%f.',length(which(pvals<=alpha)),alpha))

  if (showPlot) {
    significant=rep(0,maxi)
    if (length(w)>0) significant[1:w]=1
    res = cbind(1:maxi,R,lam,toana[out],significant)
    colnames(res)=c("NumOutliers","TestStat", "CriticalVal","Outlier","Significant")
    print(res)
  }

  return(pvals)
}

#' Bayesian test
#' @param tree Tree
#' @param constant Whether to assume that the local population size is constant
#' @param adjust Method for adjusting p-values (default is fdr)
#' @return p-values for importation
#' @export
testBayes=function(tree,constant=FALSE,adjust='fdr')
{
  if (is.null(tree$stats)) m<-keyStats(tree)$stats else m<-tree$stats

  toana<-which(!is.na(m[,'coalint']))
  dates<-m[toana,'dates']
  coalints<-m[toana,'coalint']

  if (constant)
    mod <- cmdstan_model(file.path(find.package('DetectImports'),'stan','constant.stan'))
  else
    mod <- cmdstan_model(file.path(find.package('DetectImports'),'stan','gpmodel.stan'))
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

  coalnames <- sapply(c(1:length(coalints)),function(i) paste0("coal_means[",i,"]"))
  coal_m <- suppressWarnings(draws_df[coalnames])
  pv = rep(NA,length(coalints))
  mat=as.matrix(coal_m)
  for (i in 1:length(pv)) pv[i]=mean(1-pexp(coalints[i],mat[,i]))
  #g_med <- apply(coal_m, 2, median)
  pvals=rep(1,Ntip(tree))
  #pvals[toana]=1-pexp(coalints,g_med)
  pvals[toana]=pv

  pvals=p.adjust(pvals,adjust)
  message(sprintf('%d imports were found with p<0.05. Lowest p-value was %.2e',length(which(pvals<0.05)),min(pvals,na.rm=T)))
  return(pvals)
}
