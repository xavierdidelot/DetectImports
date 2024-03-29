#' Quick test for the presence of imports
#' @param tree Tree
#' @param constant Whether to assume that the local population size is constant
#' @param online Whether to perform the online test (ignored if contant=T)
#' @param epsilon Smoothing precision parameter (ignored if constant=T)
#' @param adjust Method for adjusting p-values (default is fdr)
#' @return Results of importation test
#' @export
detectImportsFAST=function(tree,constant=F,online=F,epsilon,adjust='fdr')
{
  if (is.null(tree$stats)) tree=keyStats(tree)
  m=tree$stats
  dates=m[1:Ntip(tree),'dates']
  coalints=m[1:Ntip(tree),'coalint']
  if (constant) {
    Ne=median(coalints,na.rm = T)/log(2)
    pvals=1-pexp(coalints,1/Ne)
    NeHat=rep(Ne,length(pvals))
  } else {
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
    pvals=1-pexp(coalints,1/NeHat)
  }
  pvals=p.adjust(pvals,adjust)
  makeOutput(tree,pvals,NeHat)
}

#' ESD test for the presence of imports
#' @param tree Tree
#' @param epsilon Smoothing precision parameter
#' @param alpha Threshold for p-value significance, default is 0.05
#' @param maxi Maximum number of outliers in dataset
#' @param showPlot Whether to show a plot of the test
#' @return Results of importation test
#' @export
detectImportsESD=function(tree,epsilon,alpha=0.05,maxi=round(Ntip(tree)/10),showPlot=F)
{
  if (is.null(tree$stats)) tree=keyStats(tree)
  m=tree$stats
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

  if (showPlot) {
    significant=rep(0,maxi)
    if (length(w)>0) significant[1:w]=1
    res = cbind(1:maxi,R,lam,toana[out],significant)
    colnames(res)=c("NumOutliers","TestStat", "CriticalVal","Outlier","Significant")
    print(res)
  }
  makeOutput(tree,pvals)
}

#' Bayesian test for the presence of imports
#' @param tree Tree
#' @param constant Whether to assume that the local population size is constant
#' @param adjust Method for adjusting p-values (default is none)
#' @param nchains Number of chains to run in parallel
#' @param verbose Whether to produce verbose output
#' @param diagnostics Whether to produce diagnostic warnings
#' @param iter Number of iterations to run, with first quarter discarded
#' @param seed Seed
#' @param extrapara Named list to pass additional parameters, such as prior hyperparameters
#' @return Results of importation test
#' @export
detectImports=function(tree,constant=FALSE,adjust='none',verbose=F,diagnostics=F,nchains=4,iter=4000,seed=NULL,extrapara=list())
{
  if (is.null(tree$stats)) tree=keyStats(tree)
  m=tree$stats

  toana<-which(!is.na(m[,'coalint']))
  dates<-m[toana,'dates']
  coalints<-m[toana,'coalint']

  if (constant) {
    mod <- cmdstan_model(system.file('stan','constant.stan',package='DetectImports',mustWork = T),compile=F)
    data_list <- list(N = length(coalints), intervals=coalints, T_s=dates)
    coalnames <- rep("coal_mean",length(coalints))
  }
  else
  {
    mod <- cmdstan_model(system.file('stan','gpmodel.stan',package='DetectImports',mustWork = T),compile=F)
    if (!hasName(extrapara,'M')) extrapara$M=20
    if (!hasName(extrapara,'c')) extrapara$c=2
    if (!hasName(extrapara,'al')) extrapara$al=5
    if (!hasName(extrapara,'bl')) extrapara$bl=5
    if (!hasName(extrapara,'sigmaalpha')) extrapara$sigmaalpha=5
    data_list <- list(N = length(coalints), intervals=coalints, T_s=dates, M=extrapara$M, c=extrapara$c, al=extrapara$al, bl=extrapara$bl, sigmaalpha=extrapara$sigmaalpha)
    coalnames <- sapply(c(1:length(coalints)),function(i) paste0("coal_means[",i,"]"))
  }

  if (verbose) {
    mod$compile()
    fit <- mod$sample(data = data_list,adapt_delta=0.9,chains = nchains,parallel_chains = nchains,refresh = round(iter*0.1),iter_sampling = round(iter*0.75),iter_warmup = round(iter*0.25),seed=seed)
  } else {
    invisible(capture.output(suppressMessages({
      mod$compile()
      fit <- mod$sample(data = data_list,adapt_delta=0.9,chains = nchains,parallel_chains = nchains,refresh = 0,show_messages=F,iter_sampling = round(iter*0.75),iter_warmup = round(iter*0.25),seed=seed)
    })))
  }

  fit_summary=NULL;sampler_diagnostics=NULL
  if (diagnostics) {
    if (!constant) {
      fit_summary <- fit$summary(c("alpha", "l", "coal_means", "f_tilde"))
    } else {
      fit_summary <- fit$summary("coal_mean")
    }

    sampler_diagnostics <- fit$sampler_diagnostics()
    sampler_diagnostics <- as_draws_df(sampler_diagnostics)

    if(!all(fit_summary$rhat < 1.05)) {
      warning("Poor convergence detected: unsatisfactory rhat")
    }
    if (!all(fit_summary$ess_bulk > 2000)) {
      warning("Poor convergence detected: unsatisfactory ESS")
    }
    if(length(which(sampler_diagnostics$divergent__ > 1e-8))>1) {
      warning(sprintf("%d/%d divergent transitions detected",length(which(sampler_diagnostics$divergent__ > 1e-8)),length(sampler_diagnostics$divergent__)))
    }
  }

  draws_array <- fit$draws()
  draws_df <- as_draws_df(draws_array)

  coal_m <- suppressWarnings(draws_df[coalnames])
  pv = rep(NA,length(coalints))
  mat=as.matrix(coal_m)
  for (i in 1:length(pv)) pv[i]=mean(1-pexp(coalints[i],mat[,i]))
  g_med <- apply(coal_m, 2, median)
  pvals=rep(1,Ntip(tree))
  #pvals[toana]=1-pexp(coalints,g_med)
  pvals[toana]=pv
  mus=mus_low=mus_high=rep(NA,Ntip(tree))
  mus[toana]=1/g_med
  mus_high[toana]=1/apply(coal_m, 2, function (x) {return(quantile(x,probs=0.025))})
  mus_low [toana]=1/apply(coal_m, 2, function (x) {return(quantile(x,probs=0.975))})

  pvals=p.adjust(pvals,adjust)
  res=makeOutput(tree,pvals,mus, fit_summary, sampler_diagnostics)
  res$mus_high=mus_high
  res$mus_low =mus_low

  return(res)
}

makeOutput=function(tree,pvals,mus=NULL, fit_summary=NULL, diagnostics=NULL)
{
  res=list()
  class(res)<-'resDetectImports'
  res$tree=tree
  res$pvals=pvals
  res$mus=mus
  res$fit_summary=fit_summary
  res$diagnostics=diagnostics
  return(res)
}
