#' Simulation of coalescent dated phylogeny given a fluctuating Ne(t) function
#' @param dates Sampling dates
#' @param NeFun Population size function Ne(t)
#' @param NeMin Minimum value of Ne(t)
#' @return A simulated dated phylogeny
#' @export
simCoal = function(dates=1990:2010,NeFun=function(x){return(10)},NeMin=NA) {
  if (is.na(NeMin)) NeMin=optimize(NeFun,c(-1e5,max(dates)))$objective
  s <- sort(dates,decreasing=TRUE,index.return = TRUE)
  tim <- s$x
  ind <- s$ix
  n <- length(tim)
  nodes <- cbind(-Inf,ind[1],-Inf)#Start with one node at time -Inf and with the most recent isolate connected to it
  i <- 2
  while (i <= n) {#Graft branches one by one, from most recent to least recent
    curt <- tim[i]#Current time:start with date of isolate and go back in time until coalescence happens
    accept=F
    while (accept==F) {
      r <- -log(runif(1)) * NeMin
      fi <- which( nodes[ ,1] < curt )[1]
      if (fi<=nrow(nodes)) for (j in (fi:nrow(nodes)))  {
        if (r > (curt-nodes[j,1]) * (i-j))  {
          r <- r-(curt-nodes[j,1]) * (i-j)
          curt <- nodes[j,1]
        } else {
          curt <- curt-r/(i-j)#Proposed coalescent time
          break
        }
      }
      #Filtering the Poisson process to obtain non-homogeneous Poisson
      if (runif(1)<NeMin/NeFun(curt)) accept=T
    }
    #Create new node
    a <- nodes[ ,2:3]
    a[a >= j + n] <- a[a >= j + n] + 1
    nodes[ ,2:3] <- a#Renumbering according to table insertion in next line
    nodes2=c(curt,ind[i],0)
    if (1<=j-1) nodes2=rbind(nodes[1:(j-1), ],nodes2)
    if (j<=nrow(nodes)) nodes2=rbind(nodes2,nodes[j:nrow(nodes),])
    nodes=unname(nodes2)
    #Now choose on which branch to coalesce among the branches alive at date curt
    no <- j
    side <- 2
    w <- 1 + floor(runif(1) * (nrow(nodes)-j))
    while (w > 0)  {
      no <- no + side-1
      side <- 3-side
      if (nodes[no,side + 1] <= n ||(nodes[no,side + 1] > n && nodes[nodes[no,side + 1]-n,1] > curt))  {
        w <- w-1
      }
    }
    nodes[j,3] <- nodes[no,side + 1]
    nodes[no,side + 1] <- n + j
    i <- i + 1
  }
  v=nrow(nodes)-1
  if (nrow(nodes)>2) v=c(v,1:(nrow(nodes)-2))
  nodes=nodes[v,,drop=F]
  m=nodes[,2:3]
  m[which(m>n)]=m[which(m>n)]+1
  nodes[,2:3]=m
  nodes <- rbind(matrix(0, nrow = n, ncol = 3),nodes)
  nodes[1:n,1] <- dates

  #Convert into phylo object from package ape
  t=list()
  t$Nnode=n-1
  t$tip.label=as.character(1:n)
  t$edge=matrix(NA,2*n-2,2)
  t$edge.length=rep(NA,n*2-2)
  t$root.time=nodes[n+1,1]
  c=1
  for (i in (n+1):nrow(nodes)) for (j in 2:3) {
    t$edge[c,1]=i
    t$edge[c,2]=nodes[i,j]
    t$edge.length[c]=nodes[nodes[i,j],1]-nodes[i,1]
    c=c+1
  }
  class(t)='phylo'
  return(t)
}

#' Plot a dated phylogeny and Ne function
#' @param tree Dated phylogeny
#' @param NeFun Population size function Ne(t)
#' @return Figure
#' @export
plotBoth = function(tree,NeFun) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow=c(2,1),mar=c(4,4,1,4))
  if (!is.null(tree$root.time)) from=tree$root.time else from=-max(dist.nodes(tree)[Ntip(tree)+1,])
  to=from+max(dist.nodes(tree)[Ntip(tree)+1,])
  xs=seq(from=from,to=to,length.out=100)
  ys=xs
  for (i in 1:length(ys)) ys[i]=NeFun(ys[i])
  plot(xs,ys,type='l',xlab='',ylab='Population size', bty='l',ylim=c(0,1.05*max(ys)))
  plot(tree,show.tip.label = F)
  axisPhylo(1,backward = F)
}

#' Simulation with imports
#' @return A simulated dated phylogeny
#' @export
simImports = function() {

  NeFunLinear=function(t,texp,k) {pmax(0,(t-texp)*k)}

  dateOrigin=2020.25
  NeFun= function(t) {NeFunLinear(t,2020,1)}
  NeFun2=function(t) {NeFunLinear(t,2020.75,1)}

  days=seq(2020,2021.5,1/365)
  dates =sample(days,900,replace=T,prob=NeFun (days))
  dates=sort(dates+runif(length(dates))/365)
  dates2=sample(days,100,replace=T,prob=NeFun2(days))
  dates2=sort(dates2+runif(length(dates2))/365)
  t2=simCoal(dates2,NeFun2,1e-2)
  if (t2$root.time<dateOrigin) stop('here')
  toAdd=t2$root.time-dateOrigin
  t1=simCoal(c(dates,dateOrigin),NeFun,1e-2)

  w=which(t1$edge[,2]==length(dates)+1)
  t1$edge.length[w]=t1$edge.length[w]+toAdd
  t2$tip.label=as.numeric(t2$tip.label)+900
  t=bind.tree(t1,t2,where=length(dates)+1,position=0)
  t$imports=901
  return(t)
}
