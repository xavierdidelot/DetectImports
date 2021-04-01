#' Compute dates of nodes, number of lineages through time, terminal branch lengths and coalescent intersection
#' @param phy Tree
#' @return Tree including a matrix with dates of nodes, number of lineages through time and coalescent intersection
#' @export
keyStats = function (phy) {
  n=Ntip(phy)
  fathers=NA
  fathers[phy$edge[,2]]=phy$edge[,1]
  rootdate=phy$root.time
  if (is.null(rootdate)) rootdate=0
  dates=rootdate+unname(dist.nodes(phy)[Ntip(phy)+1,])
  l=length(dates)
  lineages=rep(NA,l)
  for (i in 1:l) {
    t=dates[i]
    lineages[i]=length(which(dates[1:n]>=t))-length(which(dates[(n+1):l]>=t))
  }
  bralen=rep(NA,l)
  for (i in 1:l) {
    t1=dates[i]
    t2=dates[fathers[i]]
    bralen[i]=t1-t2
  }
  coal=rep(NA,l)
  r=rank(dates,ties.method='first')
  orderedDates=dates;orderedDates[r]=dates
  orderedLineages=lineages;orderedLineages[r]=lineages
  di=diff(orderedDates)
  for (i in 1:n) {
    t1=dates[i]
    t2=dates[fathers[i]]
    i2=which(orderedDates==t1)[1]
    i1=which(orderedDates==t2)[1]
    coal[i]=sum(di[i1:(i2-1)]*(orderedLineages[(i1+1):i2]-1))
  }
  coalint=c(coalIntervals(phy),rep(NA,phy$Nnode))
  mat=cbind(dates,lineages,bralen,coal,coalint)
  phy$stats=mat
  return(phy)
}

coalIntervals=function(phy)  {
  int=rep(NA,Ntip(phy))
  fathers=NA
  fathers[phy$edge[,2]]=phy$edge[,1]
  rootdate=phy$root.time
  if (is.null(rootdate)) rootdate=0
  dates=rootdate+unname(dist.nodes(phy)[Ntip(phy)+1,])
  fathers[Ntip(phy)+1]=length(fathers)+1
  fathers=c(fathers,0)
  dates=c(dates,-Inf)
  tab=cbind(dates,fathers)
  #tab[,1]=dates at bottom;tab[,2]=father
  tab[ ,1] <- max(tab[ ,1])-tab[ ,1]#convert dates to ages
  ex <- rep(0, nrow(tab))#Keep track of which nodes are active
  MySort <- sort(tab[ ,1],index.return = TRUE); ind <- MySort$ix
  cur <- 1#Start with oldest leaf
  while (tab[cur,2] > 0) {ex[cur] <- 1;cur <- tab[cur,2]}#Activate path to root
  ex[length(ex)] <- 1#Activate root
  for (l in 2:Ntip(phy)) {#For all leaves in increasing order of date
    isanc <- rep(0, nrow(tab))#Ancestors of the current node
    anc <- l
    while (tab[anc,2] > 0)  {
      isanc[anc] <- 1
      anc <- tab[anc,2]
    }
    bra1 <- 0
    start <- FALSE
    found <- FALSE
    curage <- 0
    k <- 0
    for (i in ind) {
      if (i == l)  {
        start <- TRUE
        curage <- tab[l,1]
        next
      }
      if (ex[i]==0) next #Ignore non-existent nodes
      if (start)  {
        if (!found)  {
          bra1 <- bra1 + k*(tab[i,1]-curage)
          if (isanc[i])  found <- TRUE
        }
      }
      curage <- tab[i,1]
      if (i<=Ntip(phy))  k <- k + 1 else k <- k-ex[i] + 1
    }
    int[l]=bra1
    if (k != 1) warning('k!=1')
    #Activate all ancestors of current node
    cur <- l
    while (TRUE) {
      if (ex[cur] == 1) {ex[cur] <- 2;break}
      ex[cur] <- 1
      cur <- tab[cur,2]
    }
  }
  if (min(ex)==0) warning('Warning: some nodes have not been activated')
  if (length(which(ex==2))!=(Ntip(phy)-1)) warning('Warning: some internal nodes have not been activated twice')
  return(int)
}
