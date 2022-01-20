rm(list=ls())

library(ape)
library(treeio)
library(DetectImports)

countMigs=function(phy,lo)
{
  parents=c()
  parents[phy$edge[,2]]=phy$edge[,1]
  n=Ntip(phy)
  w=rep(NA,n)
  for (i in 1:n) {
    k=i
    while (!is.na(k) && lo[k]==0) k=parents[k]
    if (is.na(k)) k=0
    w[i]=k
  }
  return(length(unique(w))-1)
}


totrep=300
res=c()
for (rep in 1:totrep) {
  set.seed(rep)
  print(rep)

  #Generate XML file for Master
  system('cp start.xml master.xml')
  localcoalrate=1
  if (rep<=totrep/3) globalcoalrate=localcoalrate*0.1 else if (rep<=totrep*2/3) globalcoalrate=localcoalrate*0.25 else globalcoalrate=localcoalrate*0.5
  forwardmigrate=runif(1,min=0,max=0.5)#forward in time migration rate, symmetric between local and global demes
  forwardmigrate=forwardmigrate/(1/globalcoalrate)*(1/localcoalrate)#scaling so that backward in time migration from local to global is unif(0,0.5)
  migrate1=forwardmigrate*(1/globalcoalrate)/(1/localcoalrate)#backward in time migration from local to global
  migrate2=forwardmigrate*(1/localcoalrate)/(1/globalcoalrate)#backward in time migration from global to local
  system(sprintf('perl -i -wpe"s/localcoalrate/%f/g" master.xml',localcoalrate))
  system(sprintf('perl -i -wpe"s/globalcoalrate/%f/g" master.xml',globalcoalrate))
  system(sprintf('perl -i -wpe"s/migrate1/%f/g" master.xml',migrate1))
  system(sprintf('perl -i -wpe"s/migrate2/%f/g" master.xml',migrate2))

  dates=seq(from=0,to=1,length.out=500)
  locs=rep(0,length(dates))#c(rep(0,50),rep(1,50))
  f=file('master.xml',open='a')
  for (i in 1:length(dates)) {
    l=sprintf('\t\t\t<lineageSeedMultiple spec="MultipleIndividuals" copies="1" time="%.4f"><population spec="Population" type="@L" location="%d"/></lineageSeedMultiple>\n',dates[i],locs[i])
    writeLines(l,f)
  }
  close(f)
  system('cat end.xml >> master.xml')

  #Generate tree with Master, read tree, plot it and count number of migrations
  system(sprintf('beast2 -seed %d master.xml > /dev/null 2> /dev/null',rep))
  t=read.beast('simu_master.tree')
  phy=as.phylo(t)
  #plot(phy,show.tip.label=F)
  #axisPhylo(1)
  lo=as.numeric(t@data$location)
  lo[as.numeric(t@data$node)]=lo
  #mi=sum(lo[phy$edge[,1]]==1 & lo[phy$edge[,2]]==0)
  mi=countMigs(phy,lo)

  #Detect imports and compare with number of migrations
  phy2=collapse.singles(phy)
  set.seed(rep)
  r=detectImports(phy2,verbose = F)
  set.seed(rep)
  c=detectImports(phy2,verbose = F,constant=T)
  res=rbind(res,(c(mi,length(which(r$pvals<0.05)),length(which(r$pvals<0.01)),length(which(r$pvals<0.001)),
                      length(which(c$pvals<0.05)),length(which(c$pvals<0.01)),length(which(c$pvals<0.001)))))

}

save.image('res.RData')

pdf('figMaster.pdf',4,12)
par(mfrow=c(3,1),mar=c(5,5,2,2))
ind=1:(totrep/3)
plot(res[ind,1]+runif(length(ind)),res[ind,4]+runif(length(ind)),xlab='Correct number of migrations',ylab='Inferred number of imports',xlim=c(0,20),ylim=c(0,20),yaxs="i",xaxs="i")
lines(c(0,20),c(0,20))
r=lm(res[ind,4]~res[ind,1]-1);abline(r,lty='dashed')
ind=ind+totrep/3
plot(res[ind,1]+runif(length(ind)),res[ind,4]+runif(length(ind)),xlab='Correct number of migrations',ylab='Inferred number of imports',xlim=c(0,20),ylim=c(0,20),yaxs="i",xaxs="i")
lines(c(0,20),c(0,20))
r=lm(res[ind,4]~res[ind,1]-1);abline(r,lty='dashed')
ind=ind+totrep/3
plot(res[ind,1]+runif(length(ind)),res[ind,4]+runif(length(ind)),xlab='Correct number of migrations',ylab='Inferred number of imports',xlim=c(0,20),ylim=c(0,20),yaxs="i",xaxs="i")
lines(c(0,20),c(0,20))
r=lm(res[ind,4]~res[ind,1]-1);abline(r,lty='dashed')
par(xpd=NA)
text(-3,73,'A',cex=2)
text(-3,47,'B',cex=2)
text(-3,21,'C',cex=2)
dev.off()
system('open figMaster.pdf')
