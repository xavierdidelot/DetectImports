library(ape)
rm(list=ls())
set.seed(3)
t=rtree(10)
t$root.time=2019
dates=t$root.time+unname(dist.nodes(t)[Ntip(t)+1,])
t$tip.label=sprintf('s%d',rank(dates[1:10]))
o=rep(NA,9)
o[c(3,2,1,6,7,4,8,9,5)]=1:9
t$node.label=sprintf('c%d',o)
pdf('/tmp/fig1.pdf',7,7)
plot(t,show.tip.label = F,show.node.label = F,font=1,label.offset = 0.01,cex=1.5)
par(xpd=NA)
vpos = get("last_plot.phylo", envir = .PlotPhyloEnv)$yy
lpos = get("last_plot.phylo", envir = .PlotPhyloEnv)$xx
text(lpos[1:10]+0.05,vpos[1:10],substitute(paste(italic('s'))),cex=1.3)
text(lpos[1:10]+0.08,vpos[1:10]-0.1,substr(t$tip.label,2,10),cex=0.8,adj=0)
text(lpos[11:19]+0.05,vpos[11:19],substitute(paste(italic('c'))),cex=1.3)
text(lpos[11:19]+0.08,vpos[11:19]-0.1,substr(t$node.label,2,10),cex=0.8,adj=0)

axis(1,at=c(0,1,2),labels=c(2019,2020,2021),cex=1.5)
w=which(t$edge[,2]==which(t$tip.label=='s10'))
dates=dates-t$root.time
d=dates[t$edge[w,1]]
fathers=NA
fathers[t$edge[,2]]=t$edge[,1]
w=which(t$tip.label=='s7');lines(c(max(d,dates[fathers[w]]),dates[w]),c(3,3),col='red',lwd=2)
w=which(t$tip.label=='s8');lines(c(max(d,dates[fathers[w]]),dates[w]),c(7,7),col='red',lwd=2)
w=which(t$tip.label=='s6');lines(c(max(d,dates[fathers[w]]),dates[w]),c(6,6),col='red',lwd=2)
w=which(t$tip.label=='s9');lines(c(max(d,dates[fathers[w]]),dates[w]),c(9,9),col='red',lwd=2)
w=which(t$tip.label=='s5');lines(c(max(d,dates[fathers[w]]),dates[w]),c(8,8),col='red',lwd=2)
w=which(t$tip.label=='s4');lines(c(max(d,dates[fathers[w]]),dates[w]),c(10,10),col='red',lwd=2)
w=which(t$node.label=='c7');lines(c(max(d,dates[fathers[w]]),dates[10+w]),c(6.5,6.5),col='red',lwd=2)
dev.off()
system('open /tmp/fig1.pdf')
