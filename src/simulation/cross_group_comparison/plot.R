#!/usr/bin/env Rscript
library(phyloseq)
library(LTN)
library(ape)
argv=commandArgs(TRUE)
WORK_DIR=argv[1]
K=50
png(paste0(WORK_DIR,'/results/cross_group_comparison/ROC.png'),height = 500,width = 1000)
par(mfrow=c(1,2),pty='s')
colPal=c(2,3,1,6)
for (i in 1:2){
  s=c('single_otu','multi_otu')[i]
  resultdir=paste0(WORK_DIR,'/results/cross_group_comparison/',s)
  roclist=sapply(c('r2.','r5.','sparse_lambda10','diagonal'),function(x){readRDS(grep(x,list.files(resultdir,full.names = T),value = T))})
  for (j in 1:4){
    plot(roclist[[j]],add=(j!=1),col=colPal[j],main=paste0('scenario ',i))
  }
  legend('bottomright',
         legend = c('DirFactor (r=2)','DirFactor (r=5)','LTN (lambda=10)','LTN (diagonal)'),
         fill=colPal)
}
dev.off()

tree=phy_tree(readRDS(paste0(WORK_DIR,'/cache/ps_otu50.RData')))
png(paste0(WORK_DIR,'/results/cross_group_comparison/PMAP_simulation_single.png'),height = 600,width = 700)
i=29
PMAP=readRDS(paste0(WORK_DIR,'/cache/cross_group_comparison/single_otu/LTN/pmap/i',i,'H1_sparse_lambda10.RData'))
sim0=readRDS(paste0(WORK_DIR,"/cache/cross_group_comparison/single_otu/sim",i,"H0.RData"))
sim1=readRDS(paste0(WORK_DIR,"/cache/cross_group_comparison/single_otu/sim",i,"H1.RData"))
otus=which(apply(sim0$cnt[sim0$group2,]!=sim1$cnt[sim1$group2,],2,sum)>0)
PMAP=readRDS(paste0(WORK_DIR,'/cache/cross_group_comparison/single_otu/LTN/pmap/i',i,'H1_sparse_lambda10.RData'))
tips=rep('',K)
tips[otus]='OTU1'
plot_pmap(PMAP,tree,' ',tip_label = tips)
dev.off()
png(paste0(WORK_DIR,'/results/cross_group_comparison/PMAP_simulation_multi.png'),height = 600,width = 700)
i=22
sim0=readRDS(paste0(WORK_DIR,"/cache/cross_group_comparison/multi_otu/sim",i,"H0.RData"))
sim1=readRDS(paste0(WORK_DIR,"/cache/cross_group_comparison/multi_otu/sim",i,"H1.RData"))
otus=which(apply(sim0$cnt[sim0$group2,]!=sim1$cnt[sim1$group2,],2,sum)>0)
PMAP=readRDS(paste0(WORK_DIR,'/cache/cross_group_comparison/multi_otu/LTN/pmap/i',i,'H1_sparse_lambda10.RData'))
tips=rep('',K)
tips[otus]=paste0('OTU',1:8)
plot_pmap(PMAP,tree,' ',tip_label=tips,label_nodes = c(3,38))
dev.off()

