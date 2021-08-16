#!/usr/bin/env Rscript

library(LTN)
library(phyloseq)
library(data.tree)
library(ape)
library(ggplot2)
library(reshape2)


argv=commandArgs(TRUE)
WORK_DIR=argv[1]
resdir=paste0(WORK_DIR,"/results/application/")
system(paste0('mkdir -p ',resdir,'figures'))
system(paste0('mkdir -p ',resdir,'tables'))
lambda=10
K=50
ps=readRDS(paste0(WORK_DIR,"/cache/ps_otu50.RData"))
tree=phy_tree(ps)
tree$node.label=1:(K-1)
covariates=c("Case_Control","post_seroconversion","BF","Solid_Food","Eggs","Fish","Soy_Prod","Rye","Barley","Buckwheat_Millet")
covariates_name=c("Case_Control","post_seroconversion","Breastfeeding","Solid_Food","Eggs","Fish","Soy_Product","Rye","Barley","Buckwheat_Millet")
# PMAP plots
for (i in seq_along(covariates)){
  covariate=covariates[i]
  filenam=paste0(covariate,'_lambda',lambda)
  pmap=readRDS(paste0(resdir,'/pmap/',filenam,'.RData'))
  alpha=readRDS(paste0(resdir,'/alpha/',filenam,'.RData'))
  png(paste0(resdir,'figures/pmap_',filenam,'.png'),width = 1000,height = 600)
  if (covariate=='BF'){
    plot_pmap(pmap=pmap,tree=tree,main.text = covariates_name[i],alpha=alpha,label_nodes = which(pmap>0.5))
  }
  else{
    if (covariate=='post_seroconversion'){
      plot_pmap(pmap=pmap,tree=tree,main.text = covariates_name[i],alpha=alpha,label_nodes = c(38,39,44,45,26))
    }else{plot_pmap(pmap=pmap,tree=tree,main.text = covariates_name[i],alpha=alpha)}
  }
  dev.off()
}

# labels (Fig. S1)
png(paste0(resdir,'figures/node_label.png'),width = 1000,height = 600)
plot_pmap(pmap=rep(0,K-1),tree=tree,main.text = 'Node labels', label = 1:(K-1))
dev.off()

# tables
# PJAPs of dietary factors
pjapmat=matrix(do.call(rbind,lapply(3:length(covariates),function(i){
  covariate=covariates[i]
  filenam=paste0(covariate,'_lambda',lambda)
  pjap=readRDS(paste0(resdir,'/pjap/',filenam,'.RData'))
  num_on=sum(sample_data(ps)[,covariate]=='true')
  num_off=nrow(sample_data(ps))-num_on
  return(c(num_on,num_off,pjap))
})),ncol=3,dimnames = list(covariates_name[-(1:2)],c('on','off','pjap')))
pjapmat=pjapmat[sort(rownames(pjapmat)),]
write.csv(pjapmat,paste0(resdir,'tables/diet_pjap.csv'))

# node tables
tree_copy=tree
tree_copy$tip.label=1:K
tree_copy$node.label=(K+1):(2*K-1)
tn=as.Node(tree_copy)
tn$Do(function(x) {
  leftchild=names(x$children[1])
  rightchild=names(x$children[2])
  x$left <-leftchild
  x$right <-rightchild
}, traversal = "post-order",filterFun = isNotLeaf)
left=tn$Get('left', traversal = "pre-order", filterFun=isNotLeaf)
right=tn$Get('right', traversal = "pre-order", filterFun=isNotLeaf)
nodedf=data.frame(node=(K+1):(2*K-1),left=left,right=right,row.names = NULL)
nodepaths=nodepath(tree_copy)
taxtab=tax_table(ps)

threshold_node=function(pmap,threshold){
  return(which(pmap>threshold)+50)
}
descendant_otu=function(node){
  leftchild=nodedf[nodedf$node==node,'left']
  rightchild=nodedf[nodedf$node==node,'right']
  leftotu=unlist(lapply(nodepaths[which(unlist(lapply(nodepaths, function(x){leftchild%in%x})))],function(x){x[length(x)]}))
  rightotu=unlist(lapply(nodepaths[which(unlist(lapply(nodepaths, function(x){rightchild%in%x})))],function(x){x[length(x)]}))
  return(list(leftotu=leftotu,rightotu=rightotu))
}
common_diff=function(taxtableft,taxtabright,taxlevel=T,showpath=F){
  taxtabmerge=rbind(taxtableft,taxtabright)
  level_common=max(which(apply(rbind(taxtableft,taxtabright), 2, function(x){length(unique(x))})==1))
  if (showpath){
    taxa_common=paste0(taxtableft[1,1:level_common],collapse = '->')
  }else{
    names_split=sapply(taxtableft[1,1:level_common], function(x){strsplit(x,'__')[[1]][2]})
    taxa_common=names(names_split[max(which(!is.na(names_split)))])
  }
  if (!taxlevel){
    taxa_common=sapply(taxa_common,function(x){strsplit(x,'__')[[1]][2]})
  }
  if (level_common<ncol(taxtab)){
    leftvec=as.character(unique(taxtableft[,level_common+1]))
    rightvec=as.character(unique(taxtabright[,level_common+1]))
    if (!taxlevel){
      leftvec=sapply(leftvec, function(x){strsplit(x,'__')[[1]][2]})
      rightvec=sapply(rightvec, function(x){strsplit(x,'__')[[1]][2]})
      leftvec[is.na(leftvec)]='unclassified'
      rightvec[is.na(rightvec)]='unclassified'
    }
    lefttaxa=paste0(leftvec,collapse = ', ')
    righttaxa=paste0(rightvec,collapse=', ')
  }else{
    lefttaxa=paste0(rownames(taxtableft),collapse = ', ')
    righttaxa=paste0(rownames(taxtabright),collapse=', ')
  }
  return(list(common=taxa_common,left=lefttaxa,right=righttaxa))
}

for (i in seq_along(covariates)){
  covariate=covariates[i]
  filenam=paste0(covariate,'_lambda',lambda)
  pmap=readRDS(paste0(resdir,'/pmap/',filenam,'.RData'))
  alpha=readRDS(paste0(resdir,'/alpha/',filenam,'.RData'))
  nodes=threshold_node(pmap,0.5)
  otulist=lapply(nodes,function(x){descendant_otu(x)})
  taxtablist=lapply(otulist,function(x){list(taxtab[x$leftotu,],taxtab[x$rightotu,])})
  nodetaxa=lapply(taxtablist,function(x){common_diff(x[[1]],x[[2]],taxlevel=F)})
  common=unlist(lapply(nodetaxa, function(x){x[[1]]}))
  left=unlist(lapply(nodetaxa, function(x){x[[2]]}))
  right=unlist(lapply(nodetaxa, function(x){x[[3]]}))
  if (length(nodes)>0){
    write.csv(data.frame(node=nodes-50,taxon_in_common=common,taxa_on_left=left,taxa_on_right=right,alpha_hat=alpha[nodes-50]),paste0(resdir,'/tables/nodes_',filenam,'.csv'))
  }
}


# other figures
# MDS
ord=ordinate(ps,method='MDS',distance = 'bray')
png(paste0(resdir,'figures/MDS.png'),width=1000)
plot_ordination(ps,ord,color='Age_at_Collection',shape='Country')+scale_shape_manual(values=c(1,4))+facet_wrap(~Case_Control)+scale_color_gradientn(colours=c('red','green','blue'))+theme_bw()
dev.off()
# most abundant phylum
data=data.frame(sample_data(ps))
ps_phylum=tax_glom(ps,taxrank = 'Phylum')
data$most_abundant_phylum=apply(otu_table(ps_phylum),1,function(x){gsub('p__','',tax_table(ps_phylum)[which(x==max(x)),'Phylum'])})
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
data$post_seroconversion=data$Age_at_Collection>data$sero_age
png(paste0(resdir,'figures/most_abundant_phylum.png'))
plot_timepoints(data,x='Age_at_Collection',y='Subject_ID',color='most_abundant_phylum',shape='post_seroconversion',facet_by = 'Country',shape_values = c(4,1),color_values = c('black','orange','blue','green','red'),xlab='age at collection (day)',ylab='subject ID',color_labs = 'most abundant phylum',shape_labs = 'post-seroconversion')
dev.off()
# boxplot of empirical log-odds (seroconversion)
yyl=readRDS(paste0(WORK_DIR,"/cache/yyl_otu50.RData"))
THETA=data.frame(yyl$YL/yyl$Y)
THETA$age=data$Age_at_Collection
THETA$sero_age=data$sero_age
THETA$post_seroconversion=data$post_seroconversion
THETA$gid=rownames(THETA)
THETA$subid=data$Subject_ID
THETA$casecontrol=data$Case_Control
df=melt(THETA,id.vars=c('gid','post_seroconversion','age','subid','casecontrol','sero_age'))
df$psi=qlogis(df$value)
df$age_month=df$age/30
df$age_month_bin=cut(df$age_month,seq(0,42,3))
df$age_bin=cut(df$age,seq(0,max(df$age)+90,90))
df$node=paste0('node ',as.numeric(as.character(gsub('X','',df$variable)))-50)
df_subset=df[df$node %in% c('node 26', 'node 38', 'node 39', 'node 44','node 45') & df$age_month_bin!='(39,42]', ]
df_subset$node_f=factor(df_subset$node,levels=c('node 38', 'node 39', 'node 44','node 45','node 26'))
png(paste0(resdir,'figures/boxplot_empirical_logodds.png'),width = 1000)
ggplot(df_subset) +
  geom_boxplot(aes(y = psi, x = age_month_bin, fill = post_seroconversion))  +
  ylab('empirical log-odds') + xlab('age at collection (month)') + facet_wrap( ~ node_f, labeller = labeller(
    node_f =
      c(
        "node 38" = "A1 (Bacteroides)",
        "node 39" = "A2 (Bacteroides)",
        "node 44" = "A3 (Bacteroides)",
        "node 45" = "A4 (Bacteroides)",
        "node 26" = "A5 (Erysipelotrichaceae)"
      )
  )) +theme_bw()+ theme(axis.text.x = element_text(angle = 90))+labs(fill='post-seroconversion')+scale_fill_manual(values=c('orange','blue'))
dev.off()
# clostridiales
ps_order=tax_glom(ps,taxrank = 'Order')
data$o__Clostridiales_absolute=as.vector(as.matrix(otu_table(ps_order)[,'265871']))
data$o__Clostridiales_relative=as.vector(as.matrix(t(apply(otu_table(ps_order),1,function(x){x/sum(x)}))[,'265871']))
data$solid_food=c('off','on')[as.numeric(as.factor(data$Solid_Food))]
png(paste0(resdir,'figures/relative_abundance_Clostridiales.png'))
plot_timepoints(data,'Age_at_Collection','Subject_ID','o__Clostridiales_relative','solid_food','Case_Control',shape_values =c(4,1),color_gradientn =rainbow(3),shape_labs = 'solid food',color_labs = 'relative abundance of Clostridiales',xlab='age at collection (day)',ylab='subject ID')
dev.off()

