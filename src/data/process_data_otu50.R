#!/usr/bin/env Rscript
library(phyloseq)
library(ape)
library(data.tree)
library(LTN)
argv=commandArgs(TRUE)
WORK_DIR=argv[1]
system(paste0('mkdir -p ',WORK_DIR,'/cache'))
if (!file.exists(paste0(WORK_DIR,"/cache/ps_otu50.RData"))){
  raw_data = read.csv(paste0(WORK_DIR,"/res/diabimmune_t1d_16s_otu_table.txt"), sep = "\t", row.names = 1 )
  otutab0=apply(as.matrix(raw_data[,-ncol(raw_data)]),1:2,as.numeric)
  rownames(otutab0)=rownames(raw_data)
  colnames(otutab0)=colnames(raw_data)[-ncol(raw_data)]
  taxtab0=do.call(rbind,sapply(as.character(raw_data[,ncol(raw_data)]),function(x){strsplit(x,';  ')}))
  rownames(taxtab0)=rownames(otutab0)
  colnames(taxtab0)=c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
  ps0=phyloseq(otu_table(otutab0,taxa_are_rows = T),tax_table(taxtab0))
  ra=apply(otutab0,2,function(x){x/sum(x)})
  otu100=names(sort(rowSums(ra),decreasing = T)[1:100])
  ps100=prune_taxa(otu100,ps0)
  tax100=tax_table(ps100)
  tree=tax_tree(tax100)
  otu_preorder=tree$tip.label
  otutab100=otu_table(ps100)[otu_preorder,]
  taxtab100=tax100[otu_preorder,]
  ps100=phyloseq(otu_table(t(otutab100),taxa_are_rows = F),tax_table(taxtab100),tree)
  # sample data
  load(paste0(WORK_DIR,"/res/diabimmune_t1d_16s_metadata.rdata")) # md_16S
  rownames(md_16S)=as.character(md_16S$G_id)
  md_16S=md_16S[rownames(otu_table(ps100)),]
  ps100=merge_phyloseq(ps100,sample_data(md_16S))
  otu50=names(sort(rowSums(ra),decreasing = T)[1:50])
  ps50=prune_taxa(otu50,ps100)
  # subject data
  mmc2=readxl::read_excel(paste0(WORK_DIR,"/res/mmc2.xlsx"))
  mmc2=data.frame(mmc2)
  rownames(mmc2)=mmc2$ID..E.Espoo..Finland..T.Tartu..Estonia
  sampdat=sample_data(ps50)
  sampdat$sero_age=mmc2[as.character(sampdat$Subject_ID),"Age.at.sero.conversion..yrs"]*365
  sampdat$sero_age[is.na(sampdat$sero_age)]=max(sampdat$Age_at_Collection)+1
  sampdat$post_seroconversion=sampdat$Age_at_Collection>sampdat$sero_age
  sample_data(ps50)=sampdat
  saveRDS(ps50,paste0(WORK_DIR,"/cache/ps_otu50.RData"))
}
yyl=seqtab2y(otu_table(ps50),phy_tree(ps50))
saveRDS(yyl,paste0(WORK_DIR,"/cache/yyl_otu50.RData"))
