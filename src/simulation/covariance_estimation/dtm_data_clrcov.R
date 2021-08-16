#!/usr/bin/env Rscript

library(LTN)
argv=commandArgs(TRUE)
WORK_DIR=argv[1]
nmc=as.numeric(argv[2])
cachedir=paste0(WORK_DIR,"/cache/covariance_estimation/dtm/")
system(paste0('mkdir -p ',cachedir))
# source(paste0(WORK_DIR,"/src/utility/utility.R"))
if(!file.exists(paste0(cachedir,"dtm_diab_MoM.RData"))){
  inputdata=readRDS(paste0(WORK_DIR,"/cache/ps_sim.RData"))
  cnt=inputdata$cnt
  tree=inputdata$tree
  est=dtm_mom(cnt,tree)
  saveRDS(est,paste0(cachedir,"dtm_diab_MoM.RData"))
}

if (!file.exists(paste0(cachedir,"clrcov_true.RData"))){
input_data1=readRDS(paste0(cachedir,"dtm_diab_MoM.RData"))
input_data2=readRDS(paste0(WORK_DIR,"/cache/ps_sim.RData"))
theta=as.vector(input_data1[[1]])
tau=as.vector(input_data1[[2]])
tree=input_data2$tree
clrcov_true=clrcov_dtm_sim_log(nmc,tree,theta,tau,SSS=1,savesamp=F,dir=NULL)
if (sum(is.na(clrcov_true))+sum(is.infinite(clrcov_true))==0){
  saveRDS(clrcov_true,paste0(cachedir,"clrcov_true.RData"))
  print('finished calculating the true clr covariance in DTM example')
} else{
  warning('NA in clr covariance')
}
}


