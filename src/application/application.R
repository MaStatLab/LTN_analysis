#!/usr/bin/env Rscript
library(LTN)
argv=commandArgs(TRUE)
WORK_DIR=argv[1]
#source(paste0(WORK_DIR,"/src/utility/utility.R"))
#source(paste0(WORK_DIR,"/src/utility/mixed_effects.R"))
niter=as.numeric(argv[2])
covariate_index=as.numeric(argv[3])
#model_index=as.numeric(argv[4]) 
model_index=2
#lambda=as.numeric(argv[5])
lambda=as.numeric(argv[4])
#gprior_m=100
#modelspec=c('diagonal','sparse')[model_index]
dietcovariates=c("BF","Solid_Food","Eggs","Fish","Soy_Prod","Rye","Barley","Buckwheat_Millet")
covariates=c("Case_Control","post_seroconversion",dietcovariates)
baseline=c("control",F,rep("false",length(dietcovariates)))
covariate=covariates[covariate_index]
datadir=paste0(WORK_DIR,"/cache/")
resdir=paste0(WORK_DIR,"/results/application/")
system(paste0('mkdir -p ',resdir,'/pmap'))
system(paste0('mkdir -p ',resdir,'/pjap'))
system(paste0('mkdir -p ',resdir,'/alpha'))
system(paste0('mkdir -p ',resdir,'/time'))
system(paste0('mkdir -p ',resdir,'/memory'))
filenam=paste0(covariate,'_lambda',lambda,'.RData')
if (!file.exists(paste0(resdir,'/pjap/',filenam))){
data=readRDS(paste0(datadir,"ps_otu50.RData"))
if (covariate=="Case_Control"){
  formula=as.formula(paste0('~',covariate,'+log(Age_at_Collection)+Country+post_seroconversion+Gender+',paste(dietcovariates,collapse = '+'),'+(1 | Subject_ID)'))
}else if (covariate=="post_seroconversion"){
  formula=as.formula(paste0('~',covariate,'+log(Age_at_Collection)+Country+Gender+Case_Control+',paste(dietcovariates,collapse = '+'),'+(1 | Subject_ID)'))
}else if (covariate %in% dietcovariates){
  formula=as.formula(paste0('~',covariate,'+log(Age_at_Collection)+Country+post_seroconversion+Gender+Case_Control+',paste(dietcovariates[-which(dietcovariates==covariate)],collapse = '+'),'+(1 | Subject_ID)'))
}else{warning('Unknown covariate')}
st<-system.time(result<-ltnme(formula=formula,
      data=data,
      test_baseline=baseline[covariate_index],
      niter=niter,
      reffcov = model_index,
      SEED = 1,
      save_alpha_only = T,
      gprior_m = 100,
      pnull = 0.5,
      lambda = lambda))
cat(paste0(covariate,': ',st[3],'s'))
saveRDS(result$PMAP,paste0(resdir,'/pmap/',filenam))
saveRDS(result$PJAP,paste0(resdir,'/pjap/',filenam))
saveRDS(result$alpha_mean,paste0(resdir,'/alpha/',filenam))
saveRDS(st,paste0(resdir,'/time/',filenam))
saveRDS(object.size(result),paste0(resdir,'/memory/',filenam))
}

