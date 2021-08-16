#!/usr/bin/env Rscript

library(LTN)
argv=commandArgs(TRUE)
WORK_DIR=argv[1]
lambda=argv[2]
nsim1=argv[3]
#source(paste0(WORK_DIR,"/src/utility/utility.R"))
if (lambda=='0'){lambda='hp'}
risk=matrix(-1,4,4,dimnames = list(c("Frobenius","L1","L_inf","spectral"),c("LN-hub","LN-block","LN-sparse","DTM")))
# DTM
cachedir=paste0(WORK_DIR,"/cache/covariance_estimation/dtm/ltn/")
clrcov0=readRDS(paste0(WORK_DIR,"/cache/covariance_estimation/dtm/clrcov_true.RData"))
files=grep(paste0('lambda',lambda,'.RData'),list.files(cachedir,full.names = T),value = T)
loss=do.call(rbind,lapply(files, function(x){matloss(cov2cor(readRDS(x))-cov2cor(clrcov0))}))
risk[,'DTM']=(colSums(loss)/nrow(loss))
# LN
cachedir=paste0(WORK_DIR,"/cache/covariance_estimation/ln/ltn/")
for (m in c('hub','block','sparse')){
  files=grep(paste0(m,'_lambda',lambda,'.RData'),list.files(cachedir,full.names = T),value = T)
  loss=do.call(rbind,lapply(1:nsim1, function(x){
    clrcov0=readRDS(paste0(WORK_DIR,"/cache/covariance_estimation/ln/",m,"clrcov_",x,".RData"))
    clrcov1=readRDS(grep(paste0('i',x,'_'),files,value=T))
    matloss(cov2cor(clrcov1)-cov2cor(clrcov0))}))
  risk[,paste0("LN-",m)]=(colSums(loss)/nrow(loss))
}
if (lambda=='hp'){lambda='GammaPrior'}
result_dir=paste0(WORK_DIR,"/results/covariance_estimation/")
system(paste0('mkdir -p ',result_dir))
saveRDS(risk,paste0(result_dir,"/loss_LTN_lambda_",lambda,".RData"))

#--------------coat--------------------
risk=matrix(-1,4,4,dimnames = list(c("Frobenius","L1","L_inf","spectral"),c("LN-hub","LN-block","LN-sparse","DTM")))
cachedir=paste0(WORK_DIR,"/cache/covariance_estimation/dtm/coat/")
clrcov0=readRDS(paste0(WORK_DIR,"/cache/covariance_estimation/dtm/clrcov_true.RData"))
files=list.files(cachedir,full.names = T)
loss=do.call(rbind,lapply(files, function(x){matloss(cov2cor(readRDS(x))-cov2cor(clrcov0))}))
risk[,'DTM']=(colSums(loss)/nrow(loss))
cachedir=paste0(WORK_DIR,"/cache/covariance_estimation/ln/coat/")
for (m in c('hub','block','sparse')){
  files=grep(m,list.files(cachedir,full.names = T),value = T)
  loss=do.call(rbind,lapply(1:nsim1, function(x){
    clrcov0=readRDS(paste0(WORK_DIR,"/cache/covariance_estimation/ln/",m,"clrcov_",x,".RData"))
    clrcov1=readRDS(grep(paste0('est',x),files,value=T))
    matloss(cov2cor(clrcov1)-cov2cor(clrcov0))}))
  risk[,paste0("LN-",m)]=(colSums(loss)/nrow(loss))
}
saveRDS(risk,paste0(result_dir,"/loss_COAT.RData"))

