##########################################################
#Protein
##########################################################
#srun -t 8:00:00 -p interact -N 1 -n 1 --mem=32gb --pty /bin/bash
cd /overflow/KernelLogrankTest-master/data
module load r/4.4.0;R
rm(list = ls())
library(MASS)
library(reticulate)
library(RANN)
library(FOCI)
library(cubature)
library(matrixStats)
library(doParallel)
library(foreach)
library(dcortools)
library(RcppArmadillo)
library(pracma)
library(survival)
library(RcppEigen)
library(Rcpp)
sourceCpp(file="/overflow/fun_indep.cpp") #18
source('/overflow/user/KernelLogrankTest-master/function.r')
use_condaenv("base")  
setwd("/overflow/user/projects/KernelLogrankTest-master") #17
py_run_string("import sys;sys.path.append('/overflow/user/projects/KernelLogrankTest-master');
sys.path.append('/overflow/user/projects/kerpy-master')") #17
py_run_string("import sys; sys.path.append('/overflow/user/projects')") #18
kernel_logrank <- import("kernel_logrank")
cph_test <- import("kernel_logrank.tests.cph_test", convert = TRUE)$cph_test
wild_bootstrap_test_logrank_covariates <- import("kernel_logrank.tests.wild_bootstrap_LR", convert = TRUE)$wild_bootstrap_test_logrank_covariates
generate_synthetic_data <- import("kernel_logrank.utils.generate_synthetic_data", convert = TRUE)$generate_synthetic_data
library(ggplot2)
setwd('./data') #or 'PseudoRealData/'
A2=read.csv('surv.csv');M1=read.csv('protein.csv')
##########################################################
Method=c('CPH','KLR','IPCW','N01','U1')
Method=Method[2] #change method
#if.log='LOG';
if.log='NoLOG';
if.scale=0
if(Method%in%c('N01','U1'))if.scale=1
if(if.log=='LOG'){Y <- log(A2$time)}
if(if.log=='NoLOG'){Y <- (A2$time)}
if(if.scale==1){Y <- (Y-mean(Y))/sd(Y)}
delta <- A2$censor
LROI=dim(M1)[2]
UROI=colnames(M1)
Coxph=as.data.frame(matrix(NA,LROI,4))
Coxph[,1]=UROI
for(ii in 1:LROI)
{
	if(ii%%100==0)print(LROI-ii)
	X=as.numeric(M1[,ii])
	X[is.na(X)]=median(X,na.rm=T)
	if(Method=='CPH')temp=cph_test(X=as.matrix(X),z=array(Y),d=array(delta))
	if(Method=='U1')temp=permu_pvalue_fast1(Y, X, delta,Fw=function(x){punif(x,min=-1,max=1)}) 
	if(Method=='N01')temp=permu_pvalue_fast1(Y, X, delta) 
	if(Method=='IPCW')temp=ipcw.dcov.test1(cbind(Y,delta), X)
	if(Method=='KLR')temp=wild_bootstrap_LR1(array(as.matrix(X),c(length(Y),1)),array(Y),array(delta),'gau','gau',1,500)
	Coxph[ii,2]=temp[1];if(length(temp)>1)Coxph[ii,4]=temp[2]
	if(temp[1]<0.0005)print(c(ii,temp))
}
Coxph[,3]=-log10(p.adjust(Coxph[,2],'fdr'))
colnames(Coxph)=c('ROI','Pval','FDR','COR')
write.csv(Coxph,file=paste0(if.log,'_',Method,'.csv'),quote=F)