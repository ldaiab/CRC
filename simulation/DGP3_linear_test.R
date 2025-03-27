rm(list = ls())
############################## DGP3: permutation test for linear model ################################
library(MASS)
library(matrixStats)
library(reticulate) # for import python
library(RANN)
library(doParallel)
library(foreach)
library(dcortools)
library(RcppEigen)
library(Rcpp)
# input cpp files
# source r functions
sourceCpp(file="/userhome/fun_indep.cpp")

source("function.R")
# set Anaconda environment
use_condaenv("base")  
# set working directory
setwd("/userhome/KernelLogrankTest-master") 
py_run_string("import sys;sys.path.append('/userhome/KernelLogrankTest-master');
sys.path.append('/userhome/kerpy-master')") 

# import kernel_logrank
kernel_logrank <- import("kernel_logrank")

# import CPH test and logrank test functions
cph_test <- import("kernel_logrank.tests.cph_test", convert = TRUE)$cph_test
wild_bootstrap_test_logrank_covariates <- import("kernel_logrank.tests.wild_bootstrap_LR", convert = TRUE)$wild_bootstrap_test_logrank_covariates

# model setting  
rho <- 0.2
Sigma <- matrix(rho,3,3)
diag(Sigma) <- 1
beta <- matrix(c(-1,-1,1),3,1) 

mu <- rep(0, 3)
ID <-1 # ID <- 0 for size
#nset <- c(50,150)
nset <- c(50,100,150,200,250,300,350)

df45 <- matrix(NA,7,10) # alpha = 0.05,Cr=45%
df65 <- matrix(NA,7,10) # alpha = 0.05,Cr=65%
B <- 500
size <- 1000
start_time<-Sys.time()
registerDoParallel(150)
for (j in 1:length(nset)){
  n <- nset[j]
  ########### CR = 45% ############
  output2<-foreach(m = seq(1,size,by=1),.combine = 'rbind') %dopar% {
    X<- mvrnorm(n, mu, Sigma)
    y <-    0.8*c(X%*%beta)*ID + rnorm(n,0,1)
    C <-  -rnorm(n,0.8,0.3)+ rnorm(n,1,1)  # censoring variable=0.45
    delta <-(y<=C)
    Cr <- 1-mean(delta) # censoring rates
    Y<- pmin(y,C) # observed variable
    
    # cox PH likelihood ratio test  
    pvalue.cox <- cph_test(X=as.matrix(X),z=array(Y),d=array(delta))
    pvalue.cox
    
    #  KLR  
    pvalue.KLR <- wild_bootstrap_test_logrank_covariates(
      X=as.matrix(X), z=array(Y),d=array(delta),kernels_x='gau',kernel_z='gau')[[2]]
    
    # IPCW  
    dis.pvalue<- ipcw.dcov.test(cbind(Y,delta), X,B=B)$pvalue
    
    # proposed method with different weight function 
    Y<- (Y-mean(Y))/sd(Y)# standardization
    # U(-1,1)  
    Fw <- function(x){punif(x, -1,1)}
    pvalue.unif1 <- permu_pvalue_fast(Y, X, delta,Fw,B)[1]
    # U(-5,5)
    Fw <- function(x){punif(x, -5,5)}
    pvalue.unif50 <-  permu_pvalue_fast(Y, X, delta,Fw,B)[1]
    # standard normal 
    Fw <- function(x) {pnorm(x,0,1)}
    pvalue.snorm <-  permu_pvalue_fast(Y, X, delta,Fw,B)[1]
    # normal 
    Fw <- function(x) {pnorm(x,-1,5)}
    pvalue.norm <- permu_pvalue_fast(Y, X, delta,Fw,B)[1]
    # t(2)  
    Fw <- function(x) {pt(x,2)}
    pvalue.t2 <- permu_pvalue_fast(Y, X, delta,Fw,B)[1]
    # t(3)  
    Fw <- function(x) {pt(x,3)}
    pvalue.t3 <- permu_pvalue_fast(Y, X, delta,Fw,B)[1]
    # summary results
    res<- c(pvalue.cox, pvalue.KLR, dis.pvalue,pvalue.snorm,pvalue.norm,
            pvalue.unif1,pvalue.unif50,pvalue.t2,pvalue.t3,Cr)
    return(res)
  }
  ########### CR = 65% ############
  output4<-foreach(m = seq(1,size,by=1),.combine = 'rbind') %dopar% {
    X<- mvrnorm(n, mu, Sigma)
    y <-   0.8*c(X%*%beta)*ID + rnorm(n,0,1)
    C <-  -rnorm(n,1.8,0.3)+ rnorm(n,1,1)  # censoring variable=0.65, ID = 1
    #C <-  -rnorm(n,1.6,0.3)+ rnorm(n,1,1)  # censoring variable=0.65, ID = 0
    delta <-(y<=C)
    Cr <- 1-mean(delta) # censoring rates
    Y<- pmin(y,C) # observed variable
    
    # cox PH likelihood ratio test 
    pvalue.cox <- cph_test(X=as.matrix(X),z=array(Y),d=array(delta))
    pvalue.cox
    
    # KLR
    pvalue.KLR <- wild_bootstrap_test_logrank_covariates(
      X=as.matrix(X), z=array(Y),d=array(delta),kernels_x='gau',kernel_z='gau')[[2]]
    
    # IPCW  
    dis.pvalue<- ipcw.dcov.test(cbind(Y,delta), X,B=B)$pvalue
  
    # proposed tests with different weight functions
    Y<- (Y-mean(Y))/sd(Y)# standardization
    # U(-1,1)  
    Fw <- function(x){punif(x, -1,1)}
    pvalue.unif1 <- permu_pvalue_fast(Y, X, delta,Fw,B)[1]
    # U(-5,5)
    Fw <- function(x){punif(x, -5,5)}
    pvalue.unif50 <-  permu_pvalue_fast(Y, X, delta,Fw,B)[1]
    # standard normal 
    Fw <- function(x) {pnorm(x,0,1)}
    pvalue.snorm <-  permu_pvalue_fast(Y, X, delta,Fw,B)[1]
    # normal 
    Fw <- function(x) {pnorm(x,-1,5)}
    pvalue.norm <- permu_pvalue_fast(Y, X, delta,Fw,B)[1]
    # t(2)  
    Fw <- function(x) {pt(x,2)}
    pvalue.t2 <- permu_pvalue_fast(Y, X, delta,Fw,B)[1]
    # t(3)  
    Fw <- function(x) {pt(x,3)}
    pvalue.t3 <- permu_pvalue_fast(Y, X, delta,Fw,B)[1]
    # summary results
    res<- c(pvalue.cox, pvalue.KLR, dis.pvalue,pvalue.snorm,pvalue.norm,
            pvalue.unif1,pvalue.unif50,pvalue.t2,pvalue.t3,Cr)
    return(res)
  }
  df45[j,] <-c(colMeans(output2[,1:9]<=0.05),mean(output2[,10]))
  df65[j,] <-c(colMeans(output4[,1:9]<=0.05),mean(output4[,10]))
 }
end_time<-Sys.time()
times <- end_time-start_time
print(times)
print(df45) # power
print(df65) # power


