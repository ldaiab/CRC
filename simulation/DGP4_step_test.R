rm(list = ls())

############################## DGP: step function ################################
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
# source r function
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

df45 <- matrix(NA,7,10)
df65 <- matrix(NA,7,10)
ID <- 1 # ID <- 0 for size
B <- 500
size <- 1000
nset <- c(50,100,150,200,250,300,350)

registerDoParallel(150)
start_time<-Sys.time()
for (j in 1:length(nset)){
  n <- nset[j]
  output45<-foreach(m = seq(1,size,by=1),.combine = 'rbind') %dopar% {
    ## Step-function  
    X <- runif(n,-2,2)
    g <- function(x){
      if ((x >= -2)& (x < -1.5)){
        y <- -2
      }
      if ((x >= -1.5)& (x < -1)){
        y <- 4
      }
      if ((x >= -1)& (x < -0.5)){
        y <- -4
      }
      if((x >= -0.5)& (x<=0)){
        y <- 2
      }
      if ((x >= 0)& (x < 0.5)){
        y <- -2
      }
      if ((x >= 0.5)& (x < 1)){
        y <- 2.5
      }
      if ((x >= 1)& (x < 1.5)){
        y <- -3
      }
      if((x >= 1.5)& (x<=2)){
        y <- 3
      }
      return(y)
    }
    y <- 0.5*sapply(X,g)*ID+ rnorm(n,0,1)
    #C <-  rnorm(n,2.3,0.5)-  rnorm(n,2.2,0.1)# censoring variable ID <-0
    C <-  rnorm(n,2.3,0.5)-  rnorm(n,2,0.1)# censoring variable
    delta <-(y<=C)
    Cr <- 1-mean(delta) # censoring rates
    Y<- pmin(y,C) # observed variable
    
    # cox PH likelihood ratio test  
    pvalue.cox <- cph_test(X=as.matrix(X),z=array(Y),d=array(delta))
    
    # KLR 
    pvalue.KLR <- wild_bootstrap_test_logrank_covariates(
      X=as.matrix(X), z=array(Y),d=array(delta),kernels_x='gau',kernel_z='gau')[[2]]
    
    # IPCW  
    dis.pvalue<- ipcw.dcov.test(Y=cbind(Y,delta), X=X,B=B)$pvalue
    
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
  output65<-foreach(m = seq(1,size,by=1),.combine = 'rbind') %dopar% {
    X <- runif(n,-2,2)
    g <- function(x){
      if ((x >= -2)& (x < -1.5)){
        y <- -2
      }
      if ((x >= -1.5)& (x < -1)){
        y <- 4
      }
      if ((x >= -1)& (x < -0.5)){
        y <- -4
      }
      if((x >= -0.5)& (x<=0)){
        y <- 2
      }
      if ((x >= 0)& (x < 0.5)){
        y <- -2
      }
      if ((x >= 0.5)& (x < 1)){
        y <- 2.5
      }
      if ((x >= 1)& (x < 1.5)){
        y <- -3
      }
      if((x >= 1.5)& (x<=2)){
        y <- 3
      }
      return(y)
    }
    y <- 0.5*sapply(X,g)*ID+ rnorm(n,0,1)
    C <-  rnorm(n,2.3,0.5)-  rnorm(n,3.2,0.1)# censoring variable
    #C <-  rnorm(n,2.3,0.5)-  rnorm(n,2.7,0.1) # censoring variable ID <- 0
    delta <-(y<=C)
    Cr <- 1-mean(delta) # censoring rates
    Y<- pmin(y,C)  
    
    # cox PH likelihood ratio test  
    pvalue.cox <- cph_test(X=as.matrix(X),z=array(Y),d=array(delta))
    
    # KLR  
    pvalue.KLR <- wild_bootstrap_test_logrank_covariates(
      X=as.matrix(X), z=array(Y),d=array(delta),kernels_x='gau',kernel_z='gau')[[2]]
    
    # IPCW 
    dis.pvalue<- ipcw.dcov.test(Y=cbind(Y,delta), X=X,B=B)$pvalue
    
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
  df45[j,] <-c(colMeans(output45[,1:9]<=0.05),mean(output45[,10]))
  df65[j,] <-c(colMeans(output65[,1:9]<=0.05),mean(output65[,10]))
 }
end_time<-Sys.time()
times <- end_time-start_time
print(times)
print(df45)
print(df65)
