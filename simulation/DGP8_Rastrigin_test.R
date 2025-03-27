rm(list = ls())
############################## DGP: Rastrigin models ################################
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
sys.path.append('/userhome/kerpy-master')") #11

# import kernel_logrank
kernel_logrank <- import("kernel_logrank")

# import CPH test and logrank test functions
cph_test <- import("kernel_logrank.tests.cph_test", convert = TRUE)$cph_test
wild_bootstrap_test_logrank_covariates <- import("kernel_logrank.tests.wild_bootstrap_LR", convert = TRUE)$wild_bootstrap_test_logrank_covariates

df45 <- matrix(NA,7,10)
df65 <- matrix(NA,7,10)
#nset <- c(50,150)
nset <- c(50,100,150,200,250,300,350)

ID <-1# ID <- 0 for size

registerDoParallel(150)
B <- 500
size <- 1000
start_time<-Sys.time()
for (j in 1:length(nset)){
  n <- nset[j]
  output45<-foreach(m = seq(1,size,by=1),.combine = 'rbind') %dopar% {
    A <- 5
    # define Rastrigin function
    rastrigin <- function(x) {
      return(A*2  + (x[1]^2 - A * cos(2 * pi * x[1])) + (x[2]^2 - A * cos(2 * pi * x[2])))
    }
    # generate data
    X1 <- runif(n,-1,1)
    X2 <- runif(n,-1,1)
    X <- cbind(X1,X2)
    y <-  apply(X,1,rastrigin)*ID +   rnorm(n,0,1)
    C <-  rnorm(n,2,4)- rnorm(n,0,1)+ 5*rnorm(n,1.9,0.1) # censoring variable=0.45
    #C <-  rnorm(n,2,4)- rnorm(n,0,1)+ 5*rnorm(n,-0.3,0.1) # censoring variable ID <-0
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
    A <- 5
    rastrigin <- function(x) {
      return(A*2  + (x[1]^2 - A * cos(2 * pi * x[1])) + (x[2]^2 - A * cos(2 * pi * x[2])))
    }
    # generate data
    X1 <- runif(n,-1,1)
    X2 <- runif(n,-1,1)
    X <- cbind(X1,X2)
    y <-  apply(X,1,rastrigin)*ID +   rnorm(n,0,1)
    C <-  rnorm(n,2,4)- rnorm(n,0,1)+ 5*rnorm(n,1.2,0.1) # censoring variable=0.65
    #C <-  rnorm(n,2,4)- rnorm(n,0,1)+ 5*rnorm(n,-0.75,0.1) # censoring variable ID <-0
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
  df45[j,] <-c(colMeans(output45[,1:9]<=0.05),mean(output45[,10]))
  df65[j,] <-c(colMeans(output65[,1:9]<=0.05),mean(output65[,10]))
}
end_time<-Sys.time()
times <- end_time-start_time
print(times)
print(df45)
print(df65)

