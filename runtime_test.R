rm(list = ls())
#registerDoParallel(100)

############################## Test for multivariate models ################################
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
sourceCpp(file="/userhome/home/fun_indep.cpp")
 
source("function.R")
# set Anaconda environment
use_condaenv("base")  
# set working directory
setwd("/userhome/home/KernelLogrankTest-master") 
py_run_string("import sys;sys.path.append('/userhome/home/KernelLogrankTest-master');
sys.path.append('/userhome/home/kerpy-master')")  

# import kernel_logrank
kernel_logrank <- import("kernel_logrank")

# import CPH test and logrank test functions
cph_test <- import("kernel_logrank.tests.cph_test", convert = TRUE)$cph_test
wild_bootstrap_test_logrank_covariates <- import("kernel_logrank.tests.wild_bootstrap_LR", convert = TRUE)$wild_bootstrap_test_logrank_covariates

#########  run time of permutation tests #################

############# the linear model #########################################
rho <- 0.4
Sigma <- matrix(rho,3,3)
diag(Sigma) <- 1
beta <- matrix(c(-1,-1,1),3,1) 
mu <- rep(0, 3)
B <- 500
size <- 1000

#################### paralell ######################
start_time<-Sys.time()
n <- 100
  registerDoParallel(20)
  ########### for cluster use ############
  output2<-foreach(m = seq(1,size,by=1),.combine = 'rbind') %dopar% {
    X<- mvrnorm(n, mu, Sigma)
    y <-    0.8*c(X%*%beta) + rnorm(n,0,1)
    C <-  rnorm(n,1.4,0.3)+ rnorm(n,1,1)  # censoring variable=0.45
    delta <-(y<=C)
    Cr <- 1-mean(delta) # censoring rates
    Y<- pmin(y,C) # observed variable
    
    ############## generate censoring ############
    delta <-(y<=C)
    Cr <- 1-mean(delta) # censoring rates
    Y<- pmin(y,C) # observed variable
    
    ########### cox PH likelihood ratio test #####
    time.cox <- system.time({cph_test(X=as.matrix(X),z=array(Y),d=array(delta))})[3]
    ########### KLR ###########
    time.KLR <-system.time({wild_bootstrap_test_logrank_covariates(
      X=as.matrix(X), z=array(Y),d=array(delta),kernels_x='gau',kernel_z='gau',seed=as.integer(1))})[3]
    ########### IPCW ###########
    time.distance <- system.time({ipcw.dcov.test(cbind(Y,delta), X,B=B)})[3]
    ########### proposed tests ###########
    Y<- (Y-mean(Y))/sd(Y)# standardization
    ############## different cumulative weight function ##########################
    # U(-1,1)  
    Fw <- function(x){punif(x, -1,1)}
    time_unif <-  system.time({permu_pvalue_fast(Y, X, delta,Fw,B)})[3]
    # U(-5,5)
    Fw <- function(x){punif(x, -5,5)}
    time_unif50 <-  system.time({permu_pvalue_fast(Y, X, delta,Fw,B)})[3]
    # standard normal 
    Fw <- function(x) {pnorm(x,0,1)}
    time_norm <-  system.time({permu_pvalue_fast(Y, X, delta,Fw,B)})[3]
    # normal 
    Fw <- function(x) {pnorm(x,-1,5)}
    time_norm50 <- system.time({permu_pvalue_fast(Y, X, delta,Fw,B)})[3]
    # t(2)  
    Fw <- function(x) {pt(x,2)}
    time_t2 <- system.time({permu_pvalue_fast(Y, X, delta,Fw,B)})[3]
    # t(3)  
    Fw <- function(x) {pt(x,3)}
    time_t3 <- system.time({permu_pvalue_fast(Y, X, delta,Fw,B)})[3]
    
    res<- c(time.cox,time.KLR,time.distance, time_norm,time_norm50,time_unif,time_unif50,
            time_t2,time_t3,Cr)
    return(res)
  }

end_time<-Sys.time()
times <- end_time-start_time
print(times)
print(round(colMeans(output2),4))
round(colSds(output2),4)

