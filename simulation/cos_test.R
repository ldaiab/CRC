rm(list = ls())
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
sourceCpp(file="/code/fun_indep.cpp") #18

source("function.R")
# set Anaconda environment
use_condaenv("base")  
# set working directory
setwd("/home/KernelLogrankTest-master")  #18
py_run_string("import sys; sys.path.append('/KernelLogrankTest-master')") #18

# import kernel_logrank
kernel_logrank <- import("kernel_logrank")

# import CPH test and logrank test functions
cph_test <- import("kernel_logrank.tests.cph_test", convert = TRUE)$cph_test
wild_bootstrap_test_logrank_covariates <- import("kernel_logrank.tests.wild_bootstrap_LR", convert = TRUE)$wild_bootstrap_test_logrank_covariates

df20<- matrix(NA,7,10)  # cos(2*pi*x), Cr=0
df30<- matrix(NA,7,10)  # cos(3*pi*x), Cr=0
df80<- matrix(NA,7,10)  # cos(8*pi*x), Cr=0

df2_30<- matrix(NA,7,10) # cos(2*pi*x), Cr=30%
df3_30<- matrix(NA,7,10) # cos(3*pi*x), Cr=30%
df8_30<- matrix(NA,7,10) # cos(8*pi*x), Cr=30%

df2_60<- matrix(NA,7,10) # cos(2*pi*x), Cr=60%
df3_60<- matrix(NA,7,10) # cos(3*pi*x), Cr=60%
df8_60<- matrix(NA,7,10) # cos(8*pi*x), Cr=60%

nset <- c(50,100,150,200,250,300,350)

registerDoParallel(150)
B <- 500
size <- 1000
start_time<-Sys.time()
# cos(2*pi*x)
for (j in 1:length(nset)){
  n <- nset[j]
  ########### for  cluster use ############
  output0<-foreach(m = seq(1,size,by=1),.combine = 'rbind') %dopar% {
    X <- runif(n,-1,1)
    y <- cos(2*pi*X)+ .5*rnorm(n,0,1)
    C <-  1* rnorm(n,1.1,0.5)- rnorm(n,1,0.1)# censoring variable
    C <- 1e7 # CR=0%
    ############## generate censoring ############
    delta <-(y<=C)
    Cr <- 1-mean(delta) # censoring rates
    Y<- pmin(y,C) # observed variable
    
    ########### cox PH likelihood ratio test #####
    pvalue.cox <- cph_test(X=as.matrix(X),z=array(Y),d=array(delta))
    ########### KLR ######
    # testing code X=X, z=Y(time), d=delta
    # dim=1
    pvalue.KLR <- wild_bootstrap_test_logrank_covariates(
      X=as.matrix(X), z=array(Y),d=array(delta),kernels_x='gau',kernel_z='gau')[[2]]
    ########### IPCW ######
    dis.pvalue<- ipcw.dcov.test(Y=cbind(Y,delta), X=X,B=B)$pvalue
    
    ########### proposed tests ###################
    Y<- (Y-mean(Y))/sd(Y)# standardization
    ############## different weight function ##########################
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
  output3<-foreach(m = seq(1,size,by=1),.combine = 'rbind') %dopar% {
    X <- runif(n,-1,1)
    y <- cos(2*pi*X)+  0.5*rnorm(n,0,1)
    C <-  1* rnorm(n,1.6,0.5)- rnorm(n,1,0.1)# censoring variable
    ############## generate censoring ############
    delta <-(y<=C)
    Cr <- 1-mean(delta) # censoring rates
    Y<- pmin(y,C) 
    
    ########### cox PH likelihood ratio test #####
    pvalue.cox <- cph_test(X=as.matrix(X),z=array(Y),d=array(delta))
    ########### KLR ######
    # testing code X=X, z=Y(time), d=delta
    # dim=1

    pvalue.KLR <- wild_bootstrap_test_logrank_covariates(
      X=as.matrix(X), z=array(Y),d=array(delta),kernels_x='gau',kernel_z='gau')[[2]]

    ########### IPCW ######
    dis.pvalue<- ipcw.dcov.test(Y=cbind(Y,delta), X=X,B=B)$pvalue
    
    ########### proposed tests ###################
    Y<- (Y-mean(Y))/sd(Y)# standardization
    ############## different weight function ##########################
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
  output6<-foreach(m = seq(1,size,by=1),.combine = 'rbind') %dopar% {
    X <- runif(n,-1,1)
    y <- cos(2*pi*X)+  0.5*rnorm(n,0,1)
    C <-  1* rnorm(n,0.75,0.5)- rnorm(n,1,0.1)# censoring variable
    ############## generate censoring ############
    delta <-(y<=C)
    Cr <- 1-mean(delta) # censoring rates
    Y<- pmin(y,C) 
    
    ########### cox PH likelihood ratio test #####
    pvalue.cox <- cph_test(X=as.matrix(X),z=array(Y),d=array(delta))
    ########### KLR ######
    # testing code X=X, z=Y(time), d=delta
    # dim=1
    
    pvalue.KLR <- wild_bootstrap_test_logrank_covariates(
      X=as.matrix(X), z=array(Y),d=array(delta),kernels_x='gau',kernel_z='gau')[[2]]
    
    ########### IPCW ######
    dis.pvalue<- ipcw.dcov.test(Y=cbind(Y,delta), X=X,B=B)$pvalue
    
    ########### proposed tests ###################
    Y<- (Y-mean(Y))/sd(Y)# standardization
    ############## different weight function ##########################
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
  df20[j,] <-c(colMeans(output0[,1:9]<=0.05),mean(output0[,10]))
  df2_30[j,] <-c(colMeans(output3[,1:9]<=0.05),mean(output3[,10]))
  df2_60[j,] <-c(colMeans(output6[,1:9]<=0.05),mean(output6[,10]))
}
end_time<-Sys.time()
times <- end_time-start_time
print(times)
print(df20)
print(df2_30)
print(df2_60)
# cos(3*pi*x)

registerDoParallel(150)
B <- 500
size <- 1000
start_time<-Sys.time()
# cos(3*pi*x)
for (j in 1:length(nset)){
  n <- nset[j]
  ########### for  cluster use ############
  output0<-foreach(m = seq(1,size,by=1),.combine = 'rbind') %dopar% {
    X <- runif(n,-1,1)
    y <- cos(3*pi*X)+ .5*rnorm(n,0,1)
    C <-  1* rnorm(n,1.1,0.5)- rnorm(n,1,0.1)# censoring variable
    C <- 1e7 # CR=0%
    ############## generate censoring ############
    delta <-(y<=C)
    Cr <- 1-mean(delta) # censoring rates
    Y<- pmin(y,C) # observed variable
    
    ########### cox PH likelihood ratio test #####
    pvalue.cox <- cph_test(X=as.matrix(X),z=array(Y),d=array(delta))
    ########### KLR ######
    # testing code X=X, z=Y(time), d=delta
    # dim=1
    pvalue.KLR <- wild_bootstrap_test_logrank_covariates(
      X=as.matrix(X), z=array(Y),d=array(delta),kernels_x='gau',kernel_z='gau')[[2]]
    ########### IPCW ######
    dis.pvalue<- ipcw.dcov.test(Y=cbind(Y,delta), X=X,B=B)$pvalue
    
    ########### proposed tests ###################
    Y<- (Y-mean(Y))/sd(Y)# standardization
    ############## different weight function ##########################
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
  output3<-foreach(m = seq(1,size,by=1),.combine = 'rbind') %dopar% {
    X <- runif(n,-1,1)
    y <- cos(3*pi*X)+  0.5*rnorm(n,0,1)
    C <-  1* rnorm(n,1.6,0.5)- rnorm(n,1,0.1)# censoring variable
    ############## generate censoring ############
    delta <-(y<=C)
    Cr <- 1-mean(delta) # censoring rates
    Y<- pmin(y,C) 
    
    ########### cox PH likelihood ratio test #####
    pvalue.cox <- cph_test(X=as.matrix(X),z=array(Y),d=array(delta))
    ########### KLR ######
    # testing code X=X, z=Y(time), d=delta
    # dim=1
    
    pvalue.KLR <- wild_bootstrap_test_logrank_covariates(
      X=as.matrix(X), z=array(Y),d=array(delta),kernels_x='gau',kernel_z='gau')[[2]]
    
    ########### IPCW ######
    dis.pvalue<- ipcw.dcov.test(Y=cbind(Y,delta), X=X,B=B)$pvalue
    
    ########### proposed tests ###################
    Y<- (Y-mean(Y))/sd(Y)# standardization
    ############## different weight function ##########################
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
  output6<-foreach(m = seq(1,size,by=1),.combine = 'rbind') %dopar% {
    X <- runif(n,-1,1)
    y <- cos(3*pi*X)+  0.5*rnorm(n,0,1)
    C <-  1* rnorm(n,0.75,0.5)- rnorm(n,1,0.1)# censoring variable
    ############## generate censoring ############
    delta <-(y<=C)
    Cr <- 1-mean(delta) # censoring rates
    Y<- pmin(y,C) 
    
    ########### cox PH likelihood ratio test #####
    pvalue.cox <- cph_test(X=as.matrix(X),z=array(Y),d=array(delta))
    ########### KLR ######
    # testing code X=X, z=Y(time), d=delta
    # dim=1
    
    pvalue.KLR <- wild_bootstrap_test_logrank_covariates(
      X=as.matrix(X), z=array(Y),d=array(delta),kernels_x='gau',kernel_z='gau')[[2]]
    
    ########### IPCW ######
    dis.pvalue<- ipcw.dcov.test(Y=cbind(Y,delta), X=X,B=B)$pvalue
    
    ########### proposed tests ###################
    Y<- (Y-mean(Y))/sd(Y)# standardization
    ############## different weight function ##########################
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
  df30[j,] <-c(colMeans(output0[,1:9]<=0.05),mean(output0[,10]))
  df3_30[j,] <-c(colMeans(output3[,1:9]<=0.05),mean(output3[,10]))
  df3_60[j,] <-c(colMeans(output6[,1:9]<=0.05),mean(output6[,10]))
}
end_time<-Sys.time()
times <- end_time-start_time
print(times)
print(df30)
print(df3_30)
print(df3_60)

registerDoParallel(150)
B <- 500
size <- 1000
start_time<-Sys.time()
# cos(8*pi*x)
for (j in 1:length(nset)){
  n <- nset[j]
  ########### for  cluster use ############
  output0<-foreach(m = seq(1,size,by=1),.combine = 'rbind') %dopar% {
    X <- runif(n,-1,1)
    y <- cos(8*pi*X)+ .5*rnorm(n,0,1)
    C <-  1* rnorm(n,1.1,0.5)- rnorm(n,1,0.1)# censoring variable
    C <- 1e7 # CR=0%
    ############## generate censoring ############
    delta <-(y<=C)
    Cr <- 1-mean(delta) # censoring rates
    Y<- pmin(y,C) # observed variable
    
    ########### cox PH likelihood ratio test #####
    pvalue.cox <- cph_test(X=as.matrix(X),z=array(Y),d=array(delta))
    ########### KLR ######
    # testing code X=X, z=Y(time), d=delta
    # dim=1
    pvalue.KLR <- wild_bootstrap_test_logrank_covariates(
      X=as.matrix(X), z=array(Y),d=array(delta),kernels_x='gau',kernel_z='gau')[[2]]
    ########### IPCW ######
    dis.pvalue<- ipcw.dcov.test(Y=cbind(Y,delta), X=X,B=B)$pvalue
    
    ########### proposed tests ###################
    Y<- (Y-mean(Y))/sd(Y)# standardization
    ############## different weight function ##########################
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
  output3<-foreach(m = seq(1,size,by=1),.combine = 'rbind') %dopar% {
    X <- runif(n,-1,1)
    y <- cos(8*pi*X)+  0.5*rnorm(n,0,1)
    C <-  1* rnorm(n,1.6,0.5)- rnorm(n,1,0.1)# censoring variable
    ############## generate censoring ############
    delta <-(y<=C)
    Cr <- 1-mean(delta) # censoring rates
    Y<- pmin(y,C) 
    
    ########### cox PH likelihood ratio test #####
    pvalue.cox <- cph_test(X=as.matrix(X),z=array(Y),d=array(delta))
    ########### KLR ######
    # testing code X=X, z=Y(time), d=delta
    # dim=1
    
    pvalue.KLR <- wild_bootstrap_test_logrank_covariates(
      X=as.matrix(X), z=array(Y),d=array(delta),kernels_x='gau',kernel_z='gau')[[2]]
    
    ########### IPCW ######
    dis.pvalue<- ipcw.dcov.test(Y=cbind(Y,delta), X=X,B=B)$pvalue
    
    ########### proposed tests ###################
    Y<- (Y-mean(Y))/sd(Y)# standardization
    ############## different weight function ##########################
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
  output6<-foreach(m = seq(1,size,by=1),.combine = 'rbind') %dopar% {
    X <- runif(n,-1,1)
    y <- cos(8*pi*X)+  0.5*rnorm(n,0,1)
    C <-  1* rnorm(n,0.75,0.5)- rnorm(n,1,0.1)# censoring variable
    ############## generate censoring ############
    delta <-(y<=C)
    Cr <- 1-mean(delta) # censoring rates
    Y<- pmin(y,C) 
    
    ########### cox PH likelihood ratio test #####
    pvalue.cox <- cph_test(X=as.matrix(X),z=array(Y),d=array(delta))
    ########### KLR ######
    # testing code X=X, z=Y(time), d=delta
    # dim=1
    
    pvalue.KLR <- wild_bootstrap_test_logrank_covariates(
      X=as.matrix(X), z=array(Y),d=array(delta),kernels_x='gau',kernel_z='gau')[[2]]
    
    ########### IPCW ######
    dis.pvalue<- ipcw.dcov.test(Y=cbind(Y,delta), X=X,B=B)$pvalue
    
    ########### proposed tests ###################
    Y<- (Y-mean(Y))/sd(Y)# standardization
    ############## different weight function ##########################
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
  df80[j,] <-c(colMeans(output0[,1:9]<=0.05),mean(output0[,10]))
  df8_30[j,] <-c(colMeans(output3[,1:9]<=0.05),mean(output3[,10]))
  df8_60[j,] <-c(colMeans(output6[,1:9]<=0.05),mean(output6[,10]))
}
end_time<-Sys.time()
times <- end_time-start_time
print(times)
print(df80)
print(df8_30)
print(df8_60)
 
