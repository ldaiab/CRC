#################### DGP1: linear model, consistency ################################
#################### linear models: Y = X*beta+epsilon ######################
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
sourceCpp(file="C:/Users/fun_indep.cpp")

source("./functionLL.R")

rho <- 0
Sigma <- matrix(rho,3,3)
diag(Sigma) <- 1
beta <- matrix(c(2,2,4),3,1) 

mu <- rep(0, 3)
sigma<- 5

# some density
dmvnorm <- function(t,mu,Sigma){
  d<- length(t)
  result <- (2*pi)^{-d/2}*abs(det(Sigma))^{-1/2}*exp(-1/2*t(t-mu)%*%solve(Sigma)%*%(t-mu))
  c(result)
}
# true correlation coefficients  
true_col <-  function(wt){
  
  # the first term  
  Qintegrand <- function(x){
    len <- length(x)
    (1-pnorm(x[1],c(x[2:len]%*%beta),sigma))^2*
      dmvnorm(x[2:len],mu,Sigma)*wt(x[1])
  }
  len <-length(beta)
  Qtint <- hcubature(Qintegrand, rep(-Inf, len+1), rep(Inf, len+1), tol=1e-3)$integral
  #  the second term  
  Sigmay <- c(t(beta)%*%Sigma%*%beta+sigma^2)
  muy <- c(mu%*%beta)
  G2integrand <- function(x){(1-pnorm(x,muy,sqrt(Sigmay)))^2*wt(x)}
  G2tint <- hcubature(G2integrand, -Inf, Inf, tol=1e-3)$integral
  #  the third term  
  Gintegrand <- function(x){(1-pnorm(x,muy,sqrt(Sigmay)))*wt(x)}
  Gtint <- hcubature(Gintegrand, -Inf, Inf, tol=1e-3)$integral
  
  # the true coefficient  
  true_col <- (Qtint-G2tint)/(Gtint-G2tint)
  return(true_col)
}
# w(t)=f_y(t)
Sigmay <- c(t(beta)%*%Sigma%*%beta+sigma^2)
muy <- c(mu%*%beta)
wt <- function(x){
  dnorm(x,muy,sqrt(Sigmay))
}
Qintegrand <- function(x){
  len <- length(x)
  (1-pnorm(x[1],c(x[2:len]%*%beta),sigma))^2*
    dmvnorm(x[2:len],mu,Sigma)*wt(x[1])
}
len <-length(beta)
Qtint <- cuhre(f = Qintegrand,
               lowerLimit = rep(-Inf, len+1),
               upperLimit = rep(Inf, len+1),
               relTol = 1e-12, absTol= 1e-12)$integral
true_value_col<- c(Qtint,1/3,1/2,6*Qtint-2)
Fw<- function(x){pnorm(x,muy,sqrt(Sigmay))} 

# Fw is the cumulative w(t)
# w(t) = exp(1)
wt <- function(x){dexp(exp(x), 1)*exp(x)}
Fw<- function(x){pexp(exp(x), 1)} 
true_value_col<- true_col(wt) 
# w(t) = chisq(3)
wt <- function(x){dchisq(exp(x), 3)*exp(x)}
Fw <- function(x){pchisq(exp(x), 3)}
true_value_col <- true_col(wt)  
true_value_col
# w(t) = Unif(-50,50) 
wt <- function(x){dunif(x, -50,50)}
Fw <- function(x){punif(x, -50,50)}
true_value_col <- true_col(wt)  
true_value_col
# w(t) = N(0,1)
wt <- function(x) {dnorm(x,0,1)}
Fw <- function(x) {pnorm(x,0,1)}
true_value_col <-  true_col(wt) 
true_value_col
# w(t) = N(-1,50)
wt <- function(x) {dnorm(x,-1,50)}
Fw <- function(x) {pnorm(x,-1,50)}
true_value_col <-  true_col(wt) 
true_value_col
# w(t) = t(2)
wt <- function(x) {dt(x,2)} 
Fw <- function(x) {pt(x,2)} 
true_value_col <-  true_col(wt) 
true_value_col
# w(t) = t(3)
wt <- function(x) {dt(x,3)}  
Fw <- function(x) {pt(x,3)}  
true_value_col <-  true_col(wt)  
true_value_col
##################### DGP1: linear model #########################
n <- 800
size<-1000
########### for cluster use ############
registerDoParallel(100)
start_time<-Sys.time()
output<-foreach(m = seq(1,size,by=1),.combine = 'rbind') %dopar% {
  # generating censored data #
  X <- mvrnorm(n, mu, Sigma)
  y <- c(X%*%beta) + rnorm(n,0,sigma)
  C <-  rnorm(n,3,1)-rnorm(n,1,1)+1.3*runif(n,0,3) # censoring variable=0.3
  delta <-(y<=C)
  Cr <- 1-mean(delta)  
  Y<- pmin(y,C) 
  # our CRC correlation
  new_correlation <- new_col_mul_fast(Y, X,delta,Fw)
  res <- c(new_correlation-true_value_col,Cr)
  return(res)
}
end_time<-Sys.time()
times <- end_time-start_time
print(times)
print(colMeans(output))
print(colSds(output))

##### normal plots #################
n <- 1000
size<-10000
########### for cluster use ############
registerDoParallel(100)
start_time<-Sys.time()
output<-foreach(m = seq(1,size,by=1),.combine = 'rbind') %dopar% {
  # generating censored data #
  X <- mvrnorm(n, mu, Sigma)
  y <- c(X%*%beta) + rnorm(n,0,sigma)
  C <-  rnorm(n,3,1)-rnorm(n,1,1)+1.3*runif(n,0,3) # censoring variable=0.3
  delta <-(y<=C)
  Cr <- 1-mean(delta)  
  Y<- pmin(y,C) 
  # our coefficients
  new_correlation <- new_col_mul_fast(Y, X,delta,Fw)
  res <- c(new_correlation-true_value_col,Cr)
  return(res)
}
end_time<-Sys.time()
times <- end_time-start_time
print(times)
print(colMeans(output))
print(colSds(output))