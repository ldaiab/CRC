################################## y = x1^2+epsilon ###################################################
rm(list = ls())
library(MASS)
library(RANN)
library(cubature)
library(matrixStats)
library(doParallel)
library(foreach)
library(RcppEigen)
library(Rcpp)
library(pracma)
sourceCpp(file="./fun_indep.cpp")


### true correlation coefficients 
true_col <-  function(wt,thres = Inf){
  # the first term 
  Qintegrand <- function(x,thres){
    len <- length(x)
    (1-pchisq(x[1]-(x[2])^2,2))^2*
      dnorm(x[2],0,1)*wt(x[1])*(x[1]>-thres)*(x[1]<thres)
  }
  len <-length(beta)
  Qtint <- hcubature(Qintegrand, c(0,-Inf), c(Inf,Inf), tol=1e-3,thres = thres)$integral
  # the second term 
  G2integrand <- function(x){(1-pchisq(x,3))^2*wt(x)}
  G2tint <- hcubature(G2integrand, 0, thres, tol=1e-6)$integral
  # the third term  
  Gintegrand <- function(x){(1-pchisq(x,3))*wt(x)}
  Gtint <- hcubature(Gintegrand, 0, thres, tol=1e-6)$integral
  # the true coefficient 
  true_col <- (Qtint-G2tint)/(Gtint-G2tint)
  return(true_col)
}

# Fw is the cumulative w(t)
# w(t) = chisq(3)
wt <- function(x){dchisq(x,3)} #fy(t)
Fw <- function(x){pchisq(x,3)} #fy(t)
true_value_col <- true_col(wt)
true_value_col
# w(t) = exp(1)
wt <- function(x){dexp(x,1)}
Fw <- function(x){pexp(x,1)}
true_value_col <- true_col(wt)
true_value_col
# w(t) = chisq(4)
wt <- function(x){dchisq(x,4)}
Fw <- function(x){pchisq(x,4)}
true_value_col <- true_col(wt)
true_value_col
# w(t) = U(-1,1)
wt <- function(x){dunif(log(x), -1,1)/x}
Fw <- function(x){punif(log(x), -1,1)}
true_value_col <- true_col(wt)  
true_value_col
# w(t) = U(-5,5)
wt <- function(x){dunif(log(x), -5,5)/x}
Fw <- function(x){punif(log(x), -5,5)}
true_value_col <- true_col(wt)  
true_value_col
# w(t) = N(0,1)
wt <- function(x) {dnorm(log(x),0,1)/x}
Fw <- function(x) {pnorm(log(x),0,1)}
true_value_col <-  true_col(wt) 
true_value_col
# w(t) = N(-1,50)
wt <- function(x) {dnorm(log(x),-1,5)/x}
Fw <- function(x) {pnorm(log(x),-1,5)}
true_value_col <-  true_col(wt) 
true_value_col
# w(t) = t(2)
wt <- function(x) {dt(log(x),2)/x}
Fw <- function(x) {pt(log(x),2)}
true_value_col <-  true_col(wt) 
true_value_col
# w(t) = t(3)
wt <- function(x) {dt(log(x),3)/x}
Fw <- function(x) {pt(log(x),3)}
true_value_col <-  true_col(wt) 
true_value_col

n <- 800
size<-1000
registerDoParallel(100)
start_time<-Sys.time()
output<-foreach(m = seq(1,size,by=1),.combine = 'rbind') %dopar% {
  # DGP2
  X1 <- rnorm(n,0,1)
  gx <- X1^2 
  y <- gx + rchisq(n,2)
  X <- X1
  #C <-  1.1*rgamma(n,5,scale = 1) + 0.5*runif(n,0,1) # censoring variable, chiq1
  C <-  1.6*rgamma(n,5,scale = 1)+0.5*runif(n,0,1) # censoring variable, chiq2
  delta <-(y<=C)
  Cr <- 1-mean(delta) # censoring rates
  Y<- pmin(y,C) # observed variable
  new_correlation <- new_col_mul_fast(Y, X,delta,Fw)
  res <- c(new_correlation-true_value_col,Cr)
  return(res)
}
end_time<-Sys.time()
times <- end_time-start_time
print(times)
round(colMeans(output),7)
print(colSds(output))

