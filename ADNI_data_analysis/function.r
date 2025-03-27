### our methods######
###functions ####
#Fw <- function(x){pnorm(x, 0,1)}
new_col_mul_fast <-function(Y,X, delta,Fw){
  #### basic functions #######
  Ym <- matIY0(Y)
  n <- length(Y)
  pY<-Fw(Y)
  cIY <-colSums(Ym) 
  rIY <-rowSums(Ym)
  # find the nearest neighbors of X 
  library(RANN)
  nn_index <- nn2(X,query = X,k=2)$nn.idx[,2]
  
  ############## The first term #################
  Qint <- mean(pmin(pY,pY[nn_index]))
   
  ############## The second term #################
  Gint2 <- mean(pY)/n+sum((cIY-1)*pY)/(n^2)+sum((n-rIY)*pY)/(n^2)
   
  ############## The third term #################
  xx<-c(n:1)*c(min(pY),diff(sort(pY)))
  sc <- 1-(1-delta)/cIY 
  newsc <- sc[order(pY)] # by the order of pY
  zz <-c(1,cumprod(newsc)[-n])
  SGint <- sum(xx*zz)/n
  ###### correlation #######
  new_coeff <- (Qint - Gint2)/(SGint - Gint2)
  c(Qint, Gint2,SGint,new_coeff)
} 


permu_pvalue_fast <- function(Y, X, delta, Fw=function(x){pnorm(x, 0,1)}, B=500) {
  n <- length(Y)
  pY<-Fw(Y)
  Ym <- matIY0(Y)
  cIY <- colSums(Ym)
  ############## The first term #################
  # find the nearest neighbors of X 
  library(RANN)
  nn_index <- nn2(X,query = X,k=2)$nn.idx[,2]
  Qint <- mean(pmin(pY,pY[nn_index]))
  
  ############## The second term #################
  Gint2 <- mean(pY)/n+sum((cIY-1)*pY)/(n^2)+sum((n-rowSums(Ym))*pY)/(n^2)
  
  ############## The third term #################
  xx<-c(n:1)*c(min(pY),diff(sort(pY)))
  sc <- 1-(1-delta)/cIY 
  newsc <- sc[order(pY)] # by the order of pY
  zz <-c(1,cumprod(newsc)[-n])
  SGint <- sum(xx*zz)/n
  ###### permutation #######
  DD=SGint-Gint2
  original.coef <- (Qint-Gint2)/DD
  coef.sim <- numeric(B)
  for(m in seq_len(B)) {
    nx.index <- sample(1:n, n)
    mx <- cbind(1:n, nx.index)
    new_nn_index <- mx[order(mx[,2]),][nn_index[nx.index],][,1]
    coef.sim[m] <- (mean(pmin(pY,pY[new_nn_index]))-Gint2)/DD
  }
  pvalue <-rank(-c(original.coef, coef.sim),ties.method ="random")[1]/(B+1)
  ret <- c(pvalue,original.coef)
  return(ret)
}
## for large dynamic B
permu_pvalue_fast1 <- function(Y, X, delta, Fw=function(x){pnorm(x, 0,1)}, B.init=500,B.max=1e7,B.iter=3)
{
	B=B.init
	while(1)
	{
		ret=permu_pvalue_fast(Y, X, delta, Fw=Fw,B=B)
		SSS=ret[1]*B
		if((SSS>100)|(B>B.max))return(ret)
		if(SSS<=100)B=B*B.iter
		if(B>100000)print(c(ret,B))
	}
}

ipcw.dcov.test1 <- function(Y,X,B.init=500,B.max=1e7,B.iter=3)
{
	B=B.init
	while(1)
	{
		ret=ipcw.dcov.test (Y, X, B=B)$pvalue
		SSS=ret*B
		if((SSS>50)|(B>B.max))return(ret)
		if(SSS<=50)B=B*B.iter
		if(B>100000)print(c(ret,B))
	}
}

wild_bootstrap_LR1<-function(X,Y,d,ker,ker1,kernel_parameters_x,kernel_parameters_z,seed=as.integer(1),B.init=as.integer(500),
B.max=1e7,B.iter=3)
{
	B=as.integer(B.init)
	while(1)
	{
		ret=wild_bootstrap_test_logrank_covariates(X,Y,d,ker,ker1,kernel_parameters_x,kernel_parameters_z,seed,B)[[2]]
		SSS=ret*B
		if((SSS>50)|(B>B.max))return(ret)
		if(SSS<=50)B=as.integer(B*B.iter)
		if(B>100000)print(c(ret,B))
	}
} 