rm(list = ls())
library(kurtICM)
setwd("~\\Dropbox\\ZhouBowen\\Kurtosis-ICM\\sim")
n=100
id=3 #select distribution
rho = 0.5 #select covariance
source("ffun.R")

sqrtm <- function(sig){
  eig = eigen(sig)
  lamb = eig$vectors%*%diag(sqrt(eig$values))%*%t(eig$vectors)
  return(lamb)
}

gg <- function(id,lamb){
  p <- nrow(lamb)
  if(id==1){Z = matrix(rnorm(n*p),nrow = n)}else if(id==2){
    Z = matrix(runif(n*p,min=-sqrt(3),max = sqrt(3)),nrow = n)}else if(id==3){
      Z = matrix(rt(n*p,df=9),nrow = n); Z = Z/sqrt(9/7)}
  X = Z%*%lamb
  obj1 <- kurtICM(X)
  ke2 <- kurtlopes(X)
  ke3 <- kurtliu(X)
  return(c(obj1$ke, ke2, ke3))
}

ff <- function(p,id){
  sigma <- toeplitz(rho^(1:p-1))
  lamb <- sqrtm(sigma)
  re <- replicate(1e3,gg(id,lamb))
  return(as.vector(t(re)))
}
pp<-c(100,200,400,800,1600)
tic()
Re <- mapply(ff, pp, id)
toc()

saveRDS(Re,file=paste('error',id,'-tab1','.rds',sep =''))

