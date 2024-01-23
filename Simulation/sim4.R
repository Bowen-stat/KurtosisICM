rm(list = ls())
library(kurtICM)
library(tictoc)
setwd("~\\Dropbox\\ZhouBowen\\Kurtosis-ICM\\sim")

sqrtm <- function(sig){
  eig = eigen(sig)
  lamb = eig$vectors%*%diag(sqrt(eig$values))%*%t(eig$vectors)
  return(lamb)
}
rho=0.5
p=200
sigma <- toeplitz(rho^(1:p-1))
lamb <- sqrtm(sigma)

ddelta <- seq(-1.5,1.5,0.1)
ggamma <- (6+ddelta+sqrt((6+ddelta)*(2+ddelta)))/2/(6+ddelta)




ff <- function(n,gamma){
  Z = matrix(rbinom(n*p,size = 1,prob = gamma),nrow = n)
  Z[Z==1]=-sqrt((1-gamma)/gamma)
  Z[Z==0]=sqrt(gamma/(1-gamma))
  X = Z%*%lamb
  obj <- Ktest(X)
  return(obj$p.value)
}
id=3
nn = c(100,200,400)
n=nn[id]
m=1e3
tic()
re = mapply(ff, n,rep(ggamma,m))
toc()

saveRDS(re,file=paste('powern100-',id,'-tab4','.rds',sep =''))



