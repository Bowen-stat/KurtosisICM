rm(list = ls())
library(kurtosisICM)
library(tictoc)
setwd("~\\Dropbox\\ZhouBowen\\Kurtosis-ICM\\sim")
n=100
id=1 #select distribution


sqrtm <- function(sig){
  eig = eigen(sig)
  lamb = eig$vectors%*%diag(sqrt(eig$values))%*%t(eig$vectors)
  return(lamb)
}
gg <- function(id,lamb){
  p <- nrow(lamb)
  if(id==1){Z = matrix(rnorm(n*p),nrow = n);k0=0}else if(id==2){
    Z = matrix(runif(n*p,min=-sqrt(3),max = sqrt(3)),nrow = n);k0=-1.2}else if(id==3){
      Z = matrix(rt(n*p,df=9),nrow = n); Z = Z/sqrt(9/7);k0=1.2}
  X = Z%*%lamb
  obj <- Ktest(X,k0)
  return(obj$p.value)
}
m=1e3
ff <- function(id,rho,p){
  sigma <- toeplitz(rho^(1:p-1))
  lamb <- sqrtm(sigma)
  re <- replicate(m,gg(id,lamb))
  return(re)
}
rrho = c(0.2,0.5,0.8)
pp = c(100,200,400,800,1600)

for (id in 1:3) {
  tic()
  Re <- mapply(ff, id, rep(rrho,length(pp)), rep(pp, each = length(rrho)))
  toc()
  saveRDS(Re,file=paste('testmid-',id,'-tab3','.rds',sep =''))
}







