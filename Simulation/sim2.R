rm(list = ls())
set.seed(2024)
library(kurtosisICM)
setwd("~\\Dropbox\\ZhouBowen\\Kurtosis-ICM\\sim")
n=100
p=200
id=1 #select distribution

source("ffun.R")


sqrtm <- function(sig){
  eig = eigen(sig)
  lamb = eig$vectors%*%diag(sqrt(eig$values))%*%t(eig$vectors)
  return(lamb)
}

gg <- function(id,lamb){
  if(id==1){Z = matrix(rnorm(n*p),nrow = n)}else if(id==2){
    Z = matrix(runif(n*p,min=-sqrt(3),max = sqrt(3)),nrow = n)}else if(id==3){
      Z = matrix(rt(n*p,df=9),nrow = n); Z = Z/sqrt(9/7)}
  X = Z%*%lamb
  obj <- kurtICM(X)
  return(obj$ke)
}

ff <- function(rho,id){
  sigma <- toeplitz(rho^(1:p-1))
  lamb <- sqrtm(sigma)
  re <- replicate(1e3,gg(id,lamb))
  return(re)
}
rrho <- seq(0,0.7,0.1)

for (id in 1:3) {
if(id==1){d0=0}else if(id==2){d0=-1.2}else if(id==3){d0=1.2}
  Re <- mapply(ff, rrho, id)
data = as.data.frame(cbind(as.vector(Re),rep(rrho,each=1e3)))
names(data) = c("ka","rho")
CairoPDF(paste('fig',id,'-','boxplot-sim2','.pdf',sep=''))
boxplot(ka~rho,data = data,xlab = expression(rho), ylab = NULL, ylim = c(d0-3,d0+4))
dev.off()
}





