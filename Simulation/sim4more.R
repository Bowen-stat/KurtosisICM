rm(list = ls())
library(kurtICM)
library(tictoc)
library(Cairo)
setwd("~\\Dropbox\\ZhouBowen\\Kurtosis-ICM\\sim")

sqrtm <- function(sig){
  eig = eigen(sig)
  lamb = eig$vectors%*%diag(sqrt(eig$values))%*%t(eig$vectors)
  return(lamb)
}
rho=0.5
n=200
ddelta <- seq(-1.5,1.5,0.1)
ggamma <- (6+ddelta+sqrt((6+ddelta)*(2+ddelta)))/2/(6+ddelta)


gg<-function(lamb,gamma){
  p=nrow(lamb)
  Z = matrix(rbinom(n*p,size = 1,prob = gamma),nrow = n)
  Z[Z==1]=-sqrt((1-gamma)/gamma)
  Z[Z==0]=sqrt(gamma/(1-gamma))
  X = Z%*%lamb
  obj <- Ktest(X,0)
  return(obj$p.value)
}

ff<-function(p,gamma){
  sigma <- toeplitz(rho^(1:p-1))
  lamb <- sqrtm(sigma)
  re <- replicate(1e3,gg(lamb,gamma))
  return(re)
}


pp = 800#c(100,200,400)

tic()
Re <- mapply(ff,rep(pp, each = length(ggamma)), rep(ggamma,length(pp)))
toc()


saveRDS(Re,file=paste('powern2',n,'-tab4','.rds',sep =''))

TT = apply(`powerp--tab4`<=0.05, 2, mean)
powers = matrix(TT, nrow = 31)
miny = 0
maxy = 1
CairoPDF(paste('fig','-','power','-sim4more','.pdf',sep=''))
plot(ddelta,powers[,1], pch = 16, lty = 1, col = "blue", type = "o", ylim = c(miny,maxy),xlab = expression(Delta), ylab = "Power",lwd=2)
lines(ddelta,powers[,2], pch = 1, lty = 1, col = "red", type = "o", ylim = c(miny,maxy),lwd=2)
lines(ddelta,powers[,3], pch = 8, lty = 1, col = "green", type = "o", ylim = c(miny,maxy),lwd=2)
legend('bottomright',legend =c('p=100','p=200','p=400'),pch=c(16,1,8),lty=c(1,1,1),col = c("blue","red","green"),lwd=c(2,2,2))
dev.off()
