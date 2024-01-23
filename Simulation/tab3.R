rm(list = ls())
dig=3
setwd("~\\Dropbox\\ZhouBowen\\Kurtosis-ICM\\sim")
library('stargazer')
n=100
TT <- array(0,dim = c(1000,15,3))
for (id in 1:3) {
  TT[,,id]=readRDS(paste('testmid-',id,'-tab3','.rds',sep =''))
}

Re = matrix(0,9,5)

ff<-function(obj){
  z = apply(obj<=0.05,2,mean)
  return(matrix(z,nrow = 3))
}

for (id in 1:3) {
  if(id==1){d0=0}else if(id==2){d0=-1.2}else if(id==3){d0=1.2}
  Re[1:3+(id-1)*3,]=ff(TT[,,id])
}

colnames(Re)<-c('p=100','p=200','p=400','p=800','p=1600');
Re<-cbind(rep(c("0","0.5","0.9"),3),Re)

stargazer(Re)
