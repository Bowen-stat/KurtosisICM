rm(list = ls())
dig=3
setwd("~\\Dropbox\\ZhouBowen\\Kurtosis-ICM\\sim")
library('stargazer')

TT <- array(0,dim = c(3000,5,3))
for (id in 1:3) {
  TT[,,id]=readRDS(paste('error',id,'-tab1','.rds',sep =''))
}
# TT[,,1]=(TT[,,1])
# TT[,,2]=(TT[,,2]+1.2)
# TT[,,3]=(TT[,,3]-1.2)
Re <- matrix(0,9,5)

ff <- function(obj){
  m<-format(round((apply(obj,2,mean,na.rm=TRUE)),dig),nsmall =dig)
  s<-format(round(apply(obj,2,sd,na.rm=TRUE),dig),nsmall =dig)
  l<-rep('(',length(m));
  r<-rep(')',length(m));
  R1=paste(m,l,s,r,sep='')
  return(R1)
}

for (k in 1:3) {
  for (l in 1:3) {
    Re[l+(k-1)*3,]=ff(TT[1:1000+(l-1)*1000,,k])
  }
}

colnames(Re)<-c('p=100','p=200','p=400','p=800','p=1600');

Re<-cbind(rep(c("proposed","lopes","liu"),3),Re)
stargazer(Re)
