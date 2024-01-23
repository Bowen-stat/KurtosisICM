rm(list = ls())
setwd("~\\Dropbox\\ZhouBowen\\Kurtosis-ICM\\sim")
library('stargazer')
library(Cairo)
TT <- array(0,dim = c(2,31000,3))
for (id in 1:3) {
  TT[,,id]=readRDS(paste('powern100-',id,'-tab4','.rds',sep =''))
}
ddelta <- seq(-1.5,1.5,0.1)
nn = c(100,200,400)
ff<-function(obj){
  ke = matrix(obj[1,],nrow = length(ddelta))
  taue = matrix(obj[2,],nrow = length(ddelta))
  power = apply(abs((ke)/((taue*2)*sqrt(2/n)))>qnorm(0.975), 1, mean)
  return(power)
}
power = matrix(0,3,length(ddelta))
for (id in 1:3) {
  n = nn[id]
  power[id,] = ff(TT[,,id])
}

miny = 0
maxy = 1
CairoPDF(paste('fig','-','power','-sim4','.pdf',sep=''))
plot(ddelta,power[1,], pch = 16, lty = 1, col = "blue", type = "o", ylim = c(miny,maxy),xlab = expression(Delta), ylab = "Power",lwd=2)
lines(ddelta,power[2,], pch = 1, lty = 1, col = "red", type = "o", ylim = c(miny,maxy),lwd=2)
lines(ddelta,power[3,], pch = 8, lty = 1, col = "green", type = "o", ylim = c(miny,maxy),lwd=2)
legend('bottomright',legend =c('n=100','n=200','n=400'),pch=c(16,1,8),lty=c(1,1,1),col = c("blue","red","green"))
dev.off()
