kurtlopes <- function(X){
  n = nrow(X)
  p = ncol(X)
  samsig = cov(X)
  X = scale(X,scale = F,center = T)
  
  taun = sum(samsig^2)-((sum(diag(samsig)))^2)/n
  xf = apply(X^2, 1, sum)
  vn = var(xf)
  w = sum((apply(X^2, 2, mean))^2)
  ka = (vn-2*taun)/w
  return(ka)
}
kurtliu <- function(X){
  n = nrow(X)
  p = ncol(X)
  Xs = scale(X,scale = F,center = T)
  XXt = Xs%*%t(Xs)
  Y1 = sum((diag(XXt))^2)/(n-1)
  Y2 = (sum(diag(XXt)))^2/n/(n-1)
  Y3 = 2*(sum(XXt^2)-(n-1)*Y1)/n/(n-1)
  Y4 = sum((apply(Xs^2, 2, mean))^2)
  ka = (Y1-Y2-Y3)/Y4
  return(ka)
}