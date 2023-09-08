#' Efficient algorithm for Kurtosis estimation of Independent component model
#' @name kurtosisU
#' @param X data matrix of dimension n*p
#' @return A list with components
#' \item{ka}{The esitimation value of kurtosis}
#' \item{tau}{The estimation value of tau in the reference.}
#' @export
kurtosisU = function(X){
  n = nrow(X)
  p = ncol(X)
  if(n<4){stop("sample size must be larger than 4")}
  xxt = X%*%t(X)
  Xti = t(diag(xxt)-2*xxt)+diag(xxt)

  Y1 = sum(Xti^2)/(n*(n-1))

  I1 = (sum(Xti))^2
  I2 = n*(n-1)*Y1
  I3 = sum(Xti%*%Xti)-sum(diag(Xti%*%Xti))
  Y2 = (I1-2*I2-4*I3)/(n*(n-1)*(n-2)*(n-3))

  J1 = sum(xxt^2)-sum(diag(xxt^2))
  J2 = sum(xxt%*%xxt)-sum(diag(xxt%*%xxt))-2*(sum(diag(xxt)*xxt)-sum(diag(diag(xxt)*xxt)))
  J3 = (sum(xxt)-sum(diag(xxt)))^2-2*J1-4*J2
  Y3 = 4*J1/(n*(n-1))-8*J2/(n*(n-1)*(n-2))+4*J3/(n*(n-1)*(n-2)*(n-3))

  Xh = t(X)%*%(rep(1,n)%*%t(rep(1,n))-diag(n))%*%X
  J4 = sum((X^2)%*%t(X^2))-sum(diag((X^2)%*%t(X^2)))
  Xm = (X^3)%*%t(X)
  J5 = sum(diag((Xh)*(t(X)%*%X)))-2*(sum(Xm)-sum(diag(Xm)))
  J6 = sum(diag(Xh^2))-2*J4-4*J5
  Y4 = 4*J4/(n*(n-1))-8*J5/(n*(n-1)*(n-2))+4*J6/(n*(n-1)*(n-2)*(n-3))

  T1 = Y1-Y2-2*Y3
  T2 = Y4/2
  ka = T1/T2
  tau = Y3/Y4

  return(list(ka = ka, tau = tau))
}

