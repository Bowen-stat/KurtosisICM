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
  Xti2 = Xti%*%Xti
  I3 = sum(Xti2)-sum(diag(Xti2))
  Y2 = (I1-2*I2-4*I3)/(n*(n-1)*(n-2)*(n-3))

  J1 = sum(xxt^2)-sum(diag(xxt^2))
  xxt2 = xxt%*%xxt
  xxt3 = diag(xxt)*xxt
  J2 = sum(xxt2)-sum(diag(xxt2))-2*(sum(xxt3)-sum(diag(xxt3)))
  J3 = (sum(xxt)-sum(diag(xxt)))^2-2*J1-4*J2
  Y3 = 4*J1/(n*(n-1))-8*J2/(n*(n-1)*(n-2))+4*J3/(n*(n-1)*(n-2)*(n-3))

  Xh = t(X)%*%(rep(1,n)%*%t(rep(1,n))-diag(n))%*%X
  X22 = (X^2)%*%t(X^2)
  L1 = sum(X22)-sum(diag(X22))
  Xm = (X^3)%*%t(X)
  L2 = sum(diag((Xh)*(t(X)%*%X)))-2*(sum(Xm)-sum(diag(Xm)))
  L3 = sum(diag(Xh^2))-2*L1-4*L2
  Y4 = 4*L1/(n*(n-1))-8*L2/(n*(n-1)*(n-2))+4*L3/(n*(n-1)*(n-2)*(n-3))

  T1 = Y1-Y2-2*Y3
  T2 = Y4/2
  ka = T1/T2
  tau = Y3/Y4

  return(list(ka = ka, tau = tau))
}

