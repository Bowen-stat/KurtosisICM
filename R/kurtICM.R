#' Efficient algorithm for Kurtosis estimation of Independent component model
#' @name kurtICM
#' @param X data matrix of dimension n*p
#' @return A list with components
#' \item{ke}{The esitimation value of kurtosis}
#' \item{sde}{The estimation value of standard deviation of ka}
#' @export

kurtICM <- function(X){
  n = nrow(X)
  p = ncol(X)
  if(n<4){stop("sample size must be larger than 4")}
  Xc = scale(X,scale=F)
  sig = cov(X)
  G = Xc%*%t(Xc)
  D = t(diag(G)-2*G)+diag(G)
  b = rowSums(D)

  Y1 = 2*((n-1)*(n-4)*sum(D^2)-sum(b)^2+4*sum(b^2))/n/(n-1)/(n-2)/(n-3)
  Y2 = 4*(n-1)*(n^2-3*n+2)*sum(sig^2)/n/(n-2)/(n-3)+4*(n-1)*(sum(diag(sig)))^2/n/(n-2)/(n-3)-4*sum(diag(G)^2)/(n-2)/(n-3)
  Y3 = 4*(n-1)*(n^2-3*n+3)*sum(diag(sig)^2)/n/(n-2)/(n-3)-4*sum(Xc^4)/(n-2)/(n-3)

  ke = (Y1-4*Y2)/Y3
  sde = sqrt(2)*(Y1-2*Y2)/Y3/sqrt(n)

  return(list(ke = ke, sde = sde))
}
