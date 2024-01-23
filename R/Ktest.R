#' Hypothesis test for the Kurtosis  of Independent component models
#' @name Ktest
#' @param X data matrix of dimension n*p
#' @param k0 the value of kurtosis under the null hypothesis
#' @param alternative a character string specifying the alternative hypothesis, must be one of '"two.sided"' (default), '"greater"' or '"less"'.
#' @return A list with class \item{htest}
#' containing the following components:
#' \item{statistic}{the list containing kurtosis and its standard deviation and the test statistic.}
#' \item{p.value}{the p-value for the test}
#' \item{alternative}{a character string describing the alternative hypothesis}
#' \item{data.name}{name of the data argument}
#' @export
Ktest = function(X, k0, alternative=c("two.sided","less","greater")){
  DNAME <- deparse(substitute(X))
  s <- match.arg(alternative)
  alter <- switch(s, two.sided=0, less=1, greater=2)
  n = nrow(X)
  if(n<4){stop("sample size must be larger than 4")}
  n = nrow(X)
  p = ncol(X)
  Xc = scale(X,scale=F)
  sig = cov(X)
  G = Xc%*%t(Xc)
  D = t(diag(G)-2*G)+diag(G)
  b = rowSums(D)

  Y1 = 2*((n-1)*(n-4)*sum(D^2)-sum(b)^2+4*sum(b^2))/n/(n-1)/(n-2)/(n-3)
  Y2 = 4*(n-1)*(n^2-3*n+2)*sum(sig^2)/n/(n-2)/(n-3)+4*(n-1)*(sum(diag(sig)))^2/n/(n-2)/(n-3)-4*sum(diag(G)^2)/(n-2)/(n-3)
  Y3 = 4*(n-1)*(n^2-3*n+3)*sum(diag(sig)^2)/n/(n-2)/(n-3)-4*sum(Xc^4)/(n-2)/(n-3)

  ke = (Y1-4*Y2)/Y3
  taue = Y2/Y3

  ztest = (ke-k0)/((taue*2+k0)*sqrt(2/n))
  pval = pnorm(ztest,lower.tail = F)
  if(alter == 0){
    pval <- 2*pval
    if(pval>1){pval <- 2-pval}
    alt <- paste("kurtosis is not equal to", k0)
  }else if(altet == 1){
    alt <- paste("kurtosis is larger than", k0)
  }else{
    pval <- 1-pval
    alt <- paste("kurtosis is less than", k0)
  }
  Re <- list(statistic = c(kurtosis = ke, z = ztest), p.value = pval,
             alternative = alt, data.name = DNAME, method = "test on kurtosis of ICM")
  class(Re) <- "htest"
  return(Re)
}
