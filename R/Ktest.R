#' Hypothesis test for the Kurtosis  of Independent component models
#' @name Ktest
#' @param X data matrix of dimension n*p
#' @param k0 the value of kurtosis under the null hypothesis
#' @param alpha the confidence level
#' @param alternative a character string specifying the alternative hypothesis, must be one of '"two.sided"' (default), '"greater"' or '"less"'.
#' @return A list with class \code{htest} containing the following components:
#' \item{statistic}{the list containing kurtosis and related estimators and its transformation.}
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
  ke = kurtosisU(X)
  ka = ke[[1]]
  tau = ke[[2]]
  varh = sqrt(2/n)*(k0+2*tau)
  ztest = (ka-k0)/varh
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
  Re <- list(statistic = c(kurtosis = ka, tau = tau, z = ztest), p.value = pval,
             alternative = alt, data.name = DNAME, method = "test on kurtosis of ICM")
  class(Re) <- "htest"
  return(Re)
}
