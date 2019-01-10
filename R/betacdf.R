#' @title distribution function for the Beta distribution
#' @description a function to compute a Monte Carlo estimate of the Beta cdf Using antithetic variables method
#' @param a,b non-negative parameters of the Beta distribution
#' @param x the value of quantile
#' @return the value of the Beta cdf at x
#' @examples
#' \dontrun{
#' betacdf(3,3,0.3)
#' }
#' @export
betacdf<-function(a,b,x){
  m=1000000
  u=runif(m/2,min=0,max=x)
  u1=1-u
  g=1/beta(a,b)*u^(a-1)*(1-u)^(b-1)*x
  g1=1/beta(a,b)*u1^(a-1)*(1-u1)^(b-1)*x
  theta=(mean(g)+mean(g1))/2
  return(theta)
}
