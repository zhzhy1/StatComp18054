#' @title generate samples from a Rayleigh distribution
#' @description a function to generate samples from a Rayleigh distribution
#' @param n number of observations
#' @param sigma non-negative parameter of the Rayleigh distribution
#' @param antithetic logical;if TRUE the antithetic variables method will be used
#' @return the random samples of size n
#' @examples
#' \dontrun{
#' Rayleighsample(1000,2)
#' sample=Rayleighsample(1000,2,antithetic=FALSE)
#' plot(sample)
#' }
#' @export
Rayleighsample<-function(n,sigma,antithetic=TRUE){
  u=runif(n/2)
  if(!antithetic)
  {v=runif(n/2)}else
  {v=1-u}
  u=c(u,v)
  x=sqrt(-2*sigma^2*log(1-u))
  return(x)
}

