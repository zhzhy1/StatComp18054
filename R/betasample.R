#' @title generate samples from the Beta distribution
#' @description a function to generate random samples from the Beta distribution by the acceptance-rejection method
#' @param n number of observations
#' @param a,b non-negative parameter of the Beta distribution
#' @return the random samples of size n
#' @examples
#' \dontrun{
#' sample=betasample(1000,3,2)
#' plot(sample)
#' }
#' @export
betasample<-function(n,a,b){
  m=(1-a)/(2-a-b)
  max=m^(a-1)*((1-m)^(b-1))/beta(a,b)#maximum of f(x)
  c<-max+3
  j<-k<-0
  y<-numeric(n)
  while (k < n) {
    u <- runif(1)
    j <- j + 1
    x <- runif(1)
    rho<-x^(a-1)*((1-x)^(b-1))/beta(a,b)/c
    if (rho> u) {
      k <- k + 1
      y[k] <- x
    }
  }
  return(y)
}
