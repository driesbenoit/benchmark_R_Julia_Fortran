# Written by Dries F. Benoit 
# Faculty of economics and business administration
# Ghent University - BELGIUM

# Returns one draw from the truncated normal distribution

# Algorithm based on:
# Geweke, J. (1991). Efficient Simulation From the Multivariate Normal 
# and Student t-Distributions Subject to Linear Constraints, in Computer 
# Sciences and Statistics Proceedings of the 23d Symposium on the 
# Interface, pp. 571-578.

# Input arguments:
# mu	        -	mean of trunc normal
# sd            -	sd of trunc normal
# positive      -   if true: draw positive value 
#                   if false: draw negative value 

# NOTE: truncation from both above and below is not supported

# Licence: GPLv3 (or higher)

rtnorm <- function(mu, sd, positive){
  if (positive){
    min = 0
    max = Inf
  }else{
    min = -Inf 
    max = 0
  }

  if(is.finite(max)){
    lt <- FALSE
    c <- -(max-mu)/sd
  }else{
    lt=TRUE
    c <- (min-mu)/sd
  }	
  
  if(c <= .45){
    # normal rejection sampling
    repeat{
      x <- rnorm(n=1, mean=0, sd=1)
      if(x>c)break
    }
  } else if (c > .45){
    # exponential rejection sampling
    repeat{
      x <- rexp(n=1, rate=c)
      u <- runif(n=1, min=0, max=1)
      if(u < exp(-.5*(x^2)))break
    }
    x <- x + c
  }

  if(lt){
    return(mu+sd*x)
  }else{
    return(mu-sd*x)
  }
}
