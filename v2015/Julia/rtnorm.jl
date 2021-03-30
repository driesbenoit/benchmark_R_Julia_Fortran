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
# sigma         -	sd of trunc normal
# min           -       truncation point from below
# max           -       truncation point from above
# NOTE: truncation from both above and below is not supported

# Licence: GPLv3 (or higher)

using Distributions
function rtnorm(mu::Float64, sd::Float64; min::Float64=-Inf, max::Float64=Inf)

  local x

  if isfinite(max)
    lt = false
    c = -(max-mu)/sd
  else
    lt = true
    c = (min-mu)/sd
  end

  if c <= .45 
    # normal rejection sampling
    while true
      x = rand(Normal(0.,1.))
      if x > c
        break
      end
    end
  
  else
    # exponential rejection sampling
    while true
      x = rand(Exponential(1./c)) # check parametrization!
      u = rand()
      if u < exp(-.5*(x^2))
        x = x + c
        break
      end
    end
  end

  if lt
    return mu + sd * x
  else
    return mu - sd * x
  end
  
end
