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

# Licence: GPLv3 (or higher)

using Distributions
function rtnorm(mu::Float64, sd::Float64, positive::Bool)

  if positive
    min=0.0
    max=Inf
  else
    min=0.0
    max=0.0
  end

  local x::Float64

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
      x = rand(Normal(0,1))
      if x > c
        break
      end
    end
  
  else
    # exponential rejection sampling
    while true
      x = rand(Exponential(1/c)) # check parametrization!
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
