# Written by Dries F. Benoit 
# Faculty of economics and business administration
# Ghent University - BELGIUM

# Draws from the posterior distribution of the Bayesian probit model
# with normal prior on the regression coefficients.

# Algorithm based on:
# Albert and Chib (1993). Bayesian Analysis of Binary and Polychotomous
# Response Data. Journal of the American Statistical Association, 
# 88(422), 669-679.

# Licence GLPv3 (or higher)

using Distributions
using LinearAlgebra

function binprobbayes(y::BitArray{1}, X::Matrix{Float64}, # Data
			b0::Vector{Float64}, B0::Symmetric{Float64,Matrix{Float64}}, # Prior
			R::Int64) # Nbr mcmc draws

  n = size(X,1)
  k = size(X,2)

  # 
  B = Symmetric(inv(inv(B0) + *(X',X)))
  Bb = B0\b0

  # Initialize space to save draws
  #betadraw = zeros(R,k)
  betadraw = Matrix{Float64}(undef, R,k)

  # Set starting values
  beta = zeros(k)
  ystar = zeros(n)
  
  # Start mcmc
  for i = 1:R

    # Draw new value for ystar
    #ystar = map((x,y)->rtnorm(x,1.0,y), X*beta, y)
    Xbeta = X*beta
    for ii = 1:n
      #ystar[ii] = rtnorm(dot(X[ii,:],beta),1.0,y[ii])
      ystar[ii] = rtnorm(Xbeta[ii],1.0,y[ii])
    end

    # Draw new value for beta
    beta = rand(MvNormal(B*(Bb+*(X',ystar)),B))

	# Save draw
    betadraw[i,:] = beta
  end

  return betadraw

end
