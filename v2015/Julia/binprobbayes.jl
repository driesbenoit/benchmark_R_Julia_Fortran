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
function binprobbayes(y::Array{Int64,1}, X::Array{Float64,2}, # Data
			b0::Array{Float64,1}, B0::Array{Float64,2}, # Prior
			R::Int64, thin::Int64) # Mcmc

  n = size(X,1)
  k = size(X,2)

  #B = invvarcov(invvarcov(B0) + *(X',X))
  #Bb = *(invvarcov(B0),b0)
  B = inv(inv(B0) + *(X',X))
  Bb = *(inv(B0),b0)

  # Initialize space to save draws
  #betadraw = reshape(repeat([NaN],inner=[R*k]),R,k)
  betadraw = Array(Float64,R,k)

  # Set starting values
  beta = repeat([0.],inner=[k])
  ystar = repeat([0.],inner=[n])
  
  # Start mcmc
  for i = 1:R

        # Draw new value for ystar
        for ii = 1:n
                if y[ii]==0
                        ystar[ii] = rtnorm(*(X[ii,:],beta)[1],1.,max=0.)
                else
                        ystar[ii] = rtnorm(*(X[ii,:],beta)[1],1.,min=0.)
                end
        end

        # Draw new value for beta
        beta = rand(MvNormal(B*(Bb+*(X',ystar)),B))

	# Save draw
        betadraw[i,:] = beta

  end

  return betadraw

end
