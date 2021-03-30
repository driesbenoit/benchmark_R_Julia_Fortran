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

binprobbayes <- function(y,X,b0,B0,R,thin){

	n <- nrow(X) 
	k <- ncol(X)

	B <- solve(solve(B0) + t(X)%*%X)
	Bb <- solve(B0)%*%b0

	# Initialize space to save draws
	betadraw <- matrix(NA,nrow=R,ncol=k)

	# Set starting values
	beta <- rep(0, k)
	ystar <- rep(0,n)

	# Start mcmc
	for (i in 1:R){
	
		# Draw next value for ystar
		for (ii in 1:n){
			if (y[ii]==0){
				ystar[ii] <- rtnorm(X[ii,]%*%beta,1,max=0)
			} else {
				ystar[ii] <- rtnorm(X[ii,]%*%beta,1,min=0)
			}
		}

		# Draw new value for beta
		beta <- B%*%(Bb + (t(X)%*%ystar)) + t(chol(B))%*%rnorm(n=k)

		# Save draw
		betadraw[i,] = beta
	}

	return(betadraw)
}
