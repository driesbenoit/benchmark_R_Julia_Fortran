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

binprobbayes <- function(y,X,b0,B0,R){

	n <- nrow(X) 
	k <- ncol(X)

	B <- chol2inv(chol(chol2inv(chol(B0)) + t(X)%*%X))
	Bb <- chol2inv(chol(B0))%*%b0

	# Initialize space to save draws
	betadraw <- matrix(NA,nrow=R,ncol=k)

	# Set starting values
	beta <- rep(0, k)
	ystar <- rep(0,n)

	# Start mcmc
	for (i in 1:R){
	
		# Draw next value for ystar
        ystar = mapply(FUN=rtnorm,mu=X%*%beta,sd=1,positive=df$y)

		# Draw new value for beta
		beta <- B%*%(Bb + (t(X)%*%ystar)) + t(chol(B))%*%rnorm(n=k)

		# Save draw
		betadraw[i,] = beta
	}

	return(betadraw)
}
