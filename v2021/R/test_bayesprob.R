# Set working directory
setwd("/.")

# Source the algorithms
source("binprobbayes.R")
source("rtnorm.R")

# Read in the data
df <- read.csv("probitdata.csv")

# Execute algorithm
tic <- proc.time()
out <- binprobbayes(df$y,as.matrix(df[,2:4]),c(0,0,0),diag(3)*1000,5000)
proc.time()-tic

# Check Bayes estimate
colMeans(out)

# Plot trace plots
matplot(out,typ="l",lty=1)
